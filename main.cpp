#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "bam.hpp"
#include <sstream>
#include <boost/algorithm/string/find.hpp>
#include "sam.hpp"
#include "fasta.hpp"
#include "interval.hpp"
#include "read_filter.hpp"
#include "interval.hpp"
#include "assembler.hpp"
#include "read_clipper.hpp"
#include "smithwaterman.hpp"
#include "pairhmm.hpp"
#include <random>
#include "genetyper.hpp"

using namespace biovoltron::format;
using namespace hc;

static const std::string ref_prefix = "ref/";
static const std::string output_prefix = "output/";
static const auto bam_path = "input/hg19.bam";
static const auto bai_path = "input/hg19.bam.bai";

void call_region(std::vector<SAMRecord>& reads,
                 std::string_view ref,
                 const Interval& padded_region,
                 const Interval& origin_region,
                 std::size_t max_reads_on_assembly_region,
                 std::ostream& os);

void do_work(const Interval& process_region,
             std::size_t assembly_region_size = 245,
             std::size_t assembly_region_padding_size = 85,
             std::size_t max_reads_on_assembly_region = 200)
{
    const auto chromosome = process_region.contig;
    const auto fasta_path = ref_prefix + chromosome + ".fa";
    const auto vcf_path   = output_prefix + process_region.to_string() + ".vcf";

    auto process_padded_region = process_region;
    process_padded_region.end += assembly_region_padding_size;
    if (process_region.begin > assembly_region_padding_size) process_padded_region.begin -= assembly_region_padding_size;
    else process_padded_region.begin = 0;

    auto process_region_ref = std::string{};
    {
        auto fasta_record = biovoltron::format::Fasta{};
        std::ifstream ifs(fasta_path);
        ifs >> fasta_record;
        process_region_ref = fasta_record.seq.substr(process_padded_region.begin, process_padded_region.size());
    }
    std::transform(process_region_ref.begin(), process_region_ref.end(), process_region_ref.begin(), ::toupper);
    auto ref = std::string_view{process_region_ref};

    std::size_t ref_index{};
    {
        std::ifstream ifs(bam_path);
        bam::Header bam_header(ifs);
        auto& ref_vec = bam_header.get_member<bam::HEADER_INDEX::REFERENCE>();
        for (std::size_t i = 0; i < ref_vec.size(); i++)
        {
            if (std::get<bam::REFERENCE_INDEX::REFERENCE_NAME>(ref_vec[i]) == chromosome)
            {
                ref_index = i;
                break;
            }
        }
    }

    auto windows_number = process_region.size() / assembly_region_size + (process_region.size() % assembly_region_size ? 1 : 0);

    auto origin_region = Interval{process_region.contig, process_region.begin, process_region.begin + assembly_region_size};
    auto padded_region = Interval{process_region.contig, process_padded_region.begin, origin_region.end + assembly_region_padding_size};

    std::ofstream ofs(vcf_path);
    for (std::size_t i = 0; i < windows_number; i++)
    {
        auto region_ref = ref.substr(padded_region.begin - process_padded_region.begin, padded_region.size());

        std::ifstream ifs(bam_path);
        bam::Header bam_header(ifs);
        bam::BAM bam(bam_header);
        bam::BAI bai(bai_path);
        bai.set_region({ref_index, padded_region.begin, ref_index, padded_region.end - 1});
        bam::BAM::get_obj(ifs, bam, bai);
        std::vector<SAMRecord> reads;
        while (bam)
        {
            auto bam_string = bam.to_string();
            auto clean_bam_string = bam_string.substr(0, static_cast<std::size_t>(
                std::distance(bam_string.begin(), boost::find_nth(bam_string, "\t", 10).end())) - 1);
            std::istringstream iss(clean_bam_string);
            SAMRecord record;
            iss >> record;
            reads.push_back(record);
            bam::BAM::get_obj(ifs, bam, bai);
        }
        if (reads.empty()) std::cout << "Ignore " << origin_region.to_string() << ":    (with overlap region = " << padded_region.to_string() << ")\n";
        else call_region(reads, region_ref, padded_region, origin_region, max_reads_on_assembly_region, ofs);

        origin_region.begin += assembly_region_size;
        origin_region.end   += assembly_region_size;
        padded_region.begin  = origin_region.begin - assembly_region_padding_size;
        padded_region.end    = origin_region.end   + assembly_region_padding_size;
    }
    std::cout << "HaplotypeCaller done." << '\n';
}

int main(int argc, char *argv[])
{
    do_work(argv[1]);
    //do_work("chrM:0-16571");
}

void filter_reads(std::vector<SAMRecord>& reads)
{
    reads.erase(std::remove_if(reads.begin(), reads.end(),
        MappingQualityReadFilter{}
    ), reads.end());
    reads.erase(std::remove_if(reads.begin(), reads.end(),
        DuplicateReadFilter{}
    ), reads.end());
    reads.erase(std::remove_if(reads.begin(), reads.end(),
        SecondaryAlignmentReadFilter{}
    ), reads.end());
    reads.erase(std::remove_if(reads.begin(), reads.end(),
        MateOnSameContigReadFilter{}
    ), reads.end());
}

void hard_clip_reads(std::vector<SAMRecord>& reads, const Interval& padded_region)
{
    std::for_each(reads.begin(), reads.end(), ReadClipper::hard_clip_soft_clipped_bases);
    std::for_each(reads.begin(), reads.end(), [&](auto& read){
        ReadClipper::hard_clip_to_interval(read, padded_region);
    });
    reads.erase(std::remove_if(reads.begin(), reads.end(),
        MinimumLengthReadFilter{}
    ), reads.end());
}

void sample_reads(std::vector<SAMRecord>& reads, std::size_t n)
{
    if (reads.size() > n)
    {
        std::vector<SAMRecord> temp;
        std::sample(reads.begin(), reads.end(), std::back_inserter(temp), n, std::mt19937{std::random_device{}()});
        reads.swap(temp);
    }
}

void call_region(std::vector<SAMRecord>& reads,
                 std::string_view ref,
                 const Interval& padded_region,
                 const Interval& origin_region,
                 std::size_t max_reads_on_assembly_region,
                 std::ostream& os)
{
    Assembler assembler;
    PairHMM pairhmm;
    Genetyper genetyper;

    filter_reads(reads);
    hard_clip_reads(reads, padded_region);
    sample_reads(reads, max_reads_on_assembly_region);
    if (reads.empty()) return;
    std::cout << "----------------------------------------------------------------------------------\n";
    std::cout << "Assembling " << origin_region.to_string() << " with " << reads.size() << " reads:    (with overlap region = " << padded_region.to_string() << ")\n";

    auto haplotypes = assembler.assemble(reads, ref);
    if (haplotypes.size() <= 1) return;

    auto likelihoods = pairhmm.compute_likelihoods(haplotypes, reads);
    auto variants = genetyper.assign_genotype_likelihoods(reads, haplotypes, likelihoods, ref, padded_region, origin_region);
    for (const auto& variant : variants)
        os << variant.to_string() << '\n';
}
