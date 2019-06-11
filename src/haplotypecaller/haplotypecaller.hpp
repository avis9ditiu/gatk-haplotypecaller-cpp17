#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <random>
#include "sam/sam.hpp"
#include "fasta/fasta.hpp"
#include "utils/interval.hpp"
#include "utils/read_filter.hpp"
#include "assembler/assembler.hpp"
#include "utils/read_clipper.hpp"
#include "pairhmm/intel_pairhmm.hpp"
#include "genotyper/genotyper.hpp"

namespace hc
{

class HaplotypeCaller
{
private:
    auto load_all_reads(std::size_t ref_size) 
    {
        std::ifstream ifs(in_path);
        assert(ifs);

        std::vector<std::vector<SAMRecord>> reads_map;
        reads_map.resize(ref_size);

        std::string line;
        // read header
        while (std::getline(ifs, line) && line[0] == '@');
        do {
            std::istringstream iss(line);
            SAMRecord record;
            iss >> record;
            reads_map[record.get_alignment_begin()].emplace_back(std::move(record));
        } while (std::getline(ifs, line));
        return reads_map;
    }

    auto select_one_read(const std::vector<SAMRecord>& reads)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, reads.size()-1);
        return reads[dis(gen)];
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
        std::for_each(reads.begin(), reads.end(),[](auto& read){
            // if (read.has_well_defined_fragment_size())
                ReadClipper::revert_soft_clipped_bases(read);
            // else ReadClipper::hard_clip_soft_clipped_bases(read);
        });
        std::for_each(reads.begin(), reads.end(), [&](auto& read){
            ReadClipper::hard_clip_to_interval(read, padded_region);
        });
        reads.erase(std::remove_if(reads.begin(), reads.end(),
            MinimumLengthReadFilter{}
        ), reads.end());
    }

    void call_region(std::vector<SAMRecord>& reads,
                    std::string_view ref,
                    const Interval& padded_region,
                    const Interval& origin_region,
                    std::ostream& os)
    {
        Assembler assembler;
        IntelPairHMM pairhmm;
        Genetyper genetyper;

        filter_reads(reads);
        hard_clip_reads(reads, padded_region);
        
        if (reads.empty()) return;
        std::cout << "----------------------------------------------------------------------------------\n";
        std::cout << "Assembling " << origin_region.to_string() << " with " << reads.size() << " reads:    (with overlap region = " << padded_region.to_string() << ")\n";

        auto haplotypes = assembler.assemble(reads, ref);
        if (haplotypes.size() <= 1) return;

        auto likelihoods = pairhmm.compute_likelihoods(haplotypes, reads);
        auto variants = genetyper.assign_genotype_likelihoods(reads, haplotypes, likelihoods, ref, padded_region, origin_region);
        for (const auto& variant : variants)
            variant.print(os);
    }
    
public:
    std::string in_path, out_path, ref_path;

    void do_work(std::size_t region_size = 245,
                 std::size_t padding_size = 85)
    {
        std::ifstream ifs(ref_path);
        assert(ifs);

        auto fasta = Fasta{};
        ifs >> fasta;
        ifs.close();

        std::transform(fasta.seq.begin(), fasta.seq.end(), fasta.seq.begin(), ::toupper);
        auto ref = std::string_view{fasta.seq};

        auto windows_number = (ref.size() + region_size - 1) / region_size;
        auto origin_region = Interval{fasta.name, 0, region_size};
        auto padded_region = origin_region;
        padded_region.end += padding_size;

        auto ofs = std::ofstream{out_path};
        assert(ofs);
        ofs << "##fileformat=VCFv4.2\n";
        ofs << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
        ofs << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        ofs << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878\n";

        auto reads_map = load_all_reads(ref.size());
        for (auto i = 0u; i < windows_number; i++)
        {
            std::vector<SAMRecord> reads;
            for (auto begin = padded_region.begin; begin != padded_region.end; begin++)
                if (!reads_map[begin].empty() && begin < reads_map.size())
                    reads.emplace_back(select_one_read(reads_map[begin]));
            
            if (reads.empty()) std::cout << "Ignore " << origin_region.to_string() << ":    (with overlap region = " << padded_region.to_string() << ")\n";
            else call_region(reads, ref.substr(padded_region.begin, padded_region.size()), padded_region, origin_region, ofs);

            origin_region.begin += region_size;
            origin_region.end   += region_size;
            padded_region.begin  = origin_region.begin - padding_size;
            padded_region.end    = origin_region.end   + padding_size;
        }
        std::cout << "HaplotypeCaller done." << '\n';
    }
};

}