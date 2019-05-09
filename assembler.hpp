#pragma once

#include "quality_utils.hpp"
#include "haplotype.hpp"
#include "graph_wrapper.hpp"
#include "sam.hpp"
#include "interval.hpp"
#include <iostream>

namespace hc
{

struct Assembler
{
    static constexpr std::size_t KMER_SIZE_ITERATION_INCREASE = 10;
    static constexpr std::size_t MAX_KMER_ITERATIONS_TO_ATTEMPT = 6;
    static constexpr std::size_t MAX_UNIQUE_KMERS_COUNT_TO_DISCARD = 1000;

    Assembler() = default;
    Assembler(std::size_t starting_kmer_size, char min_phred33_quality_to_use)
        : starting_kmer_size(starting_kmer_size),
          min_base_quality_to_use(min_phred33_quality_to_use) {}

private:
    std::size_t starting_kmer_size = 25;
    char min_base_quality_to_use   = 10 + QualityUtils::ASCII_OFFSET;

    std::pair<std::vector<Haplotype>, bool>
    assemble(const std::vector<SAMRecord>& reads,
                   std::string_view ref,
                   std::size_t kmer_size, 
                   bool allow_duplicate_kmers_in_ref)
    {
        if (ref.size() < kmer_size) return {};
        
        if (!allow_duplicate_kmers_in_ref && !GraphWrapper::get_dup_kmers(ref, kmer_size).empty())
        {
            std::cout << "Not using kmer size of " << kmer_size << " in read threading assembler because reference contains non-unique kmers\n";
            return {{}, false};
        };
            
        GraphWrapper graph{kmer_size, min_base_quality_to_use, allow_duplicate_kmers_in_ref};
        
        graph.set_ref(ref);
        for (const auto& read : reads)
            graph.set_read(read);

        graph.build();
        
        if (graph.unique_kmers_count() > MAX_UNIQUE_KMERS_COUNT_TO_DISCARD)
        {
            std::cout << "Not using kmer size of " << kmer_size << " in read threading assembler because it has too much unique kmers\n";
            return {{}, true};
        }
        
        if (graph.has_cycles())
        {
            std::cout << "Not using kmer size of " << kmer_size << " in read threading assembler because it contains a cycle\n";
            return {{}, false};
        }

        if (graph.is_low_complexity())
        {
            std::cout << "Not using kmer size of " << kmer_size << " in read threading assembler because it does not produce a graph with enough complexity\n";
            return {{}, false};
        }
        std::cout << "Using kmer size of " <<  kmer_size << " in assembler\n";

        return {graph.find_paths(), false};
    }

public:
    auto assemble(const std::vector<SAMRecord>& reads, std::string_view ref)
    {
        auto [haplotypes, too_much_unique_kmers] = assemble(reads, ref, starting_kmer_size, false);
        if (haplotypes.empty() && !too_much_unique_kmers)
        {
            auto kmer_size = starting_kmer_size;
            for (std::size_t i = 1; i <= MAX_KMER_ITERATIONS_TO_ATTEMPT; i++)
            {
                kmer_size += KMER_SIZE_ITERATION_INCREASE;
                bool last_attempt = (i == MAX_KMER_ITERATIONS_TO_ATTEMPT);
                std::tie(haplotypes, too_much_unique_kmers) = assemble(reads, ref, kmer_size, last_attempt);
                if (!haplotypes.empty() || too_much_unique_kmers) break;
            }
        }
        return haplotypes;
    }
};

} // hc
