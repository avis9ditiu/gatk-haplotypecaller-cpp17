#pragma once

#include "../sam/cigar.hpp"
#include "native/avx2-smithwaterman.h"

namespace hc
{

struct IntelSWAligner
{
public:
    struct SWParameters
    {
        int w_match;
        int w_mismatch;
        int w_open;
        int w_extend;
    };

    // match=1, mismatch = -1/3, gap=-(1+k/3)
    static constexpr SWParameters ORIGINAL_DEFAULT{3, -1, -4, -3};
    static constexpr SWParameters STANDARD_NGS{25, -50, -110, -6};
    static constexpr SWParameters NEW_SW_PARAMETERS{200, -150, -260, -11};
    static constexpr SWParameters ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS{10, -15, -30, -5};

    static constexpr std::size_t MINIMAL_MISMATCH_TO_TOLERANCE = 2;
public:
    // offset, Cigar
    std::pair<size_t, Cigar> align(std::string_view ref,
                                   std::string_view alt,
                                   const SWParameters& params = NEW_SW_PARAMETERS)
    {
        if (ref.empty() || alt.empty())
            throw std::invalid_argument("Non-null sequences are required for the SW aligner");

        if (is_all_match(ref, alt))
            return {0, {1, {ref.size(), CigarOperator::M}}};

        auto [match, mismatch, open, extend] = params;
        int count{};
        auto refLength = ref.length(), altLength = alt.length();
        auto cigarArray = std::make_unique<char[]>(2 * std::max(refLength, altLength));
        return {runSWOnePairBT_avx2(match, mismatch, open, extend, (uint8_t*)ref.data(), (uint8_t*)alt.data(), refLength, altLength, 9, cigarArray.get(), (int16_t*)&count), Cigar(cigarArray.get())};
    }

private:
    bool is_all_match(std::string_view ref, std::string_view alt) const
    {
        if (alt.size() == ref.size())
        {
            std::size_t mismatch = 0;
            for (std::size_t i = 0; mismatch <= MINIMAL_MISMATCH_TO_TOLERANCE && i < ref.size(); i++)
                if (alt[i] != ref[i])
                    mismatch++;
            return mismatch <= MINIMAL_MISMATCH_TO_TOLERANCE;
        }
        return false;
    }
};

}