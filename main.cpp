#include <gtest/gtest.h>
#include "smithwaterman.hpp"

SWAligner aligner;
using namespace std::string_literals;

TEST(SWAligner, degenerate_alignment_with_indels_at_both_ends)
{
    auto ref = "TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"s;
    auto alt =               "ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"s;
    auto [offset, cigar] = aligner.align(ref, alt, SWAligner::STANDARD_NGS);
    EXPECT_EQ(offset, 14);
    EXPECT_EQ(to_string(cigar), "31M20S");
}

TEST(SWAligner, sub_string_match)
{
    auto ref = "AAACCCCC"s;
    auto alt = "CCCCC"s;
    auto [offset, cigar] = aligner.align(ref, alt, SWAligner::ORIGINAL_DEFAULT);
    EXPECT_EQ(offset, 3);
    EXPECT_EQ(to_string(cigar), "5M");
}

TEST(SWAligner, sub_string_match_long)
{
    auto ref = "ATAGAAAATAGTTTTTGGAAATATGGGTGAAGAGACATCTCCTCTTATGGAAAAAGGGATTCTAGAATTTAACAATAAATATTCCCAACTTTCCCCAAGGCTTTAAAATCTACCTTGAAGGAGCAGCTGATGTATTTCTAGAACAGACTTAGGTGTCTTGGTGTGGCCTGTAAAGAGATACTGTCTTTCTCTTTTGAGTGTAAGAGAGAAAGGACAGTCTACTCAATAAAGAGTGCTGGGAAAACTGAATATCCACACACAGAATAATAAAACTAGATCCTATCTCTCACCATATACAAAGATCAACTCAAAACAAATTAAAGACCTAAATGTAAGACAAGAAATTATAAAACTACTAGAAAAAAACACAAGGGAAATGCTTCAGGACATTGGC"s;
    auto alt = "AAAAAAA"s;
    auto [offset, cigar] = aligner.align(ref, alt, SWAligner::ORIGINAL_DEFAULT);
    EXPECT_EQ(offset, 359);
    EXPECT_EQ(to_string(cigar), "7M");
}

TEST(SWAligner, complex_read_aligned_to_ref)
{
    auto ref = "AAAGGACTGACTG"s;
    auto alt =  "ACTGACTGACTG"s;
    auto [offset, cigar] = aligner.align(ref, alt, SWAligner::ORIGINAL_DEFAULT);
    EXPECT_EQ(offset, 1);
    EXPECT_EQ(to_string(cigar), "12M");
}

TEST(SWAligner, odd_no_alignment)
{
    auto ref = "AAAGACTACTG"s;
    auto alt = "AACGGACACTG"s;
    auto [offset1, cigar1] = aligner.align(ref, alt, {50, -100, -220, -12});
    EXPECT_EQ(offset1, 1);
    EXPECT_EQ(to_string(cigar1), "2M2I3M1D4M");

    auto [offset2, cigar2] = aligner.align(ref, alt, {200, -50, -300, -22});
    EXPECT_EQ(offset2, 0);
    EXPECT_EQ(to_string(cigar2), "11M");
}

TEST(SWAligner, indels_at_start_and_end)
{
    auto ref = "AAACCCCC"s;
    auto alt = "CCCCCGGG"s;
    auto [offset, cigar] = aligner.align(ref, alt, SWAligner::ORIGINAL_DEFAULT);
    EXPECT_EQ(offset, 3);
    EXPECT_EQ(to_string(cigar), "5M3S");
}

TEST(SWAligner, identical_alignments_with_differing_flank_lengths)
{
    auto padded_ref = "GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT"s;
    auto padded_alt = "GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA--GGGCC---------------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT"s;
    padded_alt.erase(std::remove(padded_alt.begin(), padded_alt.end(), '-'), padded_alt.end());

    auto not_padded_ref = "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA"s;
    auto not_padded_alt = "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA---------GGGCC--------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA"s;
    not_padded_alt.erase(std::remove(not_padded_alt.begin(), not_padded_alt.end(), '-'), not_padded_alt.end());

    auto SW_PAD = "NNNNNNNNNN"s;
    auto paddeds_ref = SW_PAD + padded_ref + SW_PAD;
    auto paddeds_alt = SW_PAD + padded_alt + SW_PAD;
    auto not_paddeds_ref = SW_PAD + not_padded_ref + SW_PAD;
    auto not_paddeds_alt = SW_PAD + not_padded_alt + SW_PAD;

    auto [padded_offset, padded_cigar] = aligner.align(paddeds_ref, paddeds_alt, SWAligner::NEW_SW_PARAMETERS);
    auto [not_padded_offset, not_padded_cigar] = aligner.align(not_paddeds_ref, not_paddeds_alt, SWAligner::NEW_SW_PARAMETERS);

    EXPECT_EQ(padded_cigar.size(), not_padded_cigar.size());
    for (std::size_t i = 0; i < padded_cigar.size(); i++)
    {
        auto [size1, op1] = padded_cigar[i];
        auto [size2, op2] = not_padded_cigar[i];
        if (op1 == 'M' && op2 == 'M') continue;
        EXPECT_EQ(size1, size2);
        EXPECT_EQ(op1, op2);
    }
}

int main()
{
    return RUN_ALL_TESTS();
}