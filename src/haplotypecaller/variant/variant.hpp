#pragma once

#include "../utils/interval.hpp"
#include <string>
#include <sstream>
#include <vector>

namespace hc
{

struct Variant
{
    Interval location;
    std::string REF;
    std::string ALT;
    std::vector<std::string> alleles;
    std::pair<std::size_t, std::size_t> GT;
    std::size_t GQ = 0;

    Variant() = default;
    Variant(Interval location, std::vector<std::string> alleles, std::pair<std::size_t, std::size_t> GT, std::size_t GQ)
        : location(std::move(location)), alleles(std::move(alleles)), GT(GT), GQ(GQ) {}


    friend bool operator< (const Variant& lhs, const Variant& rhs)
    { return std::tie(lhs.location, lhs.REF, lhs.ALT) <  std::tie(rhs.location, rhs.REF, rhs.ALT); }

    friend bool operator==(const Variant& lhs, const Variant& rhs)
    { return std::tie(lhs.location, lhs.REF, lhs.ALT) == std::tie(rhs.location, rhs.REF, rhs.ALT); }

    void print(std::ostream& os) const
    {
        os << location.contig << '\t'
           << location.begin + 1 << '\t'
           << "." << '\t'
           << alleles[0] << '\t';
        for (std::size_t i = 1; i < alleles.size(); i++)
            os << alleles[i] << (i == alleles.size() - 1 ? '\t' : ',');
        os << "." << '\t'
           << "." << '\t'
           << "." << '\t'
           << "GT:GQ" << '\t'
           << GT.first << '/' << GT.second << ':' << GQ << '\n';
    }

    auto size() const { return location.size(); }
    bool is_snp() const { return REF.size() == ALT.size(); }
    bool is_ins() const { return REF.size() <  ALT.size(); }
    bool is_del() const { return REF.size() >  ALT.size(); }
};

} // hc