#pragma once

#include "interval.hpp"
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

    friend bool operator< (const Variant& lhs, const Variant& rhs)
    { return std::tie(lhs.location, lhs.REF, lhs.ALT) <  std::tie(rhs.location, rhs.REF, rhs.ALT); }

    friend bool operator==(const Variant& lhs, const Variant& rhs)
    { return std::tie(lhs.location, lhs.REF, lhs.ALT) == std::tie(rhs.location, rhs.REF, rhs.ALT); }

    std::string to_string() const
    {
        std::ostringstream oss;
        oss << location.contig << '\t'
           << location.begin + 1 << '\t'
           << "." << '\t'
           << (alleles.empty() ? "." : alleles[0]) << '\t';
        for (std::size_t i = 1; i < alleles.size(); i++)
            oss << alleles[i] << (i == alleles.size() - 1 ? '\t' : ',');
        oss << "." << '\t'
           << "." << '\t'
           << "." << '\t'
           << "GT:GQ" << '\t'
           << GT.first << '/' << GT.second << ':' << GQ;
        return oss.str();
    }

    auto size() const { return location.size(); }
};

} // hc