#pragma once

#include <vector>
#include <string>
#include <ostream>

using Cigar = std::vector<std::pair<std::size_t, char>>;

std::string to_string(const Cigar& cigar)
{
    std::string cigar_string;
    for (auto [size, op] : cigar)
        cigar_string += std::to_string(size) += op;
    return cigar_string;
}

Cigar to_cigar(const std::string& cigar_string)
{
    Cigar cigar;
    for (std::size_t i = 0; i < cigar_string.size(); i++)
    {
        auto size = static_cast<size_t>(cigar_string[i] - '0');
        for (i++; std::isdigit(cigar_string[i]); i++)
            size = size * 10 + cigar_string[i] - '0';
        cigar.emplace_back(size, cigar_string[i]);
    }
    return cigar;
}

std::ostream& operator<<(std::ostream& os, const Cigar& cigar)
{
    for (auto [size, op] : cigar)
        os << size << op;
    return os;
}