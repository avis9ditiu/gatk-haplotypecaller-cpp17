#pragma once

#include <istream>
#include <ostream>
#include <string>
#include <string_view>
#include <stdexcept>

namespace hc
{

struct Fasta
{
    std::string name, comment, seq;

    void clear()
    {
        name.clear();
        comment.clear();
        seq.clear();
    }

    friend
    std::istream& operator>>(std::istream& is, Fasta& fa)
    {
        fa.clear();

        std::string line;
        getline(is, line);
        if (line.front() != '>')
            throw std::runtime_error("error: expected '>'");

        auto header = std::string_view(line).substr(1);
        auto pos = header.find_first_of(" \t\v\f\r");
        fa.name = header.substr(0, pos);
        if (pos != std::string::npos)
            fa.comment = header.substr(pos + 1);

        while (std::getline(is, line))
        {
            fa.seq.append(line);
            if (is.peek() == '>') break;
        }

        return is;
    }

    friend
    std::ostream& operator<<(std::ostream& os, Fasta& fa)
    {
        os << '>' << fa.name << (fa.comment.empty() ? "" : " ") << fa.comment << '\n';
        auto seq_view = std::string_view(fa.seq);
        for (std::size_t pos = 0; pos < seq_view.size(); pos += 50)
            os << seq_view.substr(pos, 50) << '\n';

        return os;
    }
};

}