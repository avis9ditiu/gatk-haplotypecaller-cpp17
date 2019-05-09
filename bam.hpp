/**
 *  @file sam_bam.hpp
 *  @brief The program parsing SAM and BAM format file
 *  @author JHHlab corp
 */
#pragma once

#include <map>
#include <tuple>
#include <vector>
#include <string>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <array>
#include <memory>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <type_traits>
#include <assert.h>
#include <zlib.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/qi_numeric.hpp>

namespace biovoltron::format {

/// Including BAM parser, BAI index, and related type declaration
    namespace bam
    {
        class Header;
        class BAM;
    };

/// Including SAM parser and related type declaration
    namespace sam
    {
        /// The abbreviation of a vector of string.
        using VecStr = std::vector<std::string>;

        /// The abbreviation of a pair of string.
        using PairStr = std::pair<std::string, std::string>;

        /// Enumerating different types of sorting.
        enum class SORT_TYPE
        {
            UNKNOWN
            , UNSORTED
            , QUERYNAME
            , COORDINATE
        };

        /// Mapping a string to a SORT_TYPE.
        std::map< std::string, SORT_TYPE > str_to_sort
            {
                  { "unknown",    SORT_TYPE::UNKNOWN      }
                , { "unsorted",   SORT_TYPE::UNSORTED     }
                , { "queryname",  SORT_TYPE::QUERYNAME    }
                , { "coordinate", SORT_TYPE::COORDINATE   }
            };

        /// Enumerating different types of grouping.
        enum class GROUP_TYPE
        {
            NONE
            , QUERY
            , REFERENCE
        };

        /// Mapping a string to a GROUP_TYPE.
        std::map<std::string, GROUP_TYPE> str_to_group
            {
                  { "none",     GROUP_TYPE::NONE      }
                , { "query",    GROUP_TYPE::QUERY     }
                , { "reference",GROUP_TYPE::REFERENCE }
            };

        /// Enumerating different types of platforms.
        enum class PLATFORM
        {
            CAPILLARY
            , LS454
            , ILLUMINA
            , SOLID
            , HELICOS
            , IONTORRENT
            , ONT
            , PABIO
        };

        /// Mapping a string to a PLATFORM type.
        std::map<std::string, PLATFORM> str_to_platform
            {
                  { "CAPILLARY",   PLATFORM::CAPILLARY    }
                , { "LS454",       PLATFORM::LS454        }
                , { "ILLUMINA",    PLATFORM::ILLUMINA     }
                , { "SOLID",       PLATFORM::SOLID        }
                , { "HELICOS",     PLATFORM::HELICOS      }
                , { "IONTORRENT",  PLATFORM::IONTORRENT   }
                , { "ONT",         PLATFORM::ONT          }
                , { "PABIO",       PLATFORM::PABIO        }
            };

        /// Indices of SQ field in header.
        enum REFERENCE_INDEX
        {
            REFERENCE_NAME
            , REFERENCE_LENGTH
            , ALTERNATE_LOCUS
            , ALTERNATE_REFERENCE_NAME
            , GENOME_ASSEM_ID
            , SPECIES
        };

        /// Indices of RG field in header.
        enum READ_GROUP_INDEX
        {
            READ_GROUP_ID
            , BARCODE
            , READ_GROUP_DESCRIPTION
            , FLOW_ORDER
            , KEY_SEQ
            , LIBRARY
            , READ_GROUP_PLATFORM
        };

        /// Indices of PG field in header.
        enum PROGRAM_INDEX
        {
            PROGRAM_ID
            , PROGRAM_NAME
            , COMMAND
            , PROGRAM_DESCRIPTION
            , PROGRAM_VERSION
        };

        /// Indices of the whole header.
        enum HEADER_INDEX
        {
            VERSION
            , ALIGNMENT_SORT_ORDER
            , ALIGNMENT_GROUPING
            , REFERENCE
            , READ_GROUP
            , PROGRAM
            , COMMENT
            , PLAIN_TEXT
        };

        /**
         *  @brief Indices of shift amount of FLAG field
         *
         *  This is a enumerator to record how many bits should
         *  FLAG field in SAM/BAM right shift to get specific FLAG in this
         *  alignment. After shifting, LSB is the FLAG you desired.<BR>
         *  ex. (sam.getMember<sam::SAM_MEMBER_INDEX::FLAG>() >>
         *  sam::FLAG::UNMAPPED) % 2 == 0
         */
        enum FLAG
        {
            MULTI_SEG
            , EXACT_MATCH
            , UNMAPPED
            , NEXT_UNMAPPED
            , REVERSED
            , NEXT_REVERSED
            , FIRST_SEG
            , LAST_SEG
            , SECONDARY
            , QT_FAILED
            , DUPLICATED
            , SUPPLEMENTARY
        };

        /// Indices of the optional field in SAM class.
        enum OPTIONAL_FIELD_INDEX
        {
            TAG
            , VALUE_TYPE
            , VALUE
        };

        /// Indices of the member in SAM class.
        enum SAM_MEMBER_INDEX
        {
            QNAME
            , FLAG
            , RNAME
            , POS
            , MAPQ
            , CIGAR
            , RNEXT
            , PNEXT
            , TLEN
            , SEQ
            , QUAL
            , OPTIONAL_FIELDS
        };

        /// Data structure of SQ field in header class.
        using ReferenceType = std::tuple<
            std::string
            , int32_t
            , std::string
            , std::string
            , std::string
            , std::string
        >;

        /// Data structure of RG field in header class.
        using ReadGroupType = std::tuple<
            std::string
            , std::string
            , std::string
            , std::string
            , std::string
            , std::string
            , PLATFORM
        >;

        /// Data structure of PG field in header class.
        using ProgramType = std::tuple<
            std::string
            , std::string
            , std::string
            , std::string
            , std::string
        >;

        /// Data structure of the whole header class.
        using HeaderType = std::tuple<
            std::string
            , SORT_TYPE
            , GROUP_TYPE
            , std::vector<ReferenceType>
            , std::vector<ReadGroupType>
            , std::vector<ProgramType>
            , std::vector<std::string>
            , std::string
        >;

        /// Data structure of the optional field in SAM class.
        using OptionalFieldType = std::tuple<
            std::array<char,2>
            , char
            , std::string
        >;

        /// Data structure of the member data in SAM class.
        using SAMMemberType = std::tuple<
            std::string                     ///< QRNAME
            , int32_t                         ///< FLAG
            , std::string                     ///< RNAME
            , int32_t                         ///< POS
            , int32_t                         ///< MAPQ
            , std::string                     ///< CIGAR
            , std::string                     ///< RNEXT
            , int32_t                         ///< PNEXT
            , int32_t                         ///< TLEN
            , std::string                     ///< SEQ
            , std::string                     ///< QUAL
            , std::vector<OptionalFieldType>  ///< OPTIONAL FIELD
        >;

        /**
         *  @brief The class of SAM header.
         *
         *  This class is essential to sam class.
         *  Which means that every sam class must have
         *  a header class.And each sam class will share
         *  the same header.
         */
        class Header
        {
            friend class SAM;
        public:
            /// Default constructor of header class.
            Header()
                : header_     ( HeaderType {
                std::string()
                , SORT_TYPE::UNKNOWN
                , GROUP_TYPE::NONE
                , std::vector<ReferenceType>()
                , std::vector<ReadGroupType>()
                , std::vector<ProgramType>()
                , std::vector<std::string>()
                , std::string()
            })
            {
                std::get<HEADER_INDEX::REFERENCE>(header_).reserve(3);
                std::get<HEADER_INDEX::COMMENT>(header_).reserve(1);
                std::get<HEADER_INDEX::PLAIN_TEXT>(header_).reserve(50);
            };

            /**
             * @brief Istream constructor of header class.
             * @param in Istream containing header information.
             *
             * Constructing Header class through istream.
             * This constuctor will parse the whole header
             * in a sam file.
             */
            Header(std::istream& in)
            {
                preparse(in);
            };

            /// Function template used to get a specific field in header.
            template <size_t n>
            const auto& get_member() const
            {
                return std::get<n>(header_);
            };

            /**
             * @brief Overloaded >> operator.
             * @param in Istream reference constaining header information.
             * @param rhs Target header instance.
             * @return Reference of istream.
             *
             * This overloaded operator is used to parse the
             * information in istream to a header.
             */
            friend std::istream& operator>>
                (std::istream& in, Header& rhs)
            {
                rhs.preparse(in);
                return in;
            }

            /**
             * @brief Overloaded << operator.
             * @param out Ostream reference outputing header content.
             * @param rhs Target header instance.
             * @return Reference of ostream.
             *
             * This overloaded operator is used to output the
             * information of a header.
             */
            friend std::ostream& operator<<
                (std::ostream& out, const Header& rhs)
            {
                out << rhs.to_string();
                return out;
            }

            /// Out put the whole header information.
            std::string to_string() const
            {
                return std::get<PLAIN_TEXT>(header_);
            }

            /**
             * @brief Public interface to pre-parse a header.
             * @param in Istream reference containing header information.
             *
             * This function is to recognize headers,classify them
             * and invoke corresponding handling function.
             */
            void preparse(std::istream& in)
            {
                VecStr lines;
                std::string line;
                int c = in.peek();
                while( c == '@')
                {
                    getline(in, line);
                    lines.emplace_back(std::move(line));
                    c = in.peek();
                }
                preparse_impl(lines);
            };

        protected:
            HeaderType header_; ///< This is a protected data member
            ///< because bam::Header has exactly
            ///< the same header contents

            /**
             * @brief This is a private function to parse header one by one.
             * @param lines Lines which splited from header
             *
             * This function is implementation details for parsing headers.
             * Set to protected member function is because
             * bam::Header also need it.
             */
            void preparse_impl(VecStr& lines)
            {
                VecStr vec;
                std::string& text = std::get<PLAIN_TEXT>(header_);
                for(size_t i = 0;i < lines.size();++i)
                {
                    boost::split(vec,lines[i], boost::is_any_of(" \t"),
                                 boost::token_compress_on);
                    if( vec[0] == "@HD" )
                        parse_HD(vec);
                    else if( vec[0] == "@SQ" )
                        parse_SQ(vec);
                    else if( vec[0] == "@RG" )
                        parse_RG(vec);
                    else if( vec[0] == "@PG" )
                        parse_PG(vec);
                    else if( vec[0] == "@CO" )
                        parse_CO(vec);
                    else
                        ;
                    text.append(lines[i]);
                    text.append("\n");
                }
            }

        private:
            /// Parsing HD field to header.
            void parse_HD(VecStr& vec)
            {
                for( size_t i(1) ; i<vec.size() ; ++i )
                {
                    PairStr tmp = split_colon(vec[i]);
                    if( tmp.first == "VN" )
                        std::get<VERSION>(header_) = tmp.second;

                    else if( tmp.first == "SO" )
                        std::get<ALIGNMENT_SORT_ORDER>(header_) =
                            str_to_sort[tmp.second];

                    else if( tmp.first == "GO" )
                        std::get<ALIGNMENT_GROUPING>(header_) =
                            str_to_group[tmp.second];

                    else
                        ;
                }
            }

            /// Parsing SQ field to header.
            void parse_SQ(VecStr& vec)
            {
                ReferenceType ret;
                for( size_t i(1) ; i<vec.size() ; ++i )
                {
                    PairStr tmp = split_colon(vec[i]);
                    if( tmp.first == "SN" )
                    {
                        std::get<REFERENCE_NAME>(ret) = tmp.second;
                    }
                    else if( tmp.first == "LN" )
                    {
                        std::get<REFERENCE_LENGTH>(ret) =
                            std::stoi(tmp.second);
                    }
                    else if( tmp.first == "AH" )
                    {
                        std::get<ALTERNATE_LOCUS>(ret) = tmp.second;
                    }
                    else if( tmp.first == "AN" )
                    {
                        std::get<ALTERNATE_REFERENCE_NAME>(ret) = tmp.second;
                    }
                    else if( tmp.first == "AS" )
                    {
                        std::get<GENOME_ASSEM_ID>(ret) = tmp.second;
                    }
                    else if( tmp.first == "SP" )
                    {
                        std::get<SPECIES>(ret) = tmp.second;
                    }
                    else
                        ;
                }
                std::get<REFERENCE>(header_).emplace_back(ret);
            }

            /// Parsing RG field to header.
            void parse_RG(VecStr& vec)
            {
                ReadGroupType ret;
                for( size_t i(1) ; i<vec.size() ; ++i )
                {
                    PairStr tmp = split_colon(vec[i]);
                    if( tmp.first == "ID" )
                    {
                        std::get<READ_GROUP_ID>(ret) = tmp.second;
                    }
                    else if( tmp.first == "BC" )
                    {
                        std::get<BARCODE>(ret) = tmp.second;
                    }
                    else if( tmp.first == "DS" )
                    {
                        std::get<READ_GROUP_DESCRIPTION>(ret) = tmp.second;
                    }
                    else if( tmp.first == "FO" )
                    {
                        std::get<FLOW_ORDER>(ret) = tmp.second;
                    }
                    else if( tmp.first == "KS" )
                    {
                        std::get<KEY_SEQ>(ret) = tmp.second;
                    }
                    else if( tmp.first == "LB" )
                    {
                        std::get<LIBRARY>(ret) = tmp.second;
                    }
                    else if( tmp.first == "PL" )
                    {
                        std::get<READ_GROUP_PLATFORM>(ret) =
                            str_to_platform[tmp.second];
                    }
                    else
                        ;
                }
                std::get<READ_GROUP>(header_).emplace_back(ret);
            }

            /// Parsing PG field to header.
            void parse_PG(VecStr& vec)
            {
                ProgramType ret;
                for( size_t i(1) ; i<vec.size() ; ++i )
                {
                    PairStr tmp = split_colon(vec[i]);
                    if( tmp.first == "ID" )
                    {
                        std::get<PROGRAM_ID>(ret) = tmp.second;
                    }
                    else if( tmp.first == "PN" )
                    {
                        std::get<PROGRAM_NAME>(ret) = tmp.second;
                    }
                    else if( tmp.first == "CL" )
                    {
                        std::get<COMMAND>(ret) = tmp.second;
                    }
                    else if( tmp.first == "DS" )
                    {
                        std::get<PROGRAM_DESCRIPTION>(ret) = tmp.second;
                    }
                    else if( tmp.first == "VN" )
                    {
                        std::get<PROGRAM_VERSION>(ret) = tmp.second;
                    }
                    else
                        ;
                }
                std::get<PROGRAM>(header_).emplace_back(ret);
            }

            /// Parsing CO field to header.
            void parse_CO(VecStr& vec)
            {
                std::string ret;
                for( size_t i(1) ; i<vec.size() ; ++i )
                {
                    ret.append(vec[i]);
                    if( i != vec.size()-1 )
                        ret.append("\t");
                }
                std::get<COMMENT>(header_).emplace_back(ret);
            }

            /**
             * @brief Function to split string to tokens by colon
             * @param A reference to string
             */
            PairStr split_colon(std::string& str)
            {
                PairStr ret;
                size_t pos = str.find(":");
                ret.first = str.substr(0, pos);
                ret.second = str.substr(pos + 1);
                return ret;
            }
        };

        /**
         * @brief The main structure of SAM parser.
         *
         * A SAM class contains a shared header and a
         * record of data.There are serveral function
         * provided for user to manipulate or to
         * check the data.
         */
        class SAM
        {
            friend class bam::BAM;
        public:
            /**
             * @brief The normal constructor of SAM.
             * @param header The reference of the header.
             *
             * Constructing a SAM class via existent header.
             * All the SAM class share the same header.
             */
            SAM(Header& header)
                : header_       (header)
                , data_members_ ( SAMMemberType {
                    std::string()
                    , 0
                    , std::string()
                    , 0
                    , 0
                    , std::string()
                    , std::string()
                    , 0
                    , 0
                    , std::string()
                    , std::string()
                    , std::vector<OptionalFieldType>()
                })
            {
                std::get<SAM_MEMBER_INDEX::RNAME>(data_members_).reserve(23);
                std::get<SAM_MEMBER_INDEX::CIGAR>(data_members_).reserve(6);
                std::get<SAM_MEMBER_INDEX::SEQ>(data_members_).reserve(35);
                std::get<SAM_MEMBER_INDEX::QUAL>(data_members_).reserve(35);
                std::get<SAM_MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_)
                    .reserve(6);
            };

            /**
             * @brief Convert constructor from BAM to SAM.
             * @param bam Source bam::BAM object to construct SAM.
             *
             * Constructing a SAM class from bam::BAM object.
             * Header in this object will reference to bam's header.
             */
            SAM(const bam::BAM& bam);

            /// Function template used to get the specific member of SAM.
            template <size_t N>
            const auto& get_member() const
            {
                return std::get<N>(data_members_);
            };

            /// Function template used to set the specific member of SAM.
            template <size_t N>
            void set_member
                (const std::tuple_element_t<N, SAMMemberType>& rhs)
            {
                std::get<N>(data_members_) = rhs;
            };

            /// Return the header of SAM class.
            const Header& get_header() const
            {
                return header_;
            };

            /// The public interface to user.
            std::string to_string() const
            {
                return to_string<0>();
            }

            /**
             * @brief Get one data from istream.
             * @param is Istream reference containing data.
             * @param obj Target SAM class to be parsed to.
             * @return Return istream reference for continuous usage.
             *
             * Static member function used to get one data from istream
             * and parse it to SAM class.
             */
            static std::istream& get_obj(std::istream& is, SAM& obj)
            {
                std::string str;
                getline(is, str);
                if (str.size() != 0)
                {
                    VecStr vec;
                    vec.reserve(26);
                    char* c_str = new char[str.size() + 1];
                    char* savestr = nullptr;
                    memcpy(c_str, str.c_str(), str.size() + 1);
                    char* split = strtok_r(c_str, "\t ", &savestr);
                    while (split != nullptr)
                    {
                        vec.emplace_back(split);
                        split = strtok_r(nullptr, "\t ", &savestr);
                    }
                    parse_one_line_data<0>(vec, obj);
                }
                return is;
            };

            /**
             * @brief Output the whole information of the SAM file.
             * @param os Ostream reference to output data.
             * @param obj Vector of SAM class.
             *
             * Dump the whole data,including one header and several data.
             */
            static void dump(std::ostream& os, std::vector<SAM>& obj)
            {
                os << obj[0].header_;
                for(size_t i(0) ; i<obj.size() ; ++i)
                {
                    os << obj[i];
                    if(i != obj.size()-1)
                        os << '\n';
                }
            };

            /**
             * @brief Overloaded >> operator.
             * @param in Istream reference containing data.
             * @param rhs Target SAM class to be parsed to.
             * @return Istream reference for continuous usage.
             *
             * This overloaded operator is used to get data from istream
             * to the sam class.
             */
            friend std::istream& operator>> (std::istream& in, SAM& rhs)
            {
                get_obj(in, rhs);
                return in;
            };

            /**
             * @brief Overloaded << operator.
             * @param out Ostream reference to output data.
             * @param rhs Target SAM class to be printed.
             * @return Ostream reference for continuous usage.
             *
             * This overloaded operator is used to output data from
             * the SAM class.
             */
            friend std::ostream& operator<<
                (std::ostream& out, const SAM& rhs)
            {
                out << rhs.to_string<0>();
                return out;
            };

        private:
            Header& header_;                ///< Reference to the header
            ///< which this alignment
            ///< belongs to.
            SAMMemberType data_members_;    ///< A Tuple which
            ///< contains alignment
            ///< information

            /// Function template used to recursively
            /// print the content of SAM.
            template<size_t N>
            void print_all()
            {
                if constexpr( N == OPTIONAL_FIELDS )
                {
                    std::cout << opt_to_string() << '\n';
                    return;
                }
                else
                {
                    std::cout << std::get<N>(data_members_) << " ";
                    print_all<N+1>();
                }
            }

            /// Function template used to convert optional field
            /// in data into string.
            std::string opt_to_string() const
            {
                std::string ret;
                for( auto& item : std::get<OPTIONAL_FIELDS>(data_members_) )
                {
                    std::array<char,2>tag = std::get<TAG>(item);
                    ret.append( 1, tag[0] );
                    ret.append( 1, tag[1] );
                    ret.append(":");
                    ret.append( 1, std::get<VALUE_TYPE>(item));
                    ret.append(":");
                    ret.append( std::get<VALUE>(item));
                    ret.append("\t");
                }
                if( ret.size() != 0 )
                    ret.erase(ret.size()-1);
                return ret;
            }

            /**
             * @brief Function template used to parse one line data
                      to target object.
             * @param vec Vecstor of string contsaining one line data.
             * @param obj Target SAM class to be parsed to.
             *
             * It's a private function to recursively parse data to a sam
             * class. It invokes another private function(parse_optional).
             * The reason to declare is as static is that it will be
             * called by get_obj(...) which is a static member function.
             */
            template<size_t N>
            static void parse_one_line_data(VecStr& vec, SAM& obj)
            {
                if constexpr( N == OPTIONAL_FIELDS )
                    parse_optional(vec, obj);
                else
                {
                    if constexpr( std::is_integral<
                        std::tuple_element_t<N,
                            decltype(obj.data_members_)>
                    >::value )
                    {
                        int32_t res = 0;
                        if(vec[N].size() != 0){
                            auto i = vec[N].begin();
                            boost::spirit::qi::parse(i, vec[N].end(),
                                                     boost::spirit::int_, res);
                        }
                        std::get<N>(obj.data_members_) = res;
                    }
                    else
                        std::get<N>(obj.data_members_) = vec[N];
                    parse_one_line_data<N+1>(vec, obj);
                }
                return;
            }

            /// Static function used to parse optional field
            /// from istream to target sam class.
            static void parse_optional(VecStr& vec, SAM& obj)
            {
                std::get<OPTIONAL_FIELDS>(obj.data_members_).clear();
                for(size_t i(11) ; i<vec.size() ; ++i)
                {
                    std::string tmp = vec[i];
                    std::array<char,2> tag;
                    tag[0] = tmp[0];
                    tag[1] = tmp[1];
                    char vt = tmp[3];
                    std::get<OPTIONAL_FIELDS>(obj.data_members_).
                        emplace_back(tag, vt, tmp.substr(5));
                }
            }

            /**
             * @brief Function template used to convert sam data into string.
             *
             * Recursively add data to a string a return it.
             * The returned string represents one data in the sam file.
             */
            template<size_t N>
            std::string to_string() const
            {
                std::string ret;
                if constexpr( N == OPTIONAL_FIELDS )
                {
                    ret.append(opt_to_string());
                }
                else{
                    if constexpr(
                        std::is_same<
                            std::tuple_element_t< N, SAMMemberType>
                            , int32_t
                        >::value )
                        ret.append(
                            std::to_string(std::get<N>(data_members_)));
                    else
                        ret.append(std::get<N>(data_members_));
                    ret.append("\t");
                    ret.append(to_string<N+1>());
                }
                if constexpr ( N == QUAL )
                {
                    if (ret.back() == '\t')
                        ret.pop_back();
                }
                return ret;
            };
        };
    };//namespace

    namespace bam
    {
        /// A shortcut to SORT_TYPE in namespace sam.
        using SORT_TYPE = sam::SORT_TYPE;
        /// A shortcut to GROUP_TYPE in namespace sam.
        using GROUP_TYPE = sam::GROUP_TYPE;
        /// A shortcut to PLATFORM in namespace sam.
        using PLATFORM = sam::PLATFORM;
        /// A shortcut to REFERENCE_INDEX in namespace sam.
        using REFERENCE_INDEX = sam::REFERENCE_INDEX;
        /// A shortcut to READ_GROUP_INDEX in namespace sam.
        using READ_GROUP_INDEX = sam::READ_GROUP_INDEX;
        /// A shortcut to PROGRAM_INDEX in namespace sam.
        using PROGRAM_INDEX = sam::PROGRAM_INDEX;
        /// A shortcut to HEADER_INDEX in namespace sam.
        using HEADER_INDEX = sam::HEADER_INDEX;
        /// Indices of RegionType which used in BAI object.
        enum REGION_INDEX
        {
            LEFT_REF_ID, LEFT_POS, RIGHT_REF_ID, RIGHT_POS
        };
        /// A shortcut to OPTIONAL_FIELD_INDEX in namespace sam.
        using OPTIONAL_FIELD_INDEX = sam::OPTIONAL_FIELD_INDEX;
        /// Indices of CigarType stored in BAM
        enum CIGAR_INDEX
        {
            OP, LENGTH
        };
        /// Indices of BAMMemberType to access information
        /// in a BAM alignment.
        enum BAM_MEMBER_INDEX
        {
            QNAME
            , FLAG
            , REF_ID
            , POS
            , MAPQ
            , CIGAR
            , RNEXT_ID
            , PNEXT
            , TLEN
            , SEQ
            , QUAL
            , OPTIONAL_FIELDS
            , BIN
        };

        /// A shortcut to ReferenceType in namespace sam.
        using ReferenceType = sam::ReferenceType;
        /// A shortcut to ReadGroupType in namespace sam.
        using ReadGroupType = sam::ReadGroupType;
        /// A shortcut to ProgramType in namespace sam.
        using ProgramType = sam::ProgramType;
        /// A shortcut to HeaderType in namespace sam.
        using HeaderType = sam::HeaderType;
        /// Data structure to store desired reference region in BAI object.
        using RegionType = std::tuple<
            int64_t
            , int64_t
            , int64_t
            , int64_t
        >;
        /// A shortcut to OptionalFieldType in namespace sam.
        using OptionalFieldType = sam::OptionalFieldType;
        /// Data structure to store Cigar operations in BAM alignment.
        using CigarType = std::tuple<char, uint32_t>;
        /// Data structure to store information in BAM alignment.
        using BAMMemberType = std::tuple<
            std::string                     ///< QNAME
            , uint16_t                        ///< FLAG
            , int32_t                         ///< RNAME replaced with ID
            , int32_t                         ///< POS
            , uint8_t                         ///< MAPQ
            , std::vector<CigarType>          ///< CIGAR
            , int32_t                         ///< RNEXT replaced with ID
            , int32_t                         ///< PNEXT
            , int32_t                         ///< TLEN
            , std::string                     ///< SEQ
            , std::string                     ///< QUAL
            , std::vector<OptionalFieldType>  ///< Optional fields
            , uint16_t                        ///< BIN number in
            ///< BAM binning concepts
        >;

        /**
         * @brief Class to represents BAM file format's header.
         * @warning When there has any BAM object constructed from
         *          the header object, or there has any SAM object
         *          contruct from the BAM object which
         *          reference to the header object,
         *          the header MUST NOT destruct before they destruct.
         *
         * This class stores header in BAM file format.<BR>
         * Every BAM alignments from the same BAM file MUST reference
         * to the same Header object, or might be undefined behavior.
         */
        class Header : public sam::Header
        {
            friend class BAM;
        public:
            /// Default constructor to initialize every data members.
            Header()
                : block_offset_     (   0   )
                , block_length_     (   0   )
                , block_address_    (   0   )
            {
                init_zargs();
            };

            /**
             * @brief istream preparse constructor to
             *        parse header text in istream.
             * @param in istream object which
             *        contains header field in BAM file
             * @see Header(), preparse()
             *
             * This is a constructor can help you to parse
             * header information in BAM file simultaneously
             * when you construct Header object.<BR>
             * Same as call default constructor then call preparse function.
             */
            Header(std::istream& in)
                : block_offset_     (   0   )
                , block_length_     (   0   )
                , block_address_    (   0   )
            {
                init_zargs();
                preparse(in);
            };

            /**
             * @brief Convert constructor from sam::Header to bam::Header.
             * @param rhs Source sam::Header object.
             * @see sam::Header, bam::BAM(Header&, const sam::SAM&)
             *
             * A constructor convert sam::Header to bam::Header.<BR>
             * This constructor might be called when user want to
             * convert SAM object to BAM object because of the
             * header in BAM object can not directly<BR>
             * reference to the header in SAM object. So needs to create
             * a bam::Header to help new BAM object has
             * bam::Header to reference to.
             */
            Header(const sam::Header& rhs)
                : sam::Header       (  rhs  )
                , block_offset_     (   0   )
                , block_length_     (   0   )
                , block_address_    (   0   )
            {
                init_zargs();
            }

            Header(const Header& rhs) = default;
            Header(Header&& rhs) = default;
            Header& operator=(const Header& rhs) = default;
            Header& operator=(Header&& rhs) = default;
            ~Header() = default;

            /**
             *  @brief Clear all data in the object.
             *  @warning MUST NOT call this function when any
             *           SAM/BAM object is referencing to this header.
             *
             *  Clear all data in this header object.<BR>
             *  This function might be called when user want to open a
             *  new BAM file and want to reuse this object rather than
             *  destruct it.
             */
            void reset()
            {
                block_offset_ = 0;
                block_length_ = 0;
                block_address_ = 0;
                std::get<HEADER_INDEX::VERSION>(header_).clear();
                std::get<HEADER_INDEX::ALIGNMENT_SORT_ORDER>(header_) =
                    SORT_TYPE::UNKNOWN;
                std::get<HEADER_INDEX::ALIGNMENT_GROUPING>(header_) =
                    GROUP_TYPE::NONE;
                std::get<HEADER_INDEX::REFERENCE>(header_).clear();
                std::get<HEADER_INDEX::READ_GROUP>(header_).clear();
                std::get<HEADER_INDEX::PROGRAM>(header_).clear();
                std::get<HEADER_INDEX::COMMENT>(header_).clear();
                std::get<HEADER_INDEX::PLAIN_TEXT>(header_).clear();
            }

            /**
             *  @brief Parse header field in BAM file.
             *  @param in istream which contains header in BAM file.
             *
             *  Parse header information in istream. After parsing,
             *  header become usable.
             */
            void preparse(std::istream& in)
            {
                reset();
                // magic & l_text
                size_t data_capacity = 8;
                char* data = new char[8];
                read_byte_data(in, data, 8);
                if (memcmp(data, BAM_MAGIC, 4) != 0)
                {
                    std::cerr << "ERROR: BAM Magic string not match";
                    std::cerr << std::endl;
                    return;
                }
                // text
                int32_t l_text = convert_char<int32_t>(&data[4]);
                if (l_text + 1 > data_capacity)
                {
                    delete []data;
                    data = new char[l_text + 1];
                    data_capacity = l_text + 1;
                }
                read_byte_data(in, data, l_text);
                data[l_text] = '\0';
                std::string text(data);
                std::vector<std::string> lines;
                boost::split(lines, text, boost::is_any_of("\n"),
                             boost::token_compress_on);
                lines.pop_back();
                preparse_impl(lines);
                // n_ref
                read_byte_data(in, data, sizeof(int32_t));
                int32_t n_ref = convert_char<int32_t>(data);
                // references
                int32_t l_name, l_ref;
                std::string name;
                std::vector<ReferenceType>& refs =
                    std::get<HEADER_INDEX::REFERENCE>(header_);
                std::vector<ReferenceType>::iterator it;
                for (size_t i = 0;i < n_ref;++i)
                {
                    // name
                    read_byte_data(in, data, sizeof(int32_t));
                    l_name = convert_char<int32_t>(data);
                    if (l_name + 1 > data_capacity)
                    {
                        delete[] data;
                        data = new char[l_name + 1];
                        data_capacity = l_name + 1;
                    }
                    read_byte_data(in, data, l_name);
                    data[l_name] = '\0';
                    name = data;
                    // l_ref
                    read_byte_data(in, data, sizeof(int32_t));
                    l_ref = convert_char<int32_t>(data);
                    it = refs.begin();
                    while (it != refs.end())
                    {
                        if (std::get<REFERENCE_INDEX::REFERENCE_NAME>(*it)
                            == name)
                            break;
                        std::advance(it, 1);
                    }
                    if (it == refs.end())
                        refs.emplace_back(name, l_ref, "", "", "", "");
                }
                delete[] data;
            }

            /**
             *  @brief Overload >> operator.
             *  @param in istream which data comes from
             *  @param rhs Header object which will contains data
             *  @return Parameter in for continuous data flow
             *  @see preparse()
             *
             *  Overload >> operator for user to use more intuitively.
             *  Same as call rhs.preparse().
             */
            friend std::istream& operator>>
                (std::istream& in, Header& rhs)
            {
                rhs.preparse(in);
                return in;
            };

            /**
             *  @brief Overload << operator.
             *  @param out ostream which data will direct to
             *  @param rhs Header object which want to output
             *  @return Parameter out for continuous data flow
             *  @see getMember()
             *
             *  Output header with BAM file format (binary).<BR>
             *  If want to read header in plaintext, use getMember() to
             *  get header plaintext in header data structure.
             */
            friend std::ostream& operator<<
                (std::ostream& out, const Header& rhs)
            {
                char data[CHUNK_SIZE];
                uint32_t offset = 0;
                // magic string
                memcpy(data, BAM_MAGIC, sizeof(BAM_MAGIC));
                offset += sizeof(BAM_MAGIC);
                // l_text & text
                const std::string& text =
                    std::get<HEADER_INDEX::PLAIN_TEXT>(rhs.header_);
                int32_t buf = text.size();
                if (buf > INT_MAX)
                {
                    std::cerr << "ERROR: Header text length is too long.";
                    std::cerr << std::endl;
                    return out;
                }
                parse_char<int32_t>(&data[offset], buf, &offset);
                if (offset + buf > CHUNK_SIZE)
                {
                    const char* c_str = text.c_str();
                    memcpy(&data[offset], c_str, CHUNK_SIZE - offset);
                    rhs.deflate_block(out, data, CHUNK_SIZE);
                    memcpy(data, &c_str[CHUNK_SIZE - offset],
                           buf - CHUNK_SIZE + offset);
                    offset = buf - CHUNK_SIZE + offset;
                }
                else
                {
                    memcpy(&data[offset], text.c_str(), buf);
                    offset += buf;
                }
                // n_ref & refs
                const std::vector<ReferenceType>& refs =
                    std::get<HEADER_INDEX::REFERENCE>(rhs.header_);
                buf = refs.size();
                parse_char<int32_t>(&data[offset], buf, &offset);
                for (size_t i = 0;i < refs.size();++i)
                {
                    const std::string& ref_name =
                        std::get<REFERENCE_INDEX::REFERENCE_NAME>(refs[i]);
                    buf = ref_name.size();
                    if (offset + 2 * sizeof(int32_t) + buf > CHUNK_SIZE)
                    {
                        rhs.deflate_block(out, data, offset);
                        offset = 0;
                    }
                    parse_char<int32_t>(&data[offset], buf, &offset);
                    memcpy(&data[offset], ref_name.c_str(), buf);
                    offset += buf;
                    parse_char<int32_t>
                        (&data[offset],
                         std::get<REFERENCE_INDEX::REFERENCE_LENGTH>(refs[i]),
                         &offset);
                }
                rhs.deflate_block(out, data, offset);
                return out;
            };

        private:
            // Need deflate_block()
            friend std::ostream& operator<<
                (std::ostream& out, const BAM& rhs);

            /// Initialize z_stream for future inflate_block() use.
            void init_zargs()
            {
                z_args.zalloc = Z_NULL;
                z_args.zfree = Z_NULL;
                z_args.opaque = Z_NULL;
                z_args.avail_in = 0;
                z_args.next_in = Z_NULL;
                inflateInit2(&z_args, GZIP_WINDOW_BITS);
            }

            /**
             *  @brief Convert data stream to desired data type.
             *  @tparam T The type which you want data to becomes
             *  @param data Data pointer points to the starting address
             *         of the data which you want to convert
             *  @return Converted data in desired data type
             *  @see parse_char(), determine_big_endianess()
             *
             *  A private function for internal use when read data stream
             *  from BAM file, then convert byte into variable in
             *  specific type.
             */
            template <typename T>
            static T convert_char(const char* data)
            {
                if (is_big_endian_)
                {
                    size_t size = sizeof(T);
                    std::unique_ptr<char[]> temp(new char[size]);
                    for (size_t i = 0;i < size;++i)
                        temp[i] = data[size - 1 - i];
                    return (*reinterpret_cast<T*>(temp.get()));
                }
                else
                    return (*reinterpret_cast<const T*>(data));
            }

            /**
             *  @brief Convert variable into byte (little endian).
             *  @tparam T The type which you want data output to data stream
             *  @tparam U Original type which data is (Might different from T)
             *  @param out Data stream which you want to output
             *  @param data Target variable which you want to output
             *  @param offset (Optional) Data stream index remainder,
             *                will increase by sizeof(T). Default is nullptr
             *  @see convert_char(), determine_big_endianess()
             *
             *  A private function for internal use when write variable
             *  to data stream in little endian form.<BR>
             *  Template parameters T and U might be different, but typically
             *  is the same. T is the type BAM file format provides,
             *  so U must converted to T to match the documentation.<BR>
             *  Offset is optional parameter that could provide to increase
             *  file pointer (out) outside the function for
             *  convenience reasons.
             */
            template <typename T, typename U>
            static void parse_char
                (char* out, U data, uint32_t* offset = nullptr)
            {
                size_t size = sizeof(T);
                if (is_big_endian_)
                {
                    for (size_t i = 0;i < size;++i)
                        out[i] = (char)(data >> 8 * i);
                }
                else
                    memcpy(out, &data, size);

                if (offset)
                    *offset += size;
            }

            /**
             *  @brief Check whether system is big endian
             *  @return Whether system is big endian
             *
             *  Because BAM file format is liitle endian,
             *  if system is big endian then program might work improperly.
             *  So needs this function to checking endianess
             *  in working system.
             */
            static bool determine_big_endianess()
            {
                const uint16_t x = 0x0001;
                return (*(reinterpret_cast<const char*>(&x))) == 0;
            }

            /**
             *  @brief Check whether BGZF header is valid.
             *  @param header Byte contains BGZF header information
             *  @return Whether BGZF header is valid
             *
             *  Because in BGZF file format, there has a specific
             *  gzip header format. So needs this function to determine
             *  the file is BGZF file or not.
             */
            bool check_header_is_valid(char header[18])
            {
                if (header[0]  != GZIP_ID1
                    || header[1]  != GZIP_ID2
                    || header[2]  != GZIP_CM
                    || header[3]  != GZIP_FLAG
                    || header[12] != GZIP_SI1
                    || header[13] != GZIP_SI2
                    || header[14] != GZIP_SLEN) return false;

                return true;
            }

            /**
             *  @brief Decompress BGZF block to buffer.
             *  @param in istream contains data in compressed format
             *  @return Byte count that decompresssed
             *  @see check_header_is_valid(), deflate_block()
             *
             *  Use zlib inflate() to decompress a block
             *  from BGZF file format.<BR>
             *  Decomressed data will put into block_buffer_ and update
             *  block_offset_, block_length, and block_address_.
             */
            uint32_t inflate_block(std::istream& in)
            {
                char header[18];
                block_address_ = in.tellg();
                in.read(header, sizeof(header));
                if (in.eof())
                    return -1;
                else if (!check_header_is_valid(header))
                {
                    std::cerr << "ERROR: Header format not match";
                    std::cerr << std::endl;
                    return -1;
                }

                int8_t ret;
                char in_char[CHUNK_SIZE];
                uint16_t block_size = convert_char<uint16_t>(&header[16]);
                z_args.avail_in = block_size - GZIP_XLEN - 19;
                in.read(in_char, z_args.avail_in);
                z_args.next_in = reinterpret_cast<Bytef*>(in_char);
                z_args.avail_out = CHUNK_SIZE;
                z_args.next_out = reinterpret_cast<Bytef*>(block_buffer_);
                ret = inflate(&z_args, Z_FINISH);
                if (ret != Z_STREAM_END)
                {
                    inflateEnd(&z_args);
                    std::cerr << "ERROR: inflate() failed" << std::endl;
                    return -1;
                }
                in.seekg(sizeof(uint32_t), std::ios::cur);
                uint32_t i_size;
                in.read(reinterpret_cast<char*>(&i_size), sizeof(i_size));
                if (i_size != z_args.total_out)
                {
                    std::cerr << "ERROR: ";
                    std::cerr << "ISIZE doesn't match total_out in zlib";
                    std::cerr << std::endl;
                    return -1;
                }
                ret = inflateReset(&z_args);
                if (ret != Z_OK)
                {
                    std::cerr << "ERROR: inflateReset() failed";
                    std::cerr << std::endl;
                    return -1;
                }

                block_offset_ = 0;
                block_length_ = i_size;
                return z_args.total_out;
            }

            /**
             *  @brief Compress data into BGZF block.
             *  @param outfile ostream which compressed data will forward to
             *  @param data Data stream which will be compressed
             *  @param len Data length in data stream
             *  @return Byte count after compresssion
             *  @see check_header_is_valid(), inflate_block()
             *
             *  Use zlib deflate() to compress data into
             *  a block in BGZF file format.<BR>
             *  This is a static function is because of there has no
             *  "buffer" in data members, just pass data stream which
             *  want to compress then every one can use this function.
             */
            static uint32_t deflate_block
                (std::ostream& outfile, char* data, uint32_t len)
            {
                int8_t ret, flush;
                uint32_t out_length;
                z_stream dz_args;
                char out[CHUNK_SIZE + GZIP_WRAP_DATA_SIZE];
                out[0] = GZIP_ID1;
                out[1] = GZIP_ID2;
                out[2] = GZIP_CM;
                out[3] = GZIP_FLAG;
                out[10] = GZIP_XLEN;
                out[11] = '\0';
                out[12] = GZIP_SI1;
                out[13] = GZIP_SI2;
                out[14] = GZIP_SLEN;
                out[15] = '\0';

                dz_args.zalloc = Z_NULL;
                dz_args.zfree = Z_NULL;
                dz_args.opaque = Z_NULL;
                dz_args.next_in = reinterpret_cast<Bytef*>(data);
                dz_args.avail_in = len;
                dz_args.next_out = reinterpret_cast<Bytef*>(&out[18]);
                dz_args.avail_out = CHUNK_SIZE;
                ret = deflateInit2(&dz_args
                , Z_DEFAULT_COMPRESSION
                , Z_DEFLATED
                , GZIP_WINDOW_BITS
                , GZIP_MEM_LEVEL
                , Z_DEFAULT_STRATEGY);
                if (ret != Z_OK)
                {
                    std::cerr << "ERROR: deflateinit2() failed";
                    std::cerr << std::endl;
                    return -1;
                }
                ret = deflate(&dz_args, Z_FINISH);
                if (ret != Z_STREAM_END)
                {
                    deflateEnd(&dz_args);
                    std::cerr << "ERROR: deflate() failed";
                    std::cerr << std::endl;
                    return -1;
                }
                ret = deflateEnd(&dz_args);
                if (ret != Z_OK)
                {
                    std::cerr << "ERROR: deflateEnd() failed";
                    std::cerr << std::endl;
                    return -1;
                }

                out_length = dz_args.total_out + GZIP_WRAP_DATA_SIZE;
                if (out_length > CHUNK_SIZE)
                {
                    std::cerr << "ERROR: Compressed block isze > ";
                    std::cerr << "BGZF maximum block size" << std::endl;
                    return -1;
                }
                parse_char<uint16_t>(&out[16], (uint16_t)out_length - 1);
                uint32_t gzip_crc = crc32(0, Z_NULL, 0);
                gzip_crc = crc32(gzip_crc
                    , reinterpret_cast<Bytef*>(data), len);
                parse_char<uint32_t>(&out[out_length] - 8, gzip_crc);
                parse_char<uint32_t>(&out[out_length] - 4, len);
                outfile.write(out, out_length);
                return out_length;
            }

            /**
             *  @brief Read specific number of byte from block_buffer_.
             *  @param in istream which contains file pointer to BAM file
             *  @param data A buffer to store result and return to caller
             *  @param size Byte count that caller want to read
             *  @return Whether read success
             *  @see inflate_block()
             *
             *  Read binary data from decompressed data in block_buffer_.<BR>
             *  If buffer has not enough data, function will call
             *  inflate_block() to get more data from BGZF file.
             */
            bool read_byte_data
                (std::istream& in, char* data, uint64_t size)
            {
                uint32_t min, available = block_length_ - block_offset_;
                do
                {
                    if (available == 0)
                    {
                        if (inflate_block(in) == -1)
                            return false;
                        available = block_length_;
                    }
                    min = (available > size) ? size : available;
                    memcpy(data, block_buffer_ + block_offset_, min);
                    block_offset_ += min;
                    available -= min;
                    size -= min;
                    data += min;
                } while (size != 0);
                return true;
            }

            /**
             *  @brief Seek file pointer to given virtual file offset.
             *  @param in istream which contains file pointer to BGZF file
             *  @param offset virtual file offset in BGZF documentation
             *  @see tell()
             *
             *  Seek istream's file pointer to virtual file offset.<BR>
             *  Might be called when random access is activated,
             *  then istream's file pointer will move to starting address
             *  of a block and inflate_block() will be invoked to
             *  decompress block. Finallly block_offset_ will seek to
             *  virtual file offset.
             */
            void seek(std::istream& in, uint64_t offset)
            {
                in.seekg(offset >> BAM_OFFSET_SHIFT, std::ios::beg);
                if (!in)
                {
                    std::cerr << "ERROR: ";
                    std::cerr << "coffset seek to wrong file position";
                    std::cerr << std::endl;
                    return;
                }
                else if (inflate_block(in) == -1)
                    return;
                else if ((offset & 0xffff) > block_length_)
                {
                    std::cerr << "ERROR: ";
                    std::cerr << "uoffset seek to wrong file position";
                    std::cerr << std::endl;
                    return;
                }
                block_offset_ = offset & 0xffff;
            }

            /**
             *  @brief Tell caller the current virtual file offset.
             *  @return Current virtual file offset
             *  @see seek()
             */
            uint64_t tell()
            {
                return (block_address_ << BAM_OFFSET_SHIFT | block_offset_);
            }

            /// BGZF file format's ID1 field
            static const char GZIP_ID1 = 31;
            /// BGZF file format's ID2 field
            static const char GZIP_ID2 = static_cast<char>(139);
            /// BGZF file format's CM field
            static const char GZIP_CM = 8;
            /// BGZF file format's FLAG field
            static const char GZIP_FLAG = 4;
            /// BGZF file format's XLEN field
            static const char GZIP_XLEN = 6;
            /// BGZF file format's SI1 field
            static const char GZIP_SI1 = 66;
            /// BGZF file format's SI2 field
            static const char GZIP_SI2 = 67;
            /// BGZF file format's SLEN field
            static const char GZIP_SLEN = 2;
            /// Represents raw inflate/deflate for zlib
            static const int8_t GZIP_WINDOW_BITS = -15;
            /// Default memory level in zlib
            static const int8_t GZIP_MEM_LEVEL = 8;
            /// BGZF file format's block size without compresssed data size
            static const uint8_t GZIP_WRAP_DATA_SIZE = 26;
            /// BGZF file format's maximal block size
            static const uint32_t CHUNK_SIZE = 65536;
            /// The bits need to shift for block starting address in
            /// BGZF's virtual file offset
            static const uint8_t BAM_OFFSET_SHIFT = 16;
            /// BAM file format's magic string
            static constexpr char BAM_MAGIC[4] = {'B', 'A', 'M', '\1'};

            /// Records offset in this block
            uint32_t block_offset_;
            /// Records this block's length
            uint32_t block_length_;
            /// Records this block's starting address
            uint64_t block_address_;
            /// Buffers decompressed data in this block
            char block_buffer_[CHUNK_SIZE];
            /// zlib arguments
            z_stream z_args;
            /// Records system's endianess
            static bool is_big_endian_;
        };
        bool Header::is_big_endian_ = Header::determine_big_endianess();

        // Template specification for float type
        // because of float endianess convertion is different from other type.
        template <>
        void Header::parse_char<float, float>
            (char* out, float data, uint32_t* offset)
        {
            if (is_big_endian_)
            {
                char* data_ptr = reinterpret_cast<char*>(&data);
                out[0] = data_ptr[3];
                out[1] = data_ptr[2];
                out[2] = data_ptr[1];
                out[3] = data_ptr[0];
            }
            else
                memcpy(out, &data, sizeof(float));

            if (offset)
                *offset += sizeof(float);
        }

        /**
         * @brief Class to store BAI file index information
         * @warning Only BAI file read successfully then could active
         * random access or an error might occur.
         *
         * A class to store indexing information in BAI file to
         * accomplish fast random access in BAM file.<BR>
         * Class also provides functions to set user's interested
         * reference region to help locate right file offset.
         */
        class BAI
        {
            friend class BAM;
        public:
            /**
             *  @brief Filename constructor to read BAI file.
             *  @param filename BAI file's filename
             *  @see load()
             *
             *  Construct object and load index information simultaneously.
             *  Equlas to call load() to try to load BAI file into object.<BR>
             *  This is the only way to construct BAI file. An empty
             *  BAI object is meanless.
             */
            BAI(std::string filename)
                : region_           ( RegionType{-1, -1, -1, -1} )
                , current_index_    (             -1             )
                , n_no_coor_        (              0             )
            {
                load(filename);
            };

            BAI(const BAI& rhs) = default;
            BAI(BAI&& rhs) = default;
            BAI& operator=(const BAI& rhs) = default;
            BAI& operator=(BAI&& rhs) = default;
            ~BAI() = default;

            /// Clear all data in this object
            void reset()
            {
                index_info_.clear();
                overlap_chunks_.clear();
                region_ = RegionType{-1, -1, -1, -1};
                current_index_ = -1;
                n_no_coor_ = 0;
            }

            /**
             *  @brief Load index information from BAI file.
             *  @param filename BAI file's filename
             *  @return Whether information load sucessfully
             */
            bool load(std::string filename)
            {
                reset();
                std::ifstream infile(filename, std::ios::binary);
                if (!infile)
                {
                    std::cerr << "ERROR: File not open" << std::endl;
                    return false;
                }
                // magic
                char magic[4];
                infile.read(magic, 4);
                if (memcmp(magic, BAI_MAGIC, 4) != 0)
                {
                    std::cerr << "ERROR: BAI magic string not match";
                    std::cerr << std::endl;
                    return false;
                }
                // n_ref
                int32_t n_ref, n_bin, n_chunk, n_intv;
                uint32_t bin_num;
                uint64_t chunk_beg, chunk_end, ioffset;
                BinType bin_buf;
                LinearIndexType index_buf;
                std::vector<ChunkType> chunks_buf;
                infile.read(reinterpret_cast<char*>(&n_ref)
                    , sizeof(int32_t));
                for (size_t i = 0;i < n_ref;++i)
                {
                    // n_bin
                    infile.read(reinterpret_cast<char*>(&n_bin)
                        , sizeof(int32_t));
                    for (size_t j = 0;j < n_bin;++j)
                    {
                        // bin
                        infile.read(reinterpret_cast<char*>(&bin_num)
                            , sizeof(uint32_t));
                        // n_chunk
                        infile.read(reinterpret_cast<char*>(&n_chunk)
                            , sizeof(int32_t));
                        for (size_t k = 0;k < n_chunk;++k)
                        {
                            infile.read(reinterpret_cast<char*>(&chunk_beg)
                                , sizeof(uint64_t));
                            infile.read(reinterpret_cast<char*>(&chunk_end)
                                , sizeof(uint64_t));
                            chunks_buf.emplace_back(chunk_beg, chunk_end);
                        }
                        bin_buf.emplace(bin_num, chunks_buf);
                        chunks_buf.clear();
                    }
                    // n_intv
                    infile.read(reinterpret_cast<char*>(&n_intv)
                        , sizeof(int32_t));
                    for (size_t j = 0;j < n_intv;++j)
                    {
                        // ioffset
                        infile.read(reinterpret_cast<char*>(&ioffset)
                            , sizeof(uint64_t));
                        index_buf.emplace_back(ioffset);
                    }
                    index_info_.emplace_back(bin_buf, index_buf);
                    index_buf.clear();
                    bin_buf.clear();
                }
                infile.read(reinterpret_cast<char*>(&n_no_coor_)
                    , sizeof(uint64_t));
                if (infile.peek() != EOF)
                {
                    std::cerr << "ERROR: BAI file doesn't read properly";
                    std::cerr << std::endl;
                    return false;
                }
                return true;
            };

            /**
             *  @brief Test whether index information load successfully.
             *  @return Whether information load successfully
             */
            bool is_useable() const
            {
                return !(index_info_.empty());
            };

            /**
             *  @brief Set interested reference region.
             *  @param r The region is interested
             *
             *  Set interested reference region
             *  and seek to nearest block address.<BR>
             *  Will check whether region is valid then use index information
             *  to seek to a alignment might overlapped with region.
             */
            void set_region(const RegionType& r)
            {
                if (!is_useable())
                {
                    std::cerr << "ERROR: Index information is empty.";
                    std::cerr << std::endl;
                    return;
                }
                int64_t left_ref =
                    (std::get<REGION_INDEX::LEFT_REF_ID>(r) <= -1) ?
                    0 : std::get<REGION_INDEX::LEFT_REF_ID>(r);
                int64_t left_pos =
                    (std::get<REGION_INDEX::LEFT_POS>(r) <= -1) ?
                    0 : std::get<REGION_INDEX::LEFT_POS>(r);
                int64_t right_ref =
                    (std::get<REGION_INDEX::RIGHT_REF_ID>(r) <= -1) ?
                    index_info_.size() - 1
                                                                    : std::get<REGION_INDEX::RIGHT_REF_ID>(r);
                int64_t right_pos =
                    (std::get<REGION_INDEX::RIGHT_POS>(r) <= -1) ?
                    LONGEST_REFERENCE_LENGTH
                                                                 : std::get<REGION_INDEX::RIGHT_POS>(r);
                if (left_ref >= index_info_.size())
                {
                    std::cerr << "ERROR: Left reference out of range";
                    std::cerr << std::endl;
                    return;
                }
                else if (right_ref >= index_info_.size())
                {
                    std::cerr << "ERROR: Right reference out of range";
                    std::cerr << std::endl;
                    return;
                }
                else if (left_ref > right_ref ||
                         (left_ref == right_ref && left_pos > right_pos))
                {
                    std::cerr << "ERROR: Left position ";
                    std::cerr << "is greater than right position";
                    std::cerr << std::endl;
                    return;
                }

                region_ = r;
                current_index_ = -1;
                size_t min_index;  // min_offset's index in lenear_offset
                uint64_t min_offset;  // linear offset that include beg
                // buffer region in the middle of processing
                int64_t current_left_pos, current_right_pos;
                // overall overlap bins in region
                std::vector<uint32_t> overlap_bins;
                BAI::BinType::iterator bins_iter;
                overlap_chunks_.clear();
                for (int64_t current_ref = left_ref;
                     current_ref <= right_ref; ++current_ref)
                {
                    BAI::BinType& bins =
                        std::get<BAI::REFERENCE_BIN_INFO_INDEX::BINS>(
                            index_info_[current_ref]);
                    std::vector<uint64_t>& linear_offset =
                        std::get<BAI::REFERENCE_BIN_INFO_INDEX::LINEAR_INDEX>(
                            index_info_[current_ref]);
                    current_left_pos =
                        (current_ref == left_ref) ? left_pos : 0;
                    current_right_pos =
                        (current_ref == right_ref) ?
                        right_pos : LONGEST_REFERENCE_LENGTH;
                    min_index =
                        (current_left_pos >> BAM_LINEAR_SHIFT) >
                        linear_offset.size() - 1 ?
                        linear_offset.size() - 1 :
                        current_left_pos >> BAM_LINEAR_SHIFT;
                    min_offset = linear_offset[min_index];
                    overlap_bins =
                        reg2bins(current_left_pos, current_right_pos);
                    for (size_t i = 0;i < overlap_bins.size();++i)
                    {
                        bins_iter = bins.find(overlap_bins[i]);
                        if (bins_iter == bins.end()) continue;
                        for (size_t j = 0;j < (bins_iter->second).size();++j)
                        {
                            if (std::get<BAI::CHUNK_INDEX::END>(
                                (bins_iter->second)[j]) > min_offset)
                            {
                                overlap_chunks_
                                    .emplace_back((bins_iter->second)[j]);
                            }
                        }
                    }
                }
                auto chunk_cmp = [](BAI::ChunkType a, BAI::ChunkType b)
                {
                    return (std::get<BAI::CHUNK_INDEX::BEG>(a) <
                            std::get<BAI::CHUNK_INDEX::BEG>(b));
                };
                std::sort(overlap_chunks_.begin()
                    , overlap_chunks_.end()
                    , chunk_cmp);
            };

            /**
             *  @brief Get current interested reference region.
             *  @return Interested reference region
             */
            RegionType get_region()
            {
                return region_;
            };

            /**
             *  @brief Get number of unplaced reads which record in BAI file.
             *  @return The number of unplaced reads
             */
            uint64_t get_unplaced_read_num()
            {
                return n_no_coor_;
            }

            /**
             *  @brief Jump to specific reference position.
             *  @param refID Reference ID which position located
             *  @param pos Offset from the beginning of reference
             *  @see set_region()
             *
             *  Equals to call set_region() which left bound is parameters
             *  and with opened right bound,
             *  acting like user jump to a position.
             */
            void jump(int64_t refID, int64_t pos)
            {
                set_region(RegionType {refID, pos, -1, -1});
            };

            /**
             * @brief Test whether current region has start position.
             * @return Whether current region has start position
             */
            bool is_start_region_specified()
            {
                return !(std::get<REGION_INDEX::LEFT_REF_ID>(region_) == -1 &&
                         std::get<REGION_INDEX::LEFT_POS>(region_) == -1);
            }

            /**
             * @brief Test whether current region has end position.
             * @return Whether current region has end position
             */
            bool is_end_region_specified()
            {
                return !(std::get<REGION_INDEX::RIGHT_REF_ID>(region_) == -1 &&
                         std::get<REGION_INDEX::RIGHT_POS>(region_) == -1);
            }

        private:
            /**
             * @brief Find overlapping bins with specified range.
             * @param beg Starting position
             * @param end Ending position
             * @return std::vector stores overlapping bin numbers
             *
             * Calculate the list of bins that may overlap with region
             * [beg,end) (zero-based).<BR>
             * Function declared as static is because there has no data
             * member being used. Everyone can access it without an object.
             */
            static std::vector<uint32_t> reg2bins(int64_t beg, int64_t end)
            {
                int i = 0, k;
                --end;
                std::vector<uint32_t> list;
                list.emplace_back(0);
                for (k =    1 + (beg >> 26); k <=    1 + (end >> 26); ++k)
                    list.emplace_back(k);
                for (k =    9 + (beg >> 23); k <=    9 + (end >> 23); ++k)
                    list.emplace_back(k);
                for (k =   73 + (beg >> 20); k <=   73 + (end >> 20); ++k)
                    list.emplace_back(k);
                for (k =  585 + (beg >> 17); k <=  585 + (end >> 17); ++k)
                    list.emplace_back(k);
                for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k)
                    list.emplace_back(k);
                return list;
            }

            /// Indices of ChunkType to access chunk's begin and end
            enum CHUNK_INDEX {BEG, END};
            /// Indices of ReferenceBinInfoType to access
            /// bins information ans linear index
            enum REFERENCE_BIN_INFO_INDEX {BINS, LINEAR_INDEX};
            /// Data structure to store chunk region
            using ChunkType = std::tuple<uint64_t, uint64_t>;
            /// Data structure to store all bin's chunk list in a reference
            using BinType = std::map<uint32_t, std::vector<ChunkType>>;
            /// Data structure to store linear index in a reference
            using LinearIndexType = std::vector<uint64_t>;
            /// Data structure to store summary in a reference
            /// (include bins and linear index)
            using ReferenceBinInfoType = std::tuple<BinType
                , LinearIndexType>;

            /// Bits need to shift to get linear index's index
            const static uint8_t BAM_LINEAR_SHIFT = 14;
            /// Record longest reference length that BAM could accept
            const static int64_t LONGEST_REFERENCE_LENGTH = 0x1fffffff;
            /// Record BAI file's magic string
            static constexpr char BAI_MAGIC[4] = {'B', 'A', 'I', '\1'};

            /// Records the number of unmapped reads
            uint64_t n_no_coor_;
            /// Records all reference binning and linear index information
            std::vector<ReferenceBinInfoType> index_info_;
            /// Records current interested reference range
            RegionType region_;
            /// Records index of overlap_chunks_ now reached
            int64_t current_index_;
            /// Stores chunks that might overlapped with region
            std::vector<ChunkType> overlap_chunks_;
        };

        /**
         * @brief Class stores A BAM alignment.
         * @warning A BAM object has a reference to bam::Header.
         * the header MUST NOT be distruct before this object destruct.
         *
         * The main class to maintain a BAM alignment.<BR>
         * It provides functions to get/set alignment information
         * and get alignment from BAM file format.
         */
        class BAM
        {
            friend class sam::SAM;
        public:
            /**
             * @brief Header constructor to reference to right header.
             * @param header Corresponding header with same BAM file
             * @see get_obj()
             *
             * Construct a BAM object which its header is referenced to
             * the header which is loaded from desired BAM file.<BR>
             * This time BAM object still has no information about alignments
             * until called get_obj(), but header must be intialized first.
             */
            BAM(Header& header)
                : header_       ( header )
                , has_data_     ( false  )
                , data_members_ ( BAMMemberType {
                    std::string()
                    , 0
                    , -1
                    , -1
                    , 255
                    , std::vector<CigarType>()
                    , -1
                    , -1
                    , 0
                    , std::string()
                    , std::string()
                    , std::vector<OptionalFieldType>()
                    , -1
                })
            {
                std::get<BAM_MEMBER_INDEX::QNAME>(data_members_).reserve(23);
                std::get<BAM_MEMBER_INDEX::CIGAR>(data_members_).reserve(6);
                std::get<BAM_MEMBER_INDEX::SEQ>(data_members_).reserve(35);
                std::get<BAM_MEMBER_INDEX::QUAL>(data_members_).reserve(35);
                std::get<BAM_MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_)
                    .reserve(15);
            };

            /**
             * @brief Convert constructor from SAM object.
             * @param header BAM header constructed by rhs' header
             * @param rhs Source SAM object which want to convert.
             * @warning Parameter header should be constructed by rhs' header
             *          or might mismatch between alignment and header.
             *
             * A convert constructor from SAM to BAM.<BR>
             * Because SAM is a subset of BAM, BAM needs additional
             * information to construct it. So needs a bam::Header
             * to construct it. Of course, header should available
             * until this object destruct.
             */
            BAM(Header& header, const sam::SAM& rhs)
                : header_   ( header )
                , has_data_ (  true  )
            {
                using namespace sam;
                std::get<BAM_MEMBER_INDEX::QNAME>(data_members_) =
                    std::get<SAM_MEMBER_INDEX::QNAME>(rhs.data_members_);
                std::get<BAM_MEMBER_INDEX::FLAG>(data_members_) =
                    std::get<SAM_MEMBER_INDEX::FLAG>(rhs.data_members_);
                // RNAME
                const std::string& ref_name =
                    std::get<SAM_MEMBER_INDEX::RNAME>(rhs.data_members_);
                std::vector<ReferenceType>& refs =
                    std::get<HEADER_INDEX::REFERENCE>(header_.header_);
                std::vector<ReferenceType>::iterator refs_it = refs.begin();
                int32_t& ref_id =
                    std::get<BAM_MEMBER_INDEX::REF_ID>(data_members_);
                if (ref_name == "*")
                    ref_id = -1;
                else
                {
                    while (refs_it != refs.end())
                    {
                        if (std::get<REFERENCE_INDEX::REFERENCE_NAME>(*refs_it)
                            == ref_name)
                            break;
                        std::advance(refs_it, 1);
                    }
                    ref_id = (refs_it == refs.end()) ?
                             -1 : std::distance(refs.begin(), refs_it);
                }
                int32_t& pos =
                    std::get<BAM_MEMBER_INDEX::POS>(data_members_);
                pos = std::get<SAM_MEMBER_INDEX::POS>(rhs.data_members_) - 1;
                std::get<BAM_MEMBER_INDEX::MAPQ>(data_members_) =
                    std::get<SAM_MEMBER_INDEX::MAPQ>(rhs.data_members_);
                // CIGAR
                std::vector<CigarType>& cigar_vec =
                    std::get<BAM_MEMBER_INDEX::CIGAR>(data_members_);
                std::string cigar_str =
                    std::get<SAM_MEMBER_INDEX::CIGAR>(rhs.data_members_);
                std::vector<OptionalFieldType>& bam_of =
                    std::get<BAM_MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_);
                uint32_t consume_ref_len = 0;
                cigar_vec = str2cigar_vec(cigar_str, consume_ref_len);
                if (cigar_vec.size() > MAX_CIGAR_OP_NUM)
                {
                    LightString value;
                    value.append('I');
                    size_t cigar_vec_size = cigar_vec.size();
                    Header::parse_char<int32_t>(&value.c_str[1], cigar_vec_size);
                    value.size += sizeof(int32_t);
                    for (size_t i = 0;i < cigar_vec_size;++i)
                    {
                        if (value.size + sizeof(uint32_t) > value.capacity)
                            value.extend();
                        Header::parse_char<uint32_t>(&value.c_str[value.size]
                            , std::get<CIGAR_INDEX::LENGTH>(cigar_vec[i]) << 4
                              | std::get<CIGAR_INDEX::OP>(cigar_vec[i]));
                        value.size += sizeof(uint32_t);
                    }
                    bam_of.emplace_back(std::array<char, 2>{'C', 'G'}, 'B'
                        , std::string(value.c_str, value.size));
                    cigar_vec.clear();
                    cigar_vec.emplace_back('S'
                        , std::get<SAM_MEMBER_INDEX::SEQ>(
                            rhs.data_members_).size());
                    cigar_vec.emplace_back('N', (ref_id == -1) ?
                                                0 :
                                                std::get<REFERENCE_INDEX::REFERENCE_LENGTH>(refs[ref_id]));
                }
                // RNEXT
                const std::string& nref_name =
                    std::get<SAM_MEMBER_INDEX::RNEXT>(rhs.data_members_);
                if (nref_name == "*")
                    std::get<BAM_MEMBER_INDEX::RNEXT_ID>(data_members_) = -1;
                else if (nref_name == "=")
                    std::get<BAM_MEMBER_INDEX::RNEXT_ID>(data_members_) = ref_id;
                else
                {
                    refs_it = refs.begin();
                    while (refs_it != refs.end())
                    {
                        if (std::get<REFERENCE_INDEX::REFERENCE_NAME>(*refs_it)
                            == nref_name)
                            break;
                        std::advance(refs_it, 1);
                    }
                    std::get<BAM_MEMBER_INDEX::RNEXT_ID>(data_members_) =
                        (refs_it == refs.end()) ?
                        -1 : std::distance(refs.begin(), refs_it);
                }
                std::get<BAM_MEMBER_INDEX::PNEXT>(data_members_) =
                    std::get<SAM_MEMBER_INDEX::PNEXT>(rhs.data_members_) - 1;
                std::get<BAM_MEMBER_INDEX::TLEN>(data_members_) =
                    std::get<SAM_MEMBER_INDEX::TLEN>(rhs.data_members_);
                std::string& seq =
                    std::get<BAM_MEMBER_INDEX::SEQ>(data_members_);
                seq = std::get<SAM_MEMBER_INDEX::SEQ>(rhs.data_members_);
                // QUAL
                std::string& qual =
                    std::get<BAM_MEMBER_INDEX::QUAL>(data_members_);
                const std::string& qual_33 =
                    std::get<SAM_MEMBER_INDEX::QUAL>(rhs.data_members_);
                if (qual_33 == "*")
                    qual = std::string(seq.size(), '\xFF');
                else
                {
                    qual.reserve(qual_33.size());
                    for (size_t i = 0;i < qual_33.size();++i)
                        qual.push_back(qual_33[i] - 33);
                }
                // Optional fields
                const std::vector<OptionalFieldType>& sam_of =
                    std::get<SAM_MEMBER_INDEX::OPTIONAL_FIELDS>(rhs.data_members_);
                char value_type;
                int64_t number;
                for (size_t i = 0;i < sam_of.size();++i)
                {
                    value_type =
                        std::get<OPTIONAL_FIELD_INDEX::VALUE_TYPE>(sam_of[i]);
                    const std::string& value =
                        std::get<OPTIONAL_FIELD_INDEX::VALUE>(sam_of[i]);
                    if (value_type == 'i')
                    {
                        boost::spirit::qi::parse(value.begin()
                            , value.end()
                            , boost::spirit::long_
                            , number);
                        if (number <= UCHAR_MAX && number > -1)
                            value_type = 'C';
                        else if(number <= SCHAR_MAX && number >= SCHAR_MIN)
                            value_type = 'c';
                        else if (number <= USHRT_MAX && number > -1)
                            value_type = 'S';
                        else if (number <= SHRT_MAX && number >= SHRT_MIN)
                            value_type = 's';
                        else if (number <= UINT_MAX && number > -1)
                            value_type = 'I';
                        else if (number <= INT_MAX && number >= INT_MIN)
                            value_type = 'i';
                        else
                            value_type = 'Z';
                    }
                    bam_of.emplace_back(
                        std::get<OPTIONAL_FIELD_INDEX::TAG>(sam_of[i])
                        , value_type
                        , pack_optional_field(value_type, value, number));
                }
                // bin
                std::get<BAM_MEMBER_INDEX::BIN>(data_members_) =
                    reg2bin(pos, pos + consume_ref_len);
            }

            BAM(const BAM& rhs) = default;
            BAM(BAM&& rhs) = default;
            ~BAM() = default;

            /**
             * @brief Get a specific alignment field.
             * @tparam n The field index from BAM_MEMBER_INDEX
             * @return Specific alignment field
             * @see BAM_MEMBER_INDEX
             */
            template <size_t n>
            const auto& get_member() const
            {
                return std::get<n>(data_members_);
            };

            /**
             * @brief Set a specific alignment field.
             * @tparam n The field index from BAM_MEMBER_INDEX
             * @see BAM_MEMBER_INDEX
             */
            template <size_t n>
            void set_member
                (const std::tuple_element_t<n, BAMMemberType>& rhs)
            {
                std::get<n>(data_members_) = rhs;
            };

            /**
             * @brief Get the header which this object reference to.
             * @return Header object which reference to
             */
            const Header& get_header() const
            {
                return header_;
            };

            /**
             * @brief Test this object contains valid alignment information.
             * @return Whether this object contains valid alignment information
             * @warning If this object contains invalid alignment information,
             *          then no assumptions can be made about information that
             *          get from get_member().
             * @see get_member()
             */
            bool is_valid() const
            {
                return has_data_;
            }

            /**
             * @brief Set this alignment information to valid.
             * @param valid The value user want to set to
             * @warning This function should use carefully or
             * might get invalid data leads to undefined behavior.
             */
            void set_valid(bool valid)
            {
                has_data_ = valid;
            }

            /**
             * @brief Convert BAM object to SAM string format.
             * @return String that has been converted
             * @warning If object is invalid will return empty string.
             * @see is_valid()
             */
            std::string to_string()
            {
                if (!has_data_)
                    return std::string();
                LightString tmp;
                to_string_impl<0>(tmp);
                return std::string(tmp.c_str, tmp.size);
            }

            /// Clear all data in this object.
            void reset()
            {
                has_data_ = false;
                std::get<BAM_MEMBER_INDEX::FLAG>(data_members_) = 0;
                std::get<BAM_MEMBER_INDEX::REF_ID>(data_members_) = -1;
                std::get<BAM_MEMBER_INDEX::POS>(data_members_) = -1;
                std::get<BAM_MEMBER_INDEX::MAPQ>(data_members_) = 255;
                std::get<BAM_MEMBER_INDEX::RNEXT_ID>(data_members_) = -1;
                std::get<BAM_MEMBER_INDEX::PNEXT>(data_members_) = -1;
                std::get<BAM_MEMBER_INDEX::TLEN>(data_members_) = -1;
                std::get<BAM_MEMBER_INDEX::BIN>(data_members_) = 0;
                std::get<BAM_MEMBER_INDEX::QNAME>(data_members_).clear();
                std::get<BAM_MEMBER_INDEX::CIGAR>(data_members_).clear();
                std::get<BAM_MEMBER_INDEX::SEQ>(data_members_).clear();
                std::get<BAM_MEMBER_INDEX::QUAL>(data_members_).clear();
                std::get<BAM_MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_)
                    .clear();
            }

            /**
             * @brief Read alignment information from BAM file.
             * @param in istream which contains the BAM file
             * @param obj The object which want data tpo put into
             * @return Parameter in for continuous data flow
             * @warning Parameter in should be the istream object
             * which has been used to parse into the header object,
             * which is referenced by parameter obj.
             *
             * Read an alignment information in BAM file.<BR>
             * After all field are loaded sucessfully, this object
             * will become valid state, or will marked as invalid.
             */
            static std::istream& get_obj(std::istream& in, BAM& obj)
            {
                if (obj.has_data_)
                    obj.brief_reset();
                // read alignment from file
                char* data = new char[sizeof(int32_t)];
                if (!obj.header_.read_byte_data(
                    in, data, sizeof(int32_t)))
                {
                    delete[] data;
                    return in;
                }
                int32_t block_size =
                    Header::convert_char<int32_t>(data);
                delete[] data;
                data = new char[block_size];
                obj.header_.read_byte_data(in, data, block_size);
                int32_t data_counter = 32;  // to TLEN field

                // directly parse simple data into data_members_
                std::get<BAM_MEMBER_INDEX::REF_ID>(obj.data_members_) =
                    Header::convert_char<int32_t>(&data[0]);
                std::get<BAM_MEMBER_INDEX::POS>(obj.data_members_) =
                    Header::convert_char<int32_t>(&data[4]);
                std::get<BAM_MEMBER_INDEX::MAPQ>(obj.data_members_) =
                    Header::convert_char<int32_t>(&data[9]);
                std::get<BAM_MEMBER_INDEX::BIN>(obj.data_members_) =
                    Header::convert_char<int32_t>(&data[10]);
                std::get<BAM_MEMBER_INDEX::FLAG>(obj.data_members_) =
                    Header::convert_char<uint16_t>(&data[14]);
                std::get<BAM_MEMBER_INDEX::RNEXT_ID>(obj.data_members_) =
                    Header::convert_char<int32_t>(&data[20]);
                std::get<BAM_MEMBER_INDEX::PNEXT>(obj.data_members_) =
                    Header::convert_char<int32_t>(&data[24]);
                std::get<BAM_MEMBER_INDEX::TLEN>(obj.data_members_) =
                    Header::convert_char<int32_t>(&data[28]);
                // read_name 8
                uint8_t l_read_name =
                    Header::convert_char<uint8_t>(&data[8]);
                std::get<BAM_MEMBER_INDEX::QNAME>(obj.data_members_) =
                    &data[data_counter];
                data_counter += l_read_name;
                // CIGAR
                uint16_t n_cigar_op =
                    Header::convert_char<uint16_t>(&data[12]);
                int2cigar(obj, &data[data_counter], n_cigar_op);
                data_counter += n_cigar_op * sizeof(uint32_t);
                // SEQ 16
                uint8_t int_seq;
                int32_t l_seq = Header::convert_char<int32_t>(&data[16]);
                LightString seq;
                for (size_t i = 0;i < (l_seq + 1) / 2;++i)
                {
                    int_seq =
                        Header::convert_char<uint8_t>(&data[data_counter]);
                    seq.append(SEQ_NUM_TO_CHAR[int_seq >> 4]);
                    seq.append(SEQ_NUM_TO_CHAR[int_seq & 0xf]);
                    data_counter += sizeof(uint8_t);
                }
                if (seq.c_str[seq.size - 1] == '=')
                    --seq.size;
                std::get<BAM_MEMBER_INDEX::SEQ>(obj.data_members_) =
                    std::string(seq.c_str, seq.size);
                // QUAL
                char* qual_string = new char[l_seq];
                memcpy(qual_string, &data[data_counter], l_seq);
                std::get<BAM_MEMBER_INDEX::QUAL>(obj.data_members_) =
                    std::string(qual_string, l_seq);
                data_counter += l_seq;
                // optional fields
                std::vector<OptionalFieldType>& of =
                    std::get<BAM_MEMBER_INDEX::OPTIONAL_FIELDS>(obj.data_members_);
                std::array<char, 2> tag;
                char val_type;
                while (data_counter != block_size)
                {
                    tag[0] = data[data_counter++];
                    tag[1] = data[data_counter++];
                    val_type = data[data_counter++];
                    switch (val_type)
                    {
                        case 'A':
                        case 'c':
                        case 'C':
                        case 's':
                        case 'S':
                        case 'i':
                        case 'I':
                        case 'f':
                        {
                            uint8_t size = type2Size(val_type);
                            of.emplace_back(tag
                                , val_type
                                , std::string(&data[data_counter]
                                    , size));
                            data_counter += size;
                            break;
                        }
                        case 'Z':
                        case 'H':
                        {
                            int32_t beg = data_counter;
                            while (data[data_counter++] != 0) ;
                            of.emplace_back(tag
                                , val_type
                                , std::string(&data[beg]
                                    , data_counter - beg));
                            break;
                        }
                        case 'B':
                        {
                            uint8_t size = type2Size(data[data_counter]);
                            int32_t count =
                                Header::convert_char<int32_t>(
                                    &data[data_counter + 1]);
                            of.emplace_back(tag
                                , val_type
                                , std::string(&data[data_counter]
                                    , size * count + 5));
                            data_counter += size * count + 5;
                            if (tag[0] == 'C' &&
                                tag[1] == 'G' &&
                                size == sizeof(uint32_t))
                                int2cigar(obj
                                    , &data[data_counter - size * count]
                                    , count);
                            break;
                        }
                    }
                }
                obj.has_data_ = true;
                delete[] data;
                delete[] qual_string;
                return in;
            };

            /**
             * @brief Read alignment information from BAM file with
             *        random access.
             * @param in istream which contains the BAM file
             * @param obj The object which want data tpo put into
             * @param bai The index information is related to the BAM file
             * @return Parameter in for continuous data flow
             * @warning Parameter in should be the istream object which
             *          has been used to parse into the header object,
             *          which is referenced by parameter obj.<BR>
             *          Parameter bai should contains the index
             *          information related with corresponding BAM file,
             *          or will lead to undefined behaviors.
             * @see get_obj(std::istream&, BAM&)
             *
             * Read an alignment information in the specific reference
             * region which is in the BAM file.<BR>
             * First function will check index information is usable
             * then start to find the alignment which is in the region.
             * After found alignment and all field are loaded sucessfully,
             * this object will become valid state, or will marked as invalid.
             */
            static std::istream& get_obj(std::istream& in, BAM& obj, BAI& bai)
            {
                if (bai.overlap_chunks_.empty())
                    obj.brief_reset();
                else if (!bai.is_useable())
                {
                    std::cerr << "ERROR: Index information is empty";
                    std::cerr << std::endl;
                }
                else
                {
                    // start to jump to first overlap chunk
                    if (bai.current_index_ == -1)
                    {
                        bai.current_index_ = 0;
                        obj.header_.seek(in
                            , std::get<BAI::CHUNK_INDEX::BEG>(
                                bai.overlap_chunks_[0]));
                    }
                        // file position is out of our traking, needs relocate
                    else if (obj.header_.tell() >
                             std::get<BAI::CHUNK_INDEX::END>(
                                 bai.overlap_chunks_[bai.current_index_]))
                    {
                        uint64_t current_pos = obj.header_.tell();
                        auto beg_it =
                            bai.overlap_chunks_.begin() + ++bai.current_index_;
                        auto end_it = bai.overlap_chunks_.end();
                        auto step = std::distance(beg_it, end_it) / 2;
                        do
                        {
                            if (current_pos >
                                std::get<BAI::CHUNK_INDEX::BEG>(
                                    *(beg_it + step)))
                            {
                                bai.current_index_ += step;
                                std::advance(beg_it, step);
                            }
                            else if (current_pos <
                                     std::get<BAI::CHUNK_INDEX::BEG>(
                                         *(beg_it + step)))
                                std::advance(end_it, -step);
                            else
                                break;
                            step = std::distance(beg_it, end_it) / 2;
                        } while (step != 0);
                        if (current_pos >=
                            std::get<BAI::CHUNK_INDEX::END>(*beg_it))
                        {
                            ++bai.current_index_;
                            obj.header_.seek(in
                                , std::get<BAI::CHUNK_INDEX::BEG>(
                                    *end_it));
                        }
                    }

                    ALIGNMENT_STATUS current_state;
                    do
                    {
                        if (obj.header_.tell() >=
                            std::get<BAI::CHUNK_INDEX::END>(
                                bai.overlap_chunks_[bai.current_index_])
                            &&
                            bai.current_index_ !=
                            bai.overlap_chunks_.size() - 1)
                        {
                            obj.header_.seek(in
                                , std::get<BAI::CHUNK_INDEX::BEG>(
                                    bai.overlap_chunks_
                                    [++bai.current_index_]));
                        }
                        get_obj(in, obj);
                        current_state = check_alignment_status(obj, bai);
                    } while (current_state == ALIGNMENT_STATUS::NO_OVERLAP);
                    if (current_state == ALIGNMENT_STATUS::OUT_RANGE)
                        obj.brief_reset();
                }
                return in;
            };

            /**
             * @brief Convert this alignment information to
             *        BAM binary format.
             * @return Converted binary data stream stored in std::string
             * @warning If object is invalid will return empty string.
             * @see is_valid()
             */
            std::string get_bamdata() const
            {
                char data[Header::CHUNK_SIZE + 4];
                uint32_t size = 0;
                if (!has_data_)
                    return std::string();

                size = 4;
                Header::parse_char<int32_t>(&data[size]
                    , std::get<BAM_MEMBER_INDEX::REF_ID>(
                        data_members_)
                    , &size);
                Header::parse_char<int32_t>(&data[size]
                    , std::get<BAM_MEMBER_INDEX::POS>(
                        data_members_)
                    , &size);
                const std::string& read_name =
                    std::get<BAM_MEMBER_INDEX::QNAME>(data_members_);
                Header::parse_char<uint8_t>(&data[size]
                    , read_name.size() + 1
                    , &size);
                Header::parse_char<uint8_t>(&data[size]
                    , std::get<BAM_MEMBER_INDEX::MAPQ>(
                        data_members_)
                    , &size);
                Header::parse_char<uint16_t>(&data[size]
                    , std::get<BAM_MEMBER_INDEX::BIN>(
                        data_members_)
                    , &size);
                const std::vector<CigarType>& cigar =
                    std::get<BAM_MEMBER_INDEX::CIGAR>(data_members_);
                Header::parse_char<uint16_t>(&data[size]
                    , cigar.size()
                    , &size);
                Header::parse_char<uint16_t>(&data[size]
                    , std::get<BAM_MEMBER_INDEX::FLAG>(
                        data_members_)
                    , &size);
                const std::string& seq =
                    std::get<BAM_MEMBER_INDEX::SEQ>(data_members_);
                Header::parse_char<int32_t>(&data[size]
                    , seq.size()
                    , &size);
                Header::parse_char<int32_t>(&data[size]
                    , std::get<BAM_MEMBER_INDEX::RNEXT_ID>(
                        data_members_)
                    , &size);
                Header::parse_char<int32_t>(&data[size]
                    , std::get<BAM_MEMBER_INDEX::PNEXT>(
                        data_members_)
                    , &size);
                Header::parse_char<int32_t>(&data[size]
                    , std::get<BAM_MEMBER_INDEX::TLEN>(
                        data_members_)
                    , &size);
                // read_name
                memcpy(&data[size], read_name.c_str(), read_name.size() + 1);
                size += read_name.size() + 1;
                // cigar
                uint8_t cigar_op;
                uint32_t cigar_len;
                for (size_t i = 0;i < cigar.size();++i)
                {
                    cigar_op =
                        get_index(CIGAR_NUM_TO_CHAR, sizeof(CIGAR_NUM_TO_CHAR)
                            , std::get<CIGAR_INDEX::OP>(cigar[i]));
                    cigar_len = std::get<CIGAR_INDEX::LENGTH>(cigar[i]);
                    Header::parse_char<uint32_t>(&data[size]
                        , cigar_len << 4 | cigar_op
                        , &size);
                }
                // seq
                uint8_t seq1, seq2;
                for (size_t i = 0;i < seq.size();i += 2)
                {
                    seq1 = get_index(SEQ_NUM_TO_CHAR,
                                     sizeof(SEQ_NUM_TO_CHAR), seq[i]);
                    seq2 = (i + 1 == seq.size()) ?
                           0 : get_index(SEQ_NUM_TO_CHAR,
                                         sizeof(SEQ_NUM_TO_CHAR), seq[i + 1]);
                    Header::parse_char<uint8_t>(&data[size]
                        , seq1 << 4 | seq2
                        , &size);
                }
                // qual
                const std::string& qual =
                    std::get<BAM_MEMBER_INDEX::QUAL>(data_members_);
                memcpy(&data[size], qual.c_str(), qual.size());
                size += qual.size();
                // optional fields
                const std::vector<OptionalFieldType>& of =
                    std::get<BAM_MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_);
                for (size_t i = 0;i < of.size();++i)
                {
                    data[size++] =
                        std::get<OPTIONAL_FIELD_INDEX::TAG>(of[i])[0];
                    data[size++] =
                        std::get<OPTIONAL_FIELD_INDEX::TAG>(of[i])[1];
                    data[size++] =
                        std::get<OPTIONAL_FIELD_INDEX::VALUE_TYPE>(of[i]);
                    const std::string& str =
                        std::get<OPTIONAL_FIELD_INDEX::VALUE>(of[i]);
                    memcpy(&data[size], str.c_str(), str.size());
                    size += str.size();
                }
                Header::parse_char<int32_t>(data, size - sizeof(int32_t));
                return std::string(data, size);
            }

            /**
             * @brief Convert list of BAM object to
             *        BGZF compressed format and output.
             * @param out ostream which data will output to
             * @param obj List of BAM object to convert
             * @see get_bamdata()
             *
             * Convert list of BAM object to BGZF compressed format
             * and output to the istream.<BR>
             * First convert each object to binary format then
             * compress these data to BGZF file format.
             */
            static void dump(std::ostream& out, std::vector<BAM>& obj)
            {
                if (!obj.empty())
                    out << obj[0].header_;

                char data[Header::CHUNK_SIZE];
                uint32_t size = 0;
                std::string buf;
                for (size_t i = 0;i < obj.size();++i)
                {
                    buf = obj[i].get_bamdata();
                    if (size + buf.size() > Header::CHUNK_SIZE)
                    {
                        Header::deflate_block(out, data, size);
                        size = 0;
                    }
                    memcpy(&data[size], buf.c_str(), buf.size());
                    size += buf.size();
                }
                if (size != 0)
                    Header::deflate_block(out, data, size);
                out.write(EOF_marker, sizeof(EOF_marker));
            };

            /// default copy assignment
            BAM& operator=(const BAM& rhs)
            {
                has_data_ = rhs.has_data_;
                header_ = rhs.header_;
                data_members_ = rhs.data_members_;
                return *this;
            }

            /// default move assignment
            BAM& operator=(BAM&& rhs)
            {
                has_data_ = rhs.has_data_;
                header_ = rhs.header_;
                data_members_ = std::move(rhs.data_members_);
                return *this;
            }

            /**
             * @brief overload bool() operator.
             * @return Whether alignment information is valid
             */
            operator bool() const
            {
                return has_data_;
            }

            /**
             * @brief overload >> operator.
             * @param in istream which data comes from
             * @param rhs BAM object which will contain alignment information
             * @return Parameter in for continuous data flow
             * @see get_obj(std::istream&, BAM&)
             *
             * Overload >> operator for user to use more intuitively.
             * Same as call get_obj() with sequencial access.
             */
            friend std::istream& operator>> (std::istream& in, BAM& rhs)
            {
                return get_obj(in, rhs);
            };

            /**
             * @brief overload << operator.
             * @param out ostream which data will output to
             * @param rhs BAM object which want to output
             * @return Parameter out for continuous data flow
             * @see get_bamdata()
             *
             * Output an alignment information in BGZF compressed format.<BR>
             * Same as call get_bamdata() and compress it to a block.
             */
            friend std::ostream& operator<<
                (std::ostream& out, const BAM& rhs)
            {
                std::string buf = rhs.get_bamdata();
                Header::deflate_block(out
                    , const_cast<char*>(buf.c_str())
                    , buf.size());
                return out;
            };

            /// BGZF file format's EOF marker
            constexpr static char EOF_marker[] =
                "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00";

        private:

            /// Only clear the field that is necessary
            void brief_reset()
            {
                has_data_ = false;
                std::get<BAM_MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_)
                    .resize(0);
            }

            /**
             * @brief Calculate value type to their corresponding size.
             * @param c The value type
             * @return Size of value type
             */
            constexpr static uint8_t type2Size(char c)
            {
                switch (c)
                {
                    case 'A':
                    case 'c':
                    case 'C':
                        return 1; // sizeof(char, int8_t, uint8_t)
                    case 's':
                    case 'S':
                        return 2; // sizeof(int16_t, uint16_t);
                    case 'i':
                    case 'I':
                    case 'f':
                        return 4; // sizeof(int32_t, uint32_t, float_t);
                    default:
                        return 0;
                }
            }

            /**
             * @brief Get target's index in character array.
             * @param arr The character array.
             * @param target Target element want to search.
             * @return If found, return index. Else return -1
             */
            uint8_t get_index(const char* arr, size_t size, char target) const
            {
                return std::distance(arr
                    , std::find(arr
                        , arr + size
                        , target));
            }

            /**
             * @brief Convert Cigar in BAM compressed format to CigarType.
             * @param obj BAM object which CigarType will put into
             * @param data Data stream pointer points to
             *             the start of Cigar data
             * @param count Number of Cigar operations
             */
            static void int2cigar(BAM& obj, char* data, int32_t count)
            {
                std::vector<CigarType>& cigar =
                    std::get<BAM_MEMBER_INDEX::CIGAR>(obj.data_members_);
                cigar.clear();
                uint32_t int_cigar;
                const uint32_t last_element = count * sizeof(uint32_t);
                for (size_t data_counter = 0;
                     data_counter < last_element;
                     data_counter += sizeof(uint32_t))
                {
                    int_cigar =
                        Header::convert_char<uint32_t>(&data[data_counter]);
                    cigar.emplace_back(CIGAR_NUM_TO_CHAR[int_cigar & 0xf]
                        , int_cigar >> 4);
                }
            }

            /// Enumerating different overlapping state
            /// between alignment and region.
            enum class ALIGNMENT_STATUS
            {
                NO_OVERLAP
                , OVERLAPPED
                , OUT_RANGE
            };

            /**
             * @brief Determine whether BAM alignment is
             *        overlapping with region.
             * @param obj BAM object which want to test
             * @param bai BAI object contains the interested region
             * @return Overlapping status between alignment and region
             */
            static ALIGNMENT_STATUS check_alignment_status(BAM& obj, BAI& bai)
            {
                if (!obj.has_data_)
                    return ALIGNMENT_STATUS::OUT_RANGE;
                int32_t bam_ref =
                    std::get<BAM_MEMBER_INDEX::REF_ID>(obj.data_members_);
                int32_t bam_start =
                    std::get<BAM_MEMBER_INDEX::POS>(obj.data_members_);
                int32_t bam_end =
                    bam_start + std::get<BAM_MEMBER_INDEX::QUAL>(
                        obj.data_members_).size();
                RegionType region = bai.get_region();

                if (bam_ref == -1)
                    return ALIGNMENT_STATUS::OUT_RANGE;
                else if (!bai.is_start_region_specified())
                    return ALIGNMENT_STATUS::OVERLAPPED;
                else if (bam_ref <
                         std::get<REGION_INDEX::LEFT_REF_ID>(region))
                    return ALIGNMENT_STATUS::NO_OVERLAP;
                else if (bam_ref ==
                         std::get<REGION_INDEX::LEFT_REF_ID>(region))
                {
                    if (bam_end <
                        std::get<REGION_INDEX::LEFT_POS>(region))
                        return ALIGNMENT_STATUS::NO_OVERLAP;
                    else if (bai.is_end_region_specified()
                             &&
                             bam_start >
                             std::get<REGION_INDEX::RIGHT_POS>(region)
                             &&
                             bam_ref ==
                             std::get<REGION_INDEX::RIGHT_REF_ID>(region))
                        return ALIGNMENT_STATUS::OUT_RANGE;
                    else
                        return ALIGNMENT_STATUS::OVERLAPPED;
                }
                else
                {
                    if (bai.is_end_region_specified())
                    {
                        if (bam_ref >
                            std::get<REGION_INDEX::RIGHT_REF_ID>(region))
                            return ALIGNMENT_STATUS::OUT_RANGE;
                        else if (bam_ref <
                                 std::get<REGION_INDEX::RIGHT_REF_ID>(region))
                            return ALIGNMENT_STATUS::OVERLAPPED;
                        else
                        {
                            if (bam_start >
                                std::get<REGION_INDEX::RIGHT_POS>(region))
                                return ALIGNMENT_STATUS::OUT_RANGE;
                            else
                                return ALIGNMENT_STATUS::OVERLAPPED;
                        }
                    }
                    else
                        return ALIGNMENT_STATUS::OVERLAPPED;
                }
            }

            /// Light string for internal use.
            struct LightString
            {
                LightString()
                    : c_str     ( new char[1024] )
                    , size      (       0        )
                    , capacity  (      1024      )
                {
                }

                LightString (const LightString& rhs)
                {
                    size = rhs.size;
                    capacity = rhs.capacity;
                    c_str = new char[capacity];
                    memcpy(c_str, rhs.c_str, size);
                }

                LightString (LightString&& rhs) = default;

                ~LightString()
                {
                    delete[] c_str;
                }

                void extend()
                {
                    char* tmp = new char[2 * capacity];
                    memcpy(tmp, c_str, size);
                    delete[] c_str;
                    c_str = tmp;
                    capacity = 2 * capacity;
                }

                void append(const char* c, size_t s)
                {
                    while (size + s > capacity)
                        extend();
                    memcpy(&c_str[size], c, s);
                    size += s;
                }

                void append(const std::string& str)
                {
                    append(str.c_str(), str.size());
                }

                void append(const char c)
                {
                    if (size + 1 > capacity)
                        extend();
                    c_str[size++] = c;
                }

                char* c_str;
                size_t size;
                size_t capacity;
            };

            /**
             * @brief According to value type to read value from data stream.
             * @param result String that result need to append to
             * @param c The value type
             * @param pos Data stream pointer points to the start of value
             * @param offset When value is an array,
             *               needs offset to help to locate right position
             */
            static void value_type_to_string
                (LightString& result, char c, const char* pos, size_t offset = 0)
            {
                size_t size;
                char tmp[45];
                switch (c)
                {
                    case 'c':
                        size = sprintf(tmp, "%d",
                                       Header::convert_char<int8_t>(
                                           pos + offset * sizeof(int8_t)));
                        break;
                    case 'C':
                        size = sprintf(tmp, "%u",
                                       Header::convert_char<uint8_t>(
                                           pos + offset * sizeof(uint8_t)));
                        break;
                    case 's':
                        size = sprintf(tmp, "%d",
                                       Header::convert_char<int16_t>(
                                           pos + offset * sizeof(int16_t)));
                        break;
                    case 'S':
                        size = sprintf(tmp, "%u",
                                       Header::convert_char<uint16_t>(
                                           pos + offset * sizeof(uint16_t)));
                        break;
                    case 'i':
                        size = sprintf(tmp, "%d",
                                       Header::convert_char<int32_t>(
                                           pos + offset * sizeof(int32_t)));
                        break;
                    case 'I':
                        size = sprintf(tmp, "%u",
                                       Header::convert_char<uint32_t>(
                                           pos + offset * sizeof(uint32_t)));
                        break;
                    case 'f':
                        size = sprintf(tmp, "%f",
                                       Header::convert_char<float>(
                                           pos + offset * sizeof(float)));
                        break;
                    default:
                        assert(false && "Invalid value type in optional field, file might broken.");
                }
                result.append(tmp, size);
            }

            /**
             * @brief Convert list of CigarType to
             *        Cigar string in SAM text format.
             * @param result String that result need to append to
             * @param vec List of CigarType to convert
             */
            static void cigar_vec2str
                (LightString& result, const std::vector<CigarType>& vec)
            {
                if (vec.size() == 0)
                    result.append('*');
                size_t size;
                char tmp[32];
                for (size_t i = 0;i < vec.size();++i)
                {
                    size = sprintf(tmp, "%u", std::get<CIGAR_INDEX::LENGTH>(vec[i]));
                    result.append(tmp, size);
                    result.append(std::get<CIGAR_INDEX::OP>(vec[i]));
                }
            }

            /**
             * @brief Convert Cigar string to list of CigarType.
             * @param cigar_str Target string to convert
             * @param consume_ref_len Reference counter for bin calculation
             * @return Result in vector of CigarType
             */
            static std::vector<CigarType> str2cigar_vec
                (std::string& cigar_str, uint32_t& consume_ref_len)
            {
                char cigar_op;
                uint32_t cigar_len, base;
                int64_t counter = cigar_str.size() - 1;
                std::vector<CigarType> result;
                result.reserve(5);
                if (cigar_str != "*")
                {
                    while (counter != -1)
                    {
                        cigar_op = cigar_str[counter--];
                        cigar_len = 0;
                        base = 1;
                        while (counter != -1 &&
                               cigar_str[counter] >= '0' && cigar_str[counter] <= '9')
                        {
                            cigar_len += (cigar_str[counter--] - '0') * base;
                            base *= 10;
                        }
                        if (std::find(CONSUME_REF_CIGAR
                            , CONSUME_REF_CIGAR + sizeof(CONSUME_REF_CIGAR)
                            , cigar_op) !=
                            CONSUME_REF_CIGAR + sizeof(CONSUME_REF_CIGAR))
                            consume_ref_len += cigar_len;
                        result.emplace(result.begin(), cigar_op, cigar_len);
                    }
                }
                return result;
            }

            /**
             * @brief Convert quality in BAM format to
             *        SAM visable text format.
             * @param result String that result need to append to
             * @param qual quality in BAM format
             */
            static void qual2ascii(LightString& result, const std::string& qual)
            {
                if (qual[0] == (char)255)
                {
                    result.append('*');
                    return;
                }

                size_t i = result.size;
                result.append(qual);
                for ( ;i < result.size;++i)
                    result.c_str[i] += 33;
            }

            /**
             * @brief Convert value in BAM binary format to SAM text format.
             * @param tmp String that result need to append to
             * @param value_type The value type of value
             * @param target Target binary format to convert
             * @see value_type_to_string()
             */
            static void unpack_optional_field_string
                (LightString& tmp, char value_type, const std::string& target)
            {
                switch (value_type)
                {
                    case 'A':
                        tmp.append(target.c_str(), target.size());
                        break;
                    case 'Z':
                        tmp.append(target.c_str(), target.size() - 1);
                        break;
                    case 'H':
                    {
                        const char* data_str = target.c_str();
                        const char hex_char[] = "0123456789ABCDEF";
                        for (size_t j = 1;data_str[j] != '\0';++j)
                        {
                            tmp.append(hex_char[data_str[j] >> 4]);
                            tmp.append(hex_char[data_str[j] & 0xf]);
                        }
                        break;
                    }
                    case 'B':
                    {
                        const char* data_str = target.c_str();
                        tmp.append(data_str[0]);
                        tmp.append(',');
                        int32_t count =
                            bam::Header::convert_char<int32_t>(&data_str[1]);
                        for (size_t j = 0;j < count;++j)
                        {
                            value_type_to_string(tmp, data_str[0]
                                , &data_str[5], j);
                            tmp.append(",", (j == count - 1 ? 0 : 1));
                        }
                        break;
                    }
                    default:
                        value_type_to_string(tmp, value_type, target.c_str());
                }
            }

            /**
             * @brief Convert value in SAM text format to BAM binary format.
             * @param c The value type
             * @param pos Data stream pointer points to
             *            the next byte to write
             * @param target Target value needs to convert
             * @param offset (Optional) Data stream index remainder,
             *               will pass to parse_char(). Default is nullptr
             * @see Header::parse_char()
             */
            void string_to_value_type(
                char c
                , char* pos
                , int64_t target
                , uint32_t* offset = nullptr
            )
            {
                switch (c)
                {
                    case 'c':
                        Header::parse_char<int8_t>(pos, target, offset);
                        return;
                    case 'C':
                        Header::parse_char<uint8_t>(pos, target, offset);
                        return;
                    case 's':
                        Header::parse_char<int16_t>(pos, target, offset);
                        return;
                    case 'S':
                        Header::parse_char<uint16_t>(pos, target, offset);
                        return;
                    case 'i':
                        Header::parse_char<int32_t>(pos, target, offset);
                        return;
                    case 'I':
                        Header::parse_char<uint32_t>(pos, target, offset);
                        return;
                }
            }

            /**
             * @brief Convert value in SAM text format to BAM binary format.
             * @param value_type The value type of value
             * @param target Target SAM text format to convert
             * @param converted_int Preparsed integer value
             * @return String of value in BAM binary format
             * @see string_to_value_type()
             */
            std::string pack_optional_field(
                char value_type
                , const std::string& target
                , int64_t converted_int
            )
            {
                char buffer[Header::CHUNK_SIZE];
                uint32_t buffer_size = 0;
                switch (value_type)
                {
                    case 'A':
                        buffer[buffer_size++] = target[0];
                        break;
                    case 'f':
                    {
                        float f;
                        boost::spirit::qi::parse(target.begin()
                            , target.end()
                            , boost::spirit::float_
                            , f);
                        Header::parse_char<float>(buffer
                            , f
                            , &buffer_size);
                        break;
                    }
                    case 'Z':
                    case 'H':
                        memcpy(&buffer, target.c_str(), target.size() + 1);
                        buffer_size += target.size() + 1;
                        break;
                    case 'B':
                    {
                        std::vector<std::string> split_vec;
                        boost:split(split_vec
                        , target
                        , boost::is_any_of(", ")
                        , boost::token_compress_on);
                        buffer[buffer_size++] = split_vec[0][0];
                        int32_t vec_size = split_vec.size();
                        Header::parse_char<int32_t>(&buffer[buffer_size]
                            , vec_size - 1
                            , &buffer_size);
                        for (size_t i = 1;i < vec_size;++i)
                        {
                            if (buffer[0] != 'f')
                                boost::spirit::qi::parse(split_vec[i].begin()
                                    , split_vec[i].end()
                                    , boost::spirit::long_
                                    , converted_int);
                            switch (buffer[0])
                            {
                                case 'f':
                                {
                                    float f;
                                    boost::spirit::qi::parse(split_vec[i].begin()
                                        , split_vec[i].end()
                                        , boost::spirit::float_
                                        , f);
                                    Header::parse_char<float>(&buffer[buffer_size]
                                        , f
                                        , &buffer_size);
                                    break;
                                }
                                default:
                                    string_to_value_type(buffer[0]
                                        , &buffer[buffer_size]
                                        , converted_int
                                        , &buffer_size);
                            }
                        }
                    }
                    default:
                        string_to_value_type(value_type
                            , buffer
                            , converted_int
                            , &buffer_size);
                }
                return std::string(buffer, buffer_size);
            }

            /**
             * @brief Check whether value type needs to convert to 'i' type.
             * @param value_type Original value type
             * @return Converted value type
             *
             * This function might be used when BAM convert to SAM,
             * because some value type in BAM is belongs to 'i' type in SAM.
             */
            static char check_to_i_type(char value_type)
            {
                if (std::find(VALUE_TYPE_TO_I
                    , VALUE_TYPE_TO_I + sizeof(VALUE_TYPE_TO_I)
                    , value_type) !=
                    VALUE_TYPE_TO_I + sizeof(VALUE_TYPE_TO_I))
                    return 'i';
                else
                    return value_type;
            }

            /**
             * @brief Convert optional fields to SAM text format.
             * @param tmp String that result need to append to
             * @param target List of OptionalFieldType needs to convert
             */
            void optional_fields_to_string
                (LightString& tmp, std::vector<OptionalFieldType>& target)
            {
                if (target.empty())
                    return;

                char value_type;
                for (size_t i = 0;i < target.size();++i)
                {
                    tmp.append(
                        std::get<OPTIONAL_FIELD_INDEX::TAG>(target[i])[0]);
                    tmp.append(
                        std::get<OPTIONAL_FIELD_INDEX::TAG>(target[i])[1]);
                    tmp.append(':');
                    value_type =
                        std::get<OPTIONAL_FIELD_INDEX::VALUE_TYPE>(target[i]);
                    tmp.append(check_to_i_type(value_type));
                    tmp.append(':');
                    unpack_optional_field_string(tmp, value_type
                        , std::get<
                            OPTIONAL_FIELD_INDEX::VALUE
                        >(target[i]));
                    tmp.append('\t');
                }
                tmp.size -= 1;
            }

            /**
             * @brief Implememtation details of to_string.
             * @tparam Current field index needs to convert
             * @param result String that result need to append to
             */
            template <size_t N>
            void to_string_impl(LightString& result)
            {
                std::tuple_element_t<N, BAMMemberType>& target =
                    std::get<N>(data_members_);
                if constexpr (N == BAM_MEMBER_INDEX::OPTIONAL_FIELDS)
                {
                    optional_fields_to_string(result, target);
                    if (result.c_str[result.size - 1] == '\t')
                        --result.size;
                }
                else
                {
                    if constexpr (N == BAM_MEMBER_INDEX::REF_ID)
                    {
                        if (target == -1)
                            result.append('*');
                        else
                            result.append(
                                std::get<REFERENCE_INDEX::REFERENCE_NAME>(
                                    header_.get_member<
                                        HEADER_INDEX::REFERENCE>()[target]));
                    }
                    else if constexpr (N == BAM_MEMBER_INDEX::RNEXT_ID)
                    {
                        if (target == -1)
                            result.append('*');
                        else if (target == std::get<
                            BAM_MEMBER_INDEX::REF_ID>(data_members_))
                            result.append('=');
                        else
                            result.append(std::get<
                                REFERENCE_INDEX::REFERENCE_NAME>(
                                header_.get_member<
                                    HEADER_INDEX::REFERENCE>()[target]));
                    }
                    else if constexpr (N == BAM_MEMBER_INDEX::CIGAR)
                        cigar_vec2str(result, target);
                    else if constexpr (N == BAM_MEMBER_INDEX::QUAL)
                        qual2ascii(result, target);
                    else if constexpr (std::is_same<
                        decltype(target), std::string&>())
                        result.append(target);
                    else
                    {
                        size_t size;
                        char tmp[35];
                        if constexpr (N == BAM_MEMBER_INDEX::POS ||
                                      N == BAM_MEMBER_INDEX::PNEXT)
                            size = sprintf(tmp, "%d", target + 1);
                        else
                            size = sprintf(tmp, "%d", target);
                        result.append(tmp, size);
                    }
                    result.append('\t');
                    to_string_impl<N + 1>(result);
                }
            }

            /**
             * @brief Calculate a region is belongs to which bin.
             * @param beg Starting position
             * @param end Ending position
             * @return Result bin number
             *
             * Calculate the list of bins that may overlap with region
             * [beg,end) (zero-based).
             */
            static uint16_t reg2bin(int32_t beg, int32_t end)
            {
                --end;
                if (beg >> 14 == end >> 14)
                    return ((1 << 15) - 1) / 7 + (beg >> 14);
                if (beg >> 17 == end >> 17)
                    return ((1 << 12) - 1) / 7 + (beg >> 17);
                if (beg >> 20 == end >> 20)
                    return ((1 <<  9) - 1) / 7 + (beg >> 20);
                if (beg >> 23 == end >> 23)
                    return ((1 <<  6) - 1) / 7 + (beg >> 23);
                if (beg >> 26 == end >> 26)
                    return ((1 <<  3) - 1) / 7 + (beg >> 26);
                return 0;
            }

            /// Record maximum number of Cigar operations
            /// that BAM could support.
            const static uint16_t MAX_CIGAR_OP_NUM = 65535;
            /// Record which value type in BAM should convert to 'i' in SAM.
            constexpr static char VALUE_TYPE_TO_I[] =
                {
                    'c', 'C', 's', 'S', 'i', 'I'
                };
            /// Record which Cigar operation will consumes reference position.
            constexpr static char CONSUME_REF_CIGAR[] =
                {
                    'M', 'D', 'N', '=', 'X'
                };
            /// A lookup table from Cigar operation code to character.
            constexpr static char CIGAR_NUM_TO_CHAR[] =
                {
                    'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'
                };
            /// A lookup table from sequence code to character.
            constexpr static char SEQ_NUM_TO_CHAR[] =
                {
                    '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                    'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
                };

            /// Record whether current alignment information is valid.
            bool has_data_;
            /// Reference to corresponding header.
            Header& header_;
            /// Store current alignment information.
            BAMMemberType data_members_;
        };
    };

    namespace sam
    {
        SAM::SAM(const bam::BAM& rhs)
            : header_   ( rhs.header_ )
        {
            using namespace bam;
            BAM::LightString ls;
            std::get<SAM_MEMBER_INDEX::QNAME>(data_members_) =
                std::get<BAM_MEMBER_INDEX::QNAME>(rhs.data_members_);
            std::get<SAM_MEMBER_INDEX::FLAG>(data_members_) =
                std::get<BAM_MEMBER_INDEX::FLAG>(rhs.data_members_);
            // RNAME
            int32_t ref_id =
                std::get<BAM_MEMBER_INDEX::REF_ID>(rhs.data_members_);
            std::vector<ReferenceType>& refs =
                std::get<HEADER_INDEX::REFERENCE>(header_.header_);
            std::get<SAM_MEMBER_INDEX::RNAME>(data_members_) =
                (ref_id == -1) ? "*" :
                std::get<REFERENCE_INDEX::REFERENCE_NAME>(refs[ref_id]);
            std::get<SAM_MEMBER_INDEX::POS>(data_members_) =
                std::get<BAM_MEMBER_INDEX::POS>(rhs.data_members_) + 1;
            std::get<SAM_MEMBER_INDEX::MAPQ>(data_members_) =
                std::get<BAM_MEMBER_INDEX::MAPQ>(rhs.data_members_);
            // CIGAR
            BAM::cigar_vec2str(ls,
                               std::get<BAM_MEMBER_INDEX::CIGAR>(rhs.data_members_));
            std::get<SAM_MEMBER_INDEX::CIGAR>(data_members_) =
                std::string(ls.c_str, ls.size);
            // RNEXT
            ref_id = std::get<BAM_MEMBER_INDEX::RNEXT_ID>(rhs.data_members_);
            std::get<SAM_MEMBER_INDEX::RNEXT>(data_members_) =
                (ref_id == -1) ? "*" :
                (ref_id == std::get<BAM_MEMBER_INDEX::REF_ID>(
                    rhs.data_members_)) ? "=" :
                std::get<REFERENCE_INDEX::REFERENCE_NAME>(
                    refs[ref_id]);
            std::get<SAM_MEMBER_INDEX::PNEXT>(data_members_) =
                std::get<BAM_MEMBER_INDEX::PNEXT>(rhs.data_members_) + 1;
            std::get<SAM_MEMBER_INDEX::TLEN>(data_members_) =
                std::get<BAM_MEMBER_INDEX::TLEN>(rhs.data_members_);
            std::get<SAM_MEMBER_INDEX::SEQ>(data_members_) =
                std::get<BAM_MEMBER_INDEX::SEQ>(rhs.data_members_);
            // QUAL
            ls.size = 0;
            BAM::qual2ascii(ls,
                            std::get<BAM_MEMBER_INDEX::QUAL>(rhs.data_members_));
            std::get<SAM_MEMBER_INDEX::QUAL>(data_members_) =
                std::string(ls.c_str, ls.size);
            // Optional fields
            std::vector<OptionalFieldType>& sam_of =
                std::get<SAM_MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_);
            const std::vector<OptionalFieldType>& bam_of =
                std::get<BAM_MEMBER_INDEX::OPTIONAL_FIELDS>(rhs.data_members_);
            char value_type;
            for (size_t i = 0;i < bam_of.size();++i)
            {
                value_type =
                    std::get<OPTIONAL_FIELD_INDEX::VALUE_TYPE>(bam_of[i]);
                ls.size = 0;
                BAM::unpack_optional_field_string(ls, value_type
                    , std::get<
                        OPTIONAL_FIELD_INDEX::VALUE
                    >(bam_of[i]));
                sam_of.emplace_back(
                    std::get<OPTIONAL_FIELD_INDEX::TAG>(bam_of[i])
                    , bam::BAM::check_to_i_type(value_type)
                    , std::string(ls.c_str, ls.size));
            }
        }
    };

};