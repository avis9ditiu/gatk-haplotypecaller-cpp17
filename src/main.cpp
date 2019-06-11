#include <boost/program_options.hpp>
#include "haplotypecaller/haplotypecaller.hpp"

int main(int argc, char* argv[])
{
    using namespace boost::program_options;
    options_description desc("USAGE: HaplotypeCaller [arguments]");
    desc.add_options()
        ("input,I", value<std::string>(), "SAM file containing reads. Required.")
        ("output,O", value<std::string>(), "File to which variants should be written. Required.")
        ("reference,R", value<std::string>(), "Reference sequence file. Required.")
        ("help,h", "Display the help message");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help")) {
        std::cout << desc;
        return 0;
    }

    auto input = vm["input"].as<std::string>();
    auto output = vm["output"].as<std::string>();
    auto ref = vm["reference"].as<std::string>();
    
    hc::HaplotypeCaller{input, output, ref}.do_work();
}