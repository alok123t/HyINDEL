#include "ArgsParser.hpp"

#include <cassert>

int main()
{
    const std::string inpFile = "/Users/tmp/1.bam";
    bool isTrue = isValidExtension(inpFile);
    assert(isTrue == true);
    return EXIT_SUCCESS;
}
