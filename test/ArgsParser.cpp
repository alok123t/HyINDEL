#include "ArgsParser.hpp"

#include <cassert>

int main() {
    const std::string inpPaths = "/Users/tmp/1.bam,2.bam,/Users/3.bam";
    std::vector<std::string> filePaths;
    splitFilePaths(inpPaths, filePaths);
    std::vector<std::string> expPaths = {"/Users/tmp/1.bam", "2.bam", "/Users/3.bam"};
    assert(filePaths == expPaths);
    return EXIT_SUCCESS;
}
