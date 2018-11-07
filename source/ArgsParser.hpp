#ifndef ARGSPARSER_HPP
#define ARGSPARSER_HPP

// args api
#include "args.hxx"

#include <iostream>
#include <string>
#include <vector>
#include <dirent.h>

struct InputParams
{
    InputParams(): verbose(true) {}
    int insSz;
    int stdDev;
    std::vector<std::string> filePaths;
    bool verbose;
};

void splitFilePaths(const std::string mergedFileNames, std::vector<std::string> &filePaths);

bool parseInput(int argc, char const *argv[], struct InputParams *ip);

#endif // ARGSPARSER_HPP
