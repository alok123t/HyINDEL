#ifndef ARGSPARSER_HPP
#define ARGSPARSER_HPP

#include <algorithm>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <dirent.h>

// args api
#include "args.hxx"

struct ArgsParams
{
    ArgsParams() : threads((int)std::thread::hardware_concurrency()), outFolderPath("./"), verbose(true) {}
    int insSz;
    int stdDev;
    std::vector<std::string> filePaths;
    std::string outFolderPath;
    unsigned int threads;
    bool verbose;
};

void splitFilePaths(const std::string mergedFileNames, std::vector<std::string> &filePaths);

bool parseArgs(int argc, char const *argv[], struct ArgsParams &ap);

#endif // ARGSPARSER_HPP
