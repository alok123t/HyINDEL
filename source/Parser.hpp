#ifndef PARSER_HPP
#define PARSER_HPP

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
    ArgsParams() : threads((int)std::thread::hardware_concurrency()), verbose(true) {}
    int insSz;
    int stdDev;
    int readLen;
    double cov;
    std::string inpFilePath;
    std::string outFolderPath;
    unsigned int threads;
    bool verbose;
};

bool isValidExtension(const std::string &fileName);

bool parseArgs(int argc, char const *argv[], struct ArgsParams &ap);

#endif // PARSER_HPP
