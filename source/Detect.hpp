#ifndef DETECT_HPP
#define DETECT_HPP

#include "Writer.hpp"

void processInput(const std::string &inpFilePath, const int mean, const int stdDev, const int readLen, const double coverage, const std::string &outFolderPath, const bool verbose, const unsigned int threads);

#endif // DETECT_HPP
