#ifndef INPUTREADER_HPP
#define INPUTREADER_HPP

#include "OutputWriter.hpp"

void processInput(const std::string inpFilePath, const int mean, const int stdDev, const std::string outFolderPath, const bool verbose, const unsigned int threads);

#endif // INPUTREADER_HPP
