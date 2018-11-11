#ifndef INPUTREADER_HPP
#define INPUTREADER_HPP

#include "VariantDetect.hpp"

const std::string indexExt = ".bai";

void processInput(const std::string filePath, const int mean, const int stdDev, const std::string folderPath);

#endif // INPUTREADER_HPP
