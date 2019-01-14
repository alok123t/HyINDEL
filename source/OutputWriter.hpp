#ifndef OUTPUTWRITER_HPP
#define OUTPUTWRITER_HPP

#include "Misc.hpp"

void parseOutput(const std::string outFilePrefix,
                 const std::vector<std::vector<std::string>> &output,
                 int outputType);

void headerOutput(const std::string outFilePrefix, const int type);

void parseOutputN(const std::string outFilePrefix, const std::vector<OutNode> &output, const int type);

#endif // OUTPUTWRITER_HPP
