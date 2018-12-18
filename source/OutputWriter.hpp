#ifndef OUTPUTWRITER_HPP
#define OUTPUTWRITER_HPP

#include "Misc.hpp"

void parseOutput(const std::string outFilePrefix,
                 const std::vector<std::vector<std::string>> &output,
                 int outputType);

void parseOutputN(const std::string outFilePrefix, const std::vector<OutNode> &output);

#endif // OUTPUTWRITER_HPP
