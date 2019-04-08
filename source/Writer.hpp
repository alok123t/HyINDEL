#ifndef WRITER_HPP
#define WRITER_HPP

#include "Misc.hpp"

void headerOutput(const std::string &outFilePrefix, const int type);

void parseOutput(const std::string &outFilePrefix, const std::vector<OutNode> &output, const int type);

#endif // WRITER_HPP
