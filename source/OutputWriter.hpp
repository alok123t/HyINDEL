#ifndef OUTPUTWRITER_HPP
#define OUTPUTWRITER_HPP

#include "Misc.hpp"

void parseOutput(const std::string outFilePrefix, \
	const std::vector<std::vector<std::string>> &info_dels_large, \
    const std::vector<std::vector<std::string>> &info_dels_small, \
    const std::vector<std::vector<std::string>> &info_ins, \
    const std::vector<std::vector<int>> &pred_dels_large, \
    const std::vector<std::vector<int>> &pred_dels_small, \
    const std::vector<std::vector<int>> &pred_ins);

#endif // OUTPUTWRITER_HPP
