#ifndef VARIANTDETECT_HPP
#define VARIANTDETECT_HPP

#include "OutputWriter.hpp"

void DelsParse(BamTools::BamReader &br,
			   const std::vector<Range> &dels_large,
			   const std::vector<Range> &dels_small,
			   std::vector<std::vector<std::string>> &info_dels_large,
			   std::vector<std::vector<std::string>> &info_dels_small,
			   std::vector<std::vector<std::string>> &info_ins,
			   std::vector<std::vector<int>> &pred_dels_large,
			   std::vector<std::vector<int>> &pred_dels_small,
			   std::vector<std::vector<int>> &pred_ins);

#endif // VARIANTDETECT_HPP