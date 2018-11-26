#ifndef VARIANTDETECT_HPP
#define VARIANTDETECT_HPP

#include "OutputWriter.hpp"

int getDelSize(BamTools::BamAlignment &aln);

int trimSeq(std::string &scSeq, std::string &scQual, bool scAt5);

bool Align(const std::string &s1, const std::string &s2);

void DelsParse(std::string fPath, const std::vector<DiscCluster> &dels, std::vector<std::vector<std::string>> &output);

void DelsSmallMatch(std::map<int, std::vector<SoftCluster>> &cUp, std::map<int, std::vector<SoftCluster>> &cDown, std::vector<std::vector<std::string>> &output);

#endif // VARIANTDETECT_HPP