#ifndef MISC_HPP
#define MISC_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

// bamtools api
#include "api/BamReader.h"
#include "api/BamWriter.h"

// task concurrency api
#include "transwarp.h"

const bool debug = false;

const std::string indexExt = ".bai";

const std::string excludeFileName = "remove_chr_windows.bed";

const std::string inpIntervalFileName = "chr_windows.bed";

const double RO_MERGE = 0.5;

const int MIN_MAP_QUAL = 20;

const int bpTol = 5;

const int gapPenalty = -2;
const int matchScore = 1;
const int misMatchScore = -2;
const int interAlignScore = 100; //150*2/3 * 1
const int interMinAlignScore = 40;
const int intraAlignScore = 60;
const int intraMinAlignScore = 20;

const int MIN_SC_CLUSTER_SUPPORT = 3;
const int MIN_SC_DISTANCE = 5;
const int MIN_SC_SUP = 3;
const int MIN_DEL_LEN = 45;
const int MAX_DEL_LEN = 505;

const int MIN_DISC_CLUSTER_SUPPORT = 3;
const int MAX_DISC_CLUSTER_DEL_LEN = 50000;
const int MIN_LARGE_DEL = 500;
const int MAX_LARGE_DEL_LEN = 20000;

const int MAX_SPLIT_LEN = 512345;
const int MIN_SPLIT_LEN = 51;

const std::string COMMA = ",";
const std::string SEMICOLON = ";";

struct OutNode
{
    std::string chr;
    int st, en;
    int supDisc, supSR, supSC;

    OutNode(std::string t_chr, int t_st, int t_en, int t_supDisc, int t_supSR, int t_supSC) : chr(t_chr), st(t_st), en(t_en), supDisc(t_supDisc), supSR(t_supSR), supSC(t_supSC) {}
};

struct SplitNode
{
    int st, en;
    int refID;

    SplitNode(int t_st, int t_en, int t_refID) : st(t_st), en(t_en), refID(t_refID) {}
};

struct SplitCluster
{
    SplitNode info;
    std::vector<SplitNode> nodes;

    SplitCluster(SplitNode t_sn) : info(t_sn) { nodes.emplace_back(t_sn); }
};

inline bool SplitCmp(const SplitCluster &a, const SplitCluster &b)
{
    return a.info.st < b.info.st;
}

struct DiscNode
{
    int upStart, upEnd, downStart, downEnd;
    int support;
    int refID;

    DiscNode(int t_upStart, int t_upEnd, int t_downStart, int t_downEnd, int t_refID, int t_support) : upStart(t_upStart), upEnd(t_upEnd), downStart(t_downStart), downEnd(t_downEnd), refID(t_refID), support(t_support) {}
};

struct DiscCluster
{
    DiscNode info;
    std::vector<DiscNode> nodes;

    DiscCluster(DiscNode t_dn) : info(t_dn) { nodes.emplace_back(t_dn); }
};

struct SoftNode
{
    int scPos, start, end;
    int refID;
    std::string refName;
    std::string seq, scSeq, nscSeq;
    bool down, scAtRight;

    SoftNode() : scPos() {}

    SoftNode(int t_scPos, int t_start, int t_end, int t_refID, std::string t_refName, std::string t_seq, std::string t_scSeq, std::string t_nscSeq, bool t_down, bool t_scAtRight) : scPos(t_scPos), start(t_start), end(t_end), refID(t_refID), refName(t_refName), seq(t_seq), scSeq(t_scSeq), nscSeq(t_nscSeq), down(t_down), scAtRight(t_scAtRight) {}
};

struct SoftCluster
{
    SoftNode info;
    std::vector<SoftNode> nodes;

    SoftCluster() : info() {}

    SoftCluster(SoftNode t_sn) : info(t_sn) { nodes.emplace_back(t_sn); }
};

inline bool SoftCmp(const SoftNode &a, const SoftNode &b)
{
    return a.scPos < b.scPos;
}

inline bool SoftCmpUp(const SoftNode &a, const SoftNode &b)
{
    return a.end < b.end;
}

inline bool SoftCmpDown(const SoftNode &a, const SoftNode &b)
{
    return a.start < b.start;
}

inline bool SoftClusterCmp(const SoftCluster &a, const SoftCluster &b)
{
    return a.info.scPos < b.info.scPos;
}

inline bool inBetween(const int &a, const int &b, const int &x)
{
    return a <= x && x <= b;
}

inline void getFileName(const std::string &filePath, std::string &outFilePath)
{
    std::string fileName = filePath.substr(filePath.find_last_of('/') + 1);
    std::size_t dotLen = fileName.find_last_of('.');
    std::string fileNameNoExt = fileName.substr(0, dotLen);
    outFilePath += fileNameNoExt;
}

inline bool openInput(const std::string filePath, BamTools::BamReader &br)
{
    if (!br.Open(filePath))
    {
        std::cerr << "ERROR:BamTools Could not open " << filePath << '\n';
        return false;
    }
    std::string indexPath = filePath + indexExt;
    if (!br.OpenIndex(indexPath))
    {
        std::cerr << "ERROR:BamTools Could not open index " << indexPath << '\n';
        return false;
    }
    return true;
}

inline void splitMultipleTag(const std::string &tag, std::vector<std::string> &tags)
{
    std::size_t idx = tag.find_first_not_of(SEMICOLON, 0);
    std::size_t semicolonIdx = tag.find_first_of(SEMICOLON, idx);

    while (idx != std::string::npos || semicolonIdx != std::string::npos)
    {
        tags.emplace_back(tag.substr(idx, semicolonIdx - idx));

        idx = tag.find_first_not_of(SEMICOLON, semicolonIdx);
        semicolonIdx = tag.find_first_of(SEMICOLON, idx);
    }
}

inline void splitTag(const std::vector<std::string> &tags, std::vector<std::vector<std::string>> &infoTags)
{
    for (std::string tag : tags)
    {
        std::vector<std::string> infoTag;
        std::size_t idx = tag.find_first_not_of(COMMA, 0);
        std::size_t commaIdx = tag.find_first_of(COMMA, idx);

        while (idx != std::string::npos || commaIdx != std::string::npos)
        {
            infoTag.emplace_back(tag.substr(idx, commaIdx - idx));

            idx = tag.find_first_not_of(COMMA, commaIdx);
            commaIdx = tag.find_first_of(COMMA, idx);
        }
        infoTags.emplace_back(infoTag);
    }
}

#endif // MISC_HPP