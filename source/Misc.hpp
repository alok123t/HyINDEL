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
const std::string inpIntervalFileName = "tmp/pre/chr_windows.bed";
const std::string excludeFileName = "tmp/pre/remove_chr_windows.bed";

const int MIN_MAP_QUAL = 20;

const int MATCH_SCORE = 1;
const int MIS_MATCH_SCORE = -1;
const int GAP_PENALTY = -1;
const double ALN_THRESHOLD = 0.9;

const int MAX_SC_OVERLAP_DISTANCE = 5;

const int MIN_SMALL_DEL_LEN = 45;
const int MAX_SMALL_DEL_LEN = 505;
const int MAX_LARGE_DEL_LEN = 50000;

const double DISC_RO = 0.65;

const double DEL_MERGE_RO = 0.95;
const int DEL_MERGE_DIS = 5;

const double DDEL_RO = 0.90;
const int MIN_DDEL_SUPPORT = 5;

const int MIN_INS_DIFF_CLUSTERS = -5;
const int MAX_INS_DIFF_CLUSTERS = 5;
const int MAX_INS_CLOSE = 10;
const int MIN_INS_SC_SIZE = 20;

const std::string COMMA = ",";
const std::string SEMICOLON = ";";

struct OutNode
{
    std::string chr;
    int st, en;
    int supDisc, supSR, supSC;

    OutNode(const std::string &t_chr, int t_st, int t_en, int t_supDisc, int t_supSR, int t_supSC) : chr(t_chr), st(t_st), en(t_en), supDisc(t_supDisc), supSR(t_supSR), supSC(t_supSC) {}
};

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

struct DirectDelNode
{
    int pos, delLen;
    std::string refName;

    DirectDelNode() : pos() {}

    DirectDelNode(int t_pos, int t_delLen, const std::string &t_refName) : pos(t_pos), delLen(t_delLen), refName(t_refName) {}
};

struct DirectDelCluster
{
    DirectDelNode info;
    std::vector<DirectDelNode> nodes;

    DirectDelCluster() : info() {}

    DirectDelCluster(DirectDelNode t_dn) : info(t_dn) { nodes.emplace_back(t_dn); }
};

struct SoftNode
{
    int scPos, start, end;
    int refID;
    std::string refName;
    std::string seq, readName;
    bool scAtRight, isSC;

    SoftNode() : scPos() {}

    SoftNode(int t_scPos, int t_start, int t_end, int t_refID, const std::string &t_refName, const std::string &t_seq, const std::string &t_readName, bool t_scAtRight, bool t_isSC) : scPos(t_scPos), start(t_start), end(t_end), refID(t_refID), refName(t_refName), seq(t_seq), readName(t_readName), scAtRight(t_scAtRight), isSC(t_isSC) {}
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
    return a.start < b.start;
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

inline bool cmpLen(const std::string &a, const std::string &b)
{
    return a.size() > b.size();
}

inline void getFileName(const std::string &filePath, std::string &outFilePath)
{
    std::string fileName = filePath.substr(filePath.find_last_of('/') + 1);
    std::size_t dotLen = fileName.find_last_of('.');
    std::string fileNameNoExt = fileName.substr(0, dotLen);
    outFilePath += fileNameNoExt;
}

inline bool openInput(const std::string &filePath, BamTools::BamReader &br)
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

const char REVCOMP_LOOKUP[] = {'T', 0, 'G', 'H', 0, 0, 'C', 'D', 0, 0, 0, 0, 'K',
                               'N', 0, 0, 0, 'Y', 'W', 'A', 'A', 'B', 'S', 'X', 'R', 0};

#endif // MISC_HPP