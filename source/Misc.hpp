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

const int MIN_ALN_QUAL = 0;
const int MIN_MAP_QUAL = 20;

const int minSCLen = 10;
const int bpTol = 5;
const int maxSc = 20;
const int maxCluster = 10;

const int a_gap = -2;
const int a_matchScore = 1;
const int a_misMatchScore = -2;
const int interAlignScore = 100; //150*2/3 * 1
const int interMinAlignScore = 40;
const int intraAlignScore = 60;
const int intraMinAlignScore = 20;

const int MIN_SC_CLUSTER_SUPPORT = 3;
const int MIN_SC_DISTANCE = 5;

const int MIN_LARGE_DEL = 500;
const int MAX_LARGE_DEL_LEN = 20000;

// TODO: Change to 1 and check, do more analysis on what overlap criteria should be
const int matchScore = 2, mismatchScore = -2, gapPenalty = -2;
const double overlap = 0.8;

// Maximum chromosome size / 1000
const int MAX_SZ = 1123456;

const int BUCKET_SIZE = 1000;

// if cluster has less than below (up + down), it is rejected
const int MIN_SC_SUP = 3;

const int MAX_NODES_IN_CLUSTER = 500;
const int MAX_NODES_IN_BUCKET = 500;

const int MIN_DEL_LEN = 45;
const int MAX_DEL_LEN = 505;

const int MIN_DISC_CLUSTER_SUPPORT = 3;
const int MAX_DISC_CLUSTER_DEL_LEN = 50000;

const int MAX_SPLIT_LEN = 512345;
const int MIN_SPLIT_LEN = 51;

const std::string COMMA = ",";
const std::string SEMICOLON = ";";

struct OutNode
{
    std::string chr;
    int st, en, sz;
    int supDisc, supSR, supSC;
    OutNode(std::string t_chr, int t_st, int t_en, int t_sz, int t_supDisc, int t_supSR, int t_supSC) : chr(t_chr), st(t_st), en(t_en), sz(t_sz), supDisc(t_supDisc), supSR(t_supSR), supSC(t_supSC) {}
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

const int qualOffset = 33,
          minScQual = 10;

const std::map<char, int> qual2phred = {
    {'!', 33},
    {'\"', 34},
    {'#', 35},
    {'$', 36},
    {'%', 37},
    {'&', 38},
    {'\'', 39},
    {'(', 40},
    {')', 41},
    {'*', 42},
    {'+', 43},
    {',', 44},
    {'-', 45},
    {'.', 46},
    {'/', 47},
    {'0', 48},
    {'1', 49},
    {'2', 50},
    {'3', 51},
    {'4', 52},
    {'5', 53},
    {'6', 54},
    {'7', 55},
    {'8', 56},
    {'9', 57},
    {':', 58},
    {';', 59},
    {'<', 60},
    {'=', 61},
    {'>', 62},
    {'?', 63},
    {'@', 64},
    {'A', 65},
    {'B', 66},
    {'C', 67},
    {'D', 68},
    {'E', 69},
    {'F', 70},
    {'G', 71},
    {'H', 72},
    {'I', 73},
    {'J', 74},
    {'K', 75},
    {'L', 76},
    {'M', 77},
    {'N', 78},
    {'O', 79},
    {'P', 80},
    {'Q', 81},
    {'R', 82},
    {'S', 83},
    {'T', 84},
    {'U', 85},
    {'V', 86},
    {'W', 87},
    {'X', 88},
    {'Y', 89},
    {'Z', 90},
    {'[', 91},
    {'\\', 92},
    {']', 93},
    {'^', 94},
    {'_', 95},
    {'`', 96},
    {'a', 97},
    {'b', 98},
    {'c', 99},
    {'d', 100},
    {'e', 101},
    {'f', 102},
    {'g', 103},
    {'h', 104},
    {'i', 105},
    {'j', 106},
    {'k', 107},
    {'l', 108},
    {'m', 109},
    {'n', 110},
    {'o', 111},
    {'p', 112},
    {'q', 113},
    {'r', 114},
    {'s', 115},
    {'t', 116},
    {'u', 117},
    {'v', 118},
    {'w', 119},
    {'x', 120},
    {'y', 121},
    {'z', 122},
    {'{', 123},
    {'|', 124},
    {'}', 125},
    {'~', 126},
};

// TODO: class Range, add functions addRange ...
struct Range
{
    int refID1, refID2, start1, start2, end1, end2, support;
    bool isPlus1, isPlus2, isMinus1, isMinus2, isSC1, isSC2, isShort;
};

// TODO: Remove readName
struct Node
{
    std::string readName, seq, scSeq, nscSeq;
    int start, end, len, scLen, scPos;
    bool at5, scAt5, ins;
};

#endif // MISC_HPP