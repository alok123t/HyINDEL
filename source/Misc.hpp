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
#include <vector>

// bamtools api
#include "api/BamReader.h"
#include "api/BamWriter.h"

// task concurrency api
#include "transwarp.h"

const int minSCLen = 10;
const int bpTol = 5;
const int maxSc = 20;
const int maxCluster = 10;

// TODO: Change to 1 and check, do more analysis on what overlap criteria should be
const int matchScore = 2, mismatchScore = -2, gapPenalty = -2;
const double overlap = 0.8;

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

const int qualOffset = 33, minScQual = 10;

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