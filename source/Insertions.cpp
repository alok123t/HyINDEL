#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <tuple>
#include <dirent.h>

// args api
#include "args.hxx"

// task concurrency api
#include "transwarp.h"

const int matchScore = 1;
const int misMatchScore = -1;
const int gapPenalty = -1;
const double alnThreshold = 0.85;

int MIN_ORPHAN_SIZE;

// values: true, pre, mid, suf (is total length >= 1500)
// values: false, seq, "", "" (else)
std::vector<std::tuple<bool, std::string, std::string, std::string>> orphans;

bool parse(int argc, char const *argv[], int &readLen, std::string &outFolderPath, unsigned int &threads)
{
    args::ArgumentParser parser("This program identifies insertions in NGS data");
    args::HelpFlag help(parser, "help", "Help menu", {'h', "help"});

    args::Group groupOutputFolder(parser, "Output Folder", args::Group::Validators::All);
    args::ValueFlag<int> readLenVal(groupOutputFolder, "readLen", "Read length", {'l', "readLen"});
    args::ValueFlag<std::string> outFolderName(groupOutputFolder, "outFolder", "Output Folder", {'o', "out"});
    args::ValueFlag<unsigned int> threadsVal(groupOutputFolder, "threads", "Threads", {'t', "threads"});

    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        return false;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return false;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return false;
    }

    readLen = args::get(readLenVal);
    outFolderPath = args::get(outFolderName);
    threads = args::get(threadsVal);
    if (outFolderPath.back() != '/')
    {
        outFolderPath += "/";
    }
    std::string tmpInsFolder = outFolderPath + "tmp/ins/";

    DIR *dir;
    dir = opendir(outFolderPath.c_str());
    if (dir != NULL)
    {
        closedir(dir);

        DIR *dirIns;
        dirIns = opendir(tmpInsFolder.c_str());
        if (dirIns != NULL)
        {
            closedir(dirIns);
            return true;
        }
    }
    std::cerr << "Couldn't open output folder: " << outFolderPath << " or " << tmpInsFolder << '\n';
    return false;
}

bool isValidExtension(const std::string &fileName)
{
    std::size_t foundIdx = fileName.find_last_of(".");
    if (foundIdx != std::string::npos)
    {
        return fileName.substr(foundIdx) == ".fastq";
    }
    else
    {
        return false;
    }
}

const std::map<char, char> rMap = {
    {'A', 'T'},
    {'C', 'G'},
    {'G', 'C'},
    {'T', 'A'}};

void revComp(const std::string &inp, std::string &out)
{
    for (int i = inp.size() - 1; i >= 0; i--)
    {
        char c = toupper(inp.at(i));
        if (rMap.find(c) != rMap.end())
            out += rMap.at(c);
    }
}

// return args: found, start, end (positions wrt to large)
std::tuple<bool, int, int> checkAlign(const std::string &large, const std::string &small)
{
    int rows = large.size() + 1;
    int cols = small.size() + 1;

    std::vector<std::vector<int>> m(rows, std::vector<int>(cols));
    int maxScore = 0;
    for (int i = 0; i < rows; ++i)
        m.at(i).at(0) = 0;
    for (int j = 0; j < cols; ++j)
        m.at(0).at(j) = 0;
    int en;
    for (int i = 1; i < rows; ++i)
    {
        for (int j = 1; j < cols; ++j)
        {
            int sim = large.at(i - 1) == small.at(j - 1) ? matchScore : misMatchScore;
            m.at(i).at(j) = std::max(sim + m.at(i - 1).at(j - 1),
                                     std::max(gapPenalty + m.at(i - 1).at(j), gapPenalty + m.at(i).at(j - 1)));
            if (m.at(i).at(j) > maxScore)
            {
                maxScore = m.at(i).at(j);
                en = i;
            }
        }
    }
    if (maxScore >= alnThreshold * small.size())
        return std::make_tuple(true, en - small.size(), en);
    return std::make_tuple(false, -1, -1);
}

int median(std::vector<int> v)
{
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

// return args: found, median end position, orientation
std::tuple<bool, int, bool> alignSeqToContig(const std::vector<std::string> &seq, const std::string &contig)
{
    bool ok = true;
    std::vector<int> endPos;
    for (int i = 0; i < seq.size(); ++i)
    {
        std::tuple<bool, int, int> a = checkAlign(contig, seq.at(i));
        if (!std::get<0>(a))
        {
            ok = false;
            break;
        }
        endPos.emplace_back(std::get<2>(a));
    }
    if (ok)
        return std::make_tuple(true, median(endPos), true);

    ok = true;
    endPos.clear();
    for (int i = 0; i < seq.size(); ++i)
    {
        std::string revContig;
        revComp(contig, revContig);
        std::tuple<bool, int, int> a = checkAlign(revContig, seq.at(i));
        if (!std::get<0>(a))
        {
            ok = false;
            break;
        }
        endPos.emplace_back(std::get<2>(a));
    }
    if (ok)
        return std::make_tuple(true, median(endPos), false);
    return std::make_tuple(false, -1, false);
}

// return args: found, position
std::tuple<bool, int> alignContigToOrphan(const std::string &contig, const int orphanIdx, bool type)
{
    std::tuple<bool, std::string, std::string, std::string> curOrphan = orphans.at(orphanIdx);

    std::string wContig, wOrphan;
    // is upstream contig (find end pos by reversing)
    if (type)
    {
        wContig = contig;
        reverse(wContig.begin(), wContig.end());
    }
    // is downstream contig
    else
    {
        wContig = contig;
    }

    // large orphan
    if (std::get<0>(curOrphan))
    {
        if (type)
        {
            wOrphan = std::get<1>(curOrphan);
            reverse(wOrphan.begin(), wOrphan.end());
        }
        else
        {
            wOrphan = std::get<3>(curOrphan);
        }
    }
    // small orphan
    else
    {
        if (type)
        {
            wOrphan = std::get<1>(curOrphan);
            reverse(wOrphan.begin(), wOrphan.end());
        }
        else
        {
            wOrphan = std::get<1>(curOrphan);
        }
    }

    std::tuple<bool, int, int> a = checkAlign(wOrphan, wContig);
    if (std::get<0>(a))
    {
        return std::make_tuple(true, std::get<2>(a));
    }
    return std::make_tuple(false, -1);
}

void longInsertion(const std::string &upContig, const std::string &downContig, std::tuple<bool, std::string, std::string> &output)
{
    int orphanUpId = -1, orphanDownId = -1, orphanUpPos, reversedOrphanUpPos, orphanDownPos;
    // align upContig to Orphan contig
    for (int i = 0; i < orphans.size(); ++i)
    {
        std::tuple<bool, int> pUp = alignContigToOrphan(upContig, i, true);
        if (std::get<0>(pUp))
        {
            orphanUpId = i;
            orphanUpPos = std::get<1>(pUp);
            break;
        }
    }
    // align downContig to Orphan contig
    for (int i = 0; i < orphans.size(); ++i)
    {
        std::tuple<bool, int> pDown = alignContigToOrphan(downContig, i, false);
        if (std::get<0>(pDown))
        {
            orphanDownId = i;
            orphanDownPos = std::get<1>(pDown);
            break;
        }
    }
    // mapped both upContig and downContig onto same orphan
    if (orphanUpId == orphanDownId && orphanUpId != -1)
    {
        std::tuple<bool, std::string, std::string, std::string> curOrphan = orphans.at(orphanUpId);
        // large orphan
        if (std::get<0>(curOrphan))
        {
            reversedOrphanUpPos = 500 - orphanUpPos;
            std::string ret = std::get<1>(curOrphan).substr(reversedOrphanUpPos) +
                              std::get<2>(curOrphan) +
                              std::get<3>(curOrphan).substr(0, orphanDownPos);
            output = std::make_tuple(true, ret, std::string());
        }
        // small orphan
        else
        {
            reversedOrphanUpPos = std::get<1>(curOrphan).size() - orphanUpPos;
            int tillLen = orphanDownPos - reversedOrphanUpPos;
            std::string ret = std::get<1>(curOrphan).substr(reversedOrphanUpPos, tillLen);
            output = std::make_tuple(true, ret, std::string());
        }
        return;
    }

    std::string retPre, retSuf;
    if (orphanUpId != -1)
    {
        std::tuple<bool, std::string, std::string, std::string> curOrphan = orphans.at(orphanUpId);
        // large orphan
        if (std::get<0>(curOrphan))
        {
            reversedOrphanUpPos = 500 - orphanUpPos;
            retPre = std::get<1>(curOrphan).substr(reversedOrphanUpPos) +
                     std::get<2>(curOrphan) +
                     std::get<3>(curOrphan);
        }
        // small orphan
        else
        {
            reversedOrphanUpPos = std::get<1>(curOrphan).size() - orphanUpPos;
            retPre = std::get<1>(curOrphan).substr(reversedOrphanUpPos);
        }
    }
    if (orphanDownId != -1)
    {
        std::tuple<bool, std::string, std::string, std::string> curOrphan = orphans.at(orphanDownId);
        // large orphan
        if (std::get<0>(curOrphan))
        {
            retSuf = std::get<1>(curOrphan) +
                     std::get<2>(curOrphan) +
                     std::get<3>(curOrphan).substr(0, orphanDownPos);
        }
        // small orphan
        else
        {
            retSuf = std::get<1>(curOrphan).substr(orphanDownPos);
        }
    }
    output = std::make_tuple(false, retPre, retSuf);
}

void smallInsertion(const std::vector<std::string> &upSeq, const std::vector<std::string> &downSeq, const std::vector<std::string> &contigs, const std::string &refName, const int pos, const int upSup, const int downSup, const std::string &outFolderPath)
{
    int upId = -1, downId = -1, upPos, downPos;
    bool upTy = false, downTy = false;
    // align downSeq to Contig
    for (int i = 0; i < contigs.size(); ++i)
    {
        std::tuple<bool, int, bool> pDown = alignSeqToContig(downSeq, contigs.at(i));
        if (std::get<0>(pDown))
        {
            downId = i;
            downPos = std::get<1>(pDown);
            downTy = std::get<2>(pDown);
        }
    }
    // align upSeq (reversed) to Contig (reversed)
    std::vector<std::string> reversedUpSeq;
    for (int i = 0; i < upSeq.size(); ++i)
    {
        std::string tmp = upSeq.at(i);
        reverse(tmp.begin(), tmp.end());
        reversedUpSeq.emplace_back(tmp);
    }
    for (int i = 0; i < contigs.size(); ++i)
    {
        std::string reversedContig = contigs.at(i);
        reverse(reversedContig.begin(), reversedContig.end());
        std::tuple<bool, int, bool> pUp = alignSeqToContig(reversedUpSeq, reversedContig);
        if (std::get<0>(pUp))
        {
            upId = i;
            upPos = std::get<1>(pUp);
            upTy = std::get<2>(pUp);
        }
    }

    if (upId == downId && upId != -1)
    {
        std::ofstream ofs;
        std::string outFile = outFolderPath + "output.vcf";
        ofs.open(outFile, std::ofstream::out | std::ofstream::app);
        int reversedUpPos = contigs.at(upId).size() - upPos;
        int tillLen = downPos - reversedUpPos;
        if (upTy && downTy)
        {
            std::string out = contigs.at(upId).substr(reversedUpPos, tillLen);
            if (out.size() < 50)
            {
                ofs.close();
                return;
            }
            // ofs << "Small\t" << refName << '\t' << pos << '\t' << upSup << '\t' << downSup << '\t' << out << '\n';
            ofs << refName << '\t' << pos << '\t' << "0" << '\t' << "N" << '\t' << "<INS>" << '\t' << "." << '\t' << "PASS" << '\t'
                << "SVTYPE=INS;SVLEN=" << out.size() << ";SC=" << upSup + downSup << ";SEQ=" << out << ";" << '\t'
                << "GT:SU:SC" << '\t' << "./.:" << upSup + downSup << ":" << upSup + downSup << '\n';
        }
        else if (!upTy && !downTy)
        {
            std::string revContig;
            revComp(contigs.at(upId), revContig);
            std::string out = revContig.substr(reversedUpPos, tillLen);
            if (out.size() < 50)
            {
                ofs.close();
                return;
            }
            // ofs << "Small\t" << refName << '\t' << pos << '\t' << upSup << '\t' << downSup << '\t' << out << '\n';
            ofs << refName << '\t' << pos << '\t' << "0" << '\t' << "N" << '\t' << "<INS>" << '\t' << "." << '\t' << "PASS" << '\t'
                << "SVTYPE=INS;SVLEN=" << out.size() << ";SC=" << upSup + downSup << ";SEQ=" << out << ";" << '\t'
                << "GT:SU:SC" << '\t' << "./.:" << upSup + downSup << ":" << upSup + downSup << '\n';
        }
        ofs.close();
    }
    else
    {
        std::string downContig, upContig;
        if (downId != -1)
        {
            if (downTy)
            {
                downContig = contigs.at(downId).substr(0, downPos);
            }
            else
            {
                std::string revContig;
                revComp(contigs.at(downId), revContig);
                downContig = revContig.substr(0, downPos);
            }
        }
        if (upId != -1)
        {
            int reversedPos = contigs.at(upId).size() - upPos;
            if (upTy)
            {
                upContig = contigs.at(upId).substr(reversedPos);
            }
            else
            {
                std::string revContig;
                revComp(contigs.at(upId), revContig);
                upContig = revContig.substr(reversedPos);
            }
        }
        if (upContig.size() > 0 && downContig.size() > 0)
        {
            std::tuple<bool, std::string, std::string> largeOutput;
            longInsertion(upContig, downContig, largeOutput);
            std::ofstream ofs;
            std::string outFile = outFolderPath + "output.vcf";
            ofs.open(outFile, std::ofstream::out | std::ofstream::app);
            if (std::get<0>(largeOutput))
            {
                // ofs << "Large\t" << refName << '\t' << pos << '\t' << upSup << '\t' << downSup << '\t' << std::get<1>(largeOutput) << '\n'
                //     << std::flush;
                std::string outSeq = std::get<1>(largeOutput);
                ofs << refName << '\t' << pos << '\t' << "0" << '\t' << "N" << '\t' << "<INS>" << '\t' << "." << '\t' << "PASS" << '\t'
                    << "SVTYPE=INS;SVLEN=" << outSeq.size() << ";SC=" << upSup + downSup << ";SEQ=" << outSeq << ";" << '\t'
                    << "GT:SU:SC" << '\t' << "./.:" << upSup + downSup << ":" << upSup + downSup << '\n'
                    << std::flush;
            }
            else
            {
                // imprecise output
                // replace with larger orphanUpContig, if found
                if (std::get<1>(largeOutput).size() > upContig.size())
                    upContig = std::get<1>(largeOutput);
                // replace with larger orphanDownContig, if found
                if (std::get<2>(largeOutput).size() > downContig.size())
                    downContig = std::get<2>(largeOutput);
                // ofs << "Imprecise\t" << refName << '\t' << pos << '\t' << upSup << '\t' << downSup << '\t' << upContig << '\t' << downContig << '\n'
                //     << std::flush;
                if (upContig.size() < 50 || downContig.size() < 50)
                {
                    ofs.close();
                    return;
                }
                ofs << refName << '\t' << pos << '\t' << "0" << '\t' << "N" << '\t' << "<INS>" << '\t' << "." << '\t' << "PASS" << '\t'
                    << "SVTYPE=INS;SC=" << upSup + downSup << ";SEQUP=" << upContig << ";SEQDOWN=" << downContig << ";IMPRECISE" << '\t'
                    << "GT:SU:SC" << '\t' << "./.:" << upSup + downSup << ":" << upSup + downSup << '\n'
                    << std::flush;
            }
            ofs.close();
        }
    }
}

void readOrphans(const std::string &fOrphansName)
{
    std::ifstream ifs;
    std::string head;

    ifs.open(fOrphansName, std::ifstream::in);
    while (std::getline(ifs, head))
    {
        std::string seq;
        std::getline(ifs, seq);
        if (seq.size() >= 1500)
        {
            std::string pre = seq.substr(0, 500);
            std::string mid = seq.substr(500, (int)seq.size() - 1000);
            std::string suf = seq.substr((int)seq.size() - 500);
            orphans.emplace_back(std::make_tuple(true, pre, mid, suf));

            std::string preRC, midRC, sufRC;
            revComp(pre, preRC);
            revComp(mid, midRC);
            revComp(suf, sufRC);
            orphans.emplace_back(std::make_tuple(true, sufRC, midRC, preRC));
        }
        else
        {
            if (seq.size() < MIN_ORPHAN_SIZE)
            {
                continue;
            }
            orphans.emplace_back(std::make_tuple(false, seq, std::string(), std::string()));

            std::string seqRC;
            revComp(seq, seqRC);
            orphans.emplace_back(std::make_tuple(false, seqRC, std::string(), std::string()));
        }
    }
    ifs.close();
}
void readFiles(const std::string &fName, const std::string &outFolderPath)
{

    std::size_t foundIdx = fName.find_last_of(".");
    std::string fPrefixName = fName.substr(0, foundIdx);
    std::string fContigsName = outFolderPath + "tmp/ins/" + fPrefixName + ".contigs.fa";
    std::string fSeqName = outFolderPath + "tmp/ins/" + fPrefixName + ".seq";

    std::ifstream ifsSeq;
    int upCo, downCo, pos, upSup, downSup;
    std::vector<std::string> upSeq, downSeq;
    std::string refName;
    int minSeqLen = 1000000;

    ifsSeq.open(fSeqName, std::ifstream::in);
    ifsSeq >> upCo;
    for (int i = 0; i < upCo; ++i)
    {
        std::string seq;
        ifsSeq >> seq;
        upSeq.emplace_back(seq);
        minSeqLen = std::min(minSeqLen, (int)seq.size());
    }
    ifsSeq >> downCo;
    for (int i = 0; i < downCo; ++i)
    {
        std::string seq;
        ifsSeq >> seq;
        downSeq.emplace_back(seq);
        minSeqLen = std::min(minSeqLen, (int)seq.size());
    }
    ifsSeq >> refName >> pos >> upSup >> downSup;
    ifsSeq.close();

    if (minSeqLen <= 20)
    {
        return;
    }
    std::ifstream ifsContig;
    std::vector<std::string> contigs;
    std::string head;

    ifsContig.open(fContigsName, std::ifstream::in);
    while (std::getline(ifsContig, head))
    {
        std::string seq;
        std::getline(ifsContig, seq);
        contigs.emplace_back(seq);
    }
    ifsContig.close();

    smallInsertion(upSeq, downSeq, contigs, refName, pos, upSup, downSup, outFolderPath);
}

void traverseFiles(const std::string &outFolderPath, const unsigned int threads)
{
    DIR *dir;
    struct dirent *ent;
    std::string tmpInsFolder = outFolderPath + "tmp/ins/";
    dir = opendir(tmpInsFolder.c_str());

    transwarp::parallel executor{threads};
    std::vector<std::shared_ptr<transwarp::task<void>>> tasks;

    if (dir != NULL)
    {
        std::string fOrphansName = outFolderPath + "tmp/ins/orphans.contigs.fa";
        std::vector<std::string> orphans;
        readOrphans(fOrphansName);

        while ((ent = readdir(dir)) != NULL)
        {
            std::string fName = ent->d_name;
            if (isValidExtension(fName))
            {
                std::size_t orphansFound = fName.find("orphans");
                if (orphansFound != std::string::npos)
                    continue;
                // readFiles(fName, outFolderPath, ofs);
                auto outFolderTask = transwarp::make_value_task(outFolderPath);
                auto fNameTask = transwarp::make_value_task(fName);
                auto processTask = transwarp::make_task(transwarp::consume, readFiles, fNameTask, outFolderTask);
                tasks.emplace_back(processTask);
            }
        }
        closedir(dir);
    }

    for (auto task : tasks)
    {
        task->schedule(executor);
    }

    for (auto task : tasks)
    {
        task->get();
    }
}

int main(int argc, char const *argv[])
{
    int readLen;
    unsigned int threads;
    std::string outFolderPath;
    if (!parse(argc, argv, readLen, outFolderPath, threads))
        return EXIT_FAILURE;
    std::cerr << "[Step5] Align Insertions start\n";
    std::cerr << "        Read length: " << readLen << '\n';
    std::cerr << "        Output folder: " << outFolderPath << '\n';
    std::cerr << "        Threads: " << threads << '\n';

    MIN_ORPHAN_SIZE = 2 * readLen;

    traverseFiles(outFolderPath, threads);

    std::cerr << "[Step5] Align Insertions end\n";

    return EXIT_SUCCESS;
}
