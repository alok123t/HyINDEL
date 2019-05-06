#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <tuple>
#include <dirent.h>

// args api
#include "args.hxx"

const int matchScore = 1;
const int misMatchScore = -1;
const int gapPenalty = -1;
const double alnThreshold = 0.9;

bool parse(int argc, char const *argv[], std::string &outFolderPath)
{
    args::ArgumentParser parser("This program identifies insertions in NGS data");
    args::HelpFlag help(parser, "help", "Help menu", {'h', "help"});

    args::Group groupOutputFolder(parser, "Output Folder", args::Group::Validators::All);
    args::ValueFlag<std::string> outFolderName(groupOutputFolder, "outFolder", "Output Folder", {'o', "out"});

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

    outFolderPath = args::get(outFolderName);
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

std::tuple<int, std::string> align(const std::string &contigSeq, const std::string &scSeqUp, const std::string &scSeqDown)
{
    std::string revCompContigSeq;
    revComp(contigSeq, revCompContigSeq);

    std::tuple<bool, int, int> pUp1 = checkAlign(contigSeq, scSeqUp);
    std::tuple<bool, int, int> pUp2 = checkAlign(revCompContigSeq, scSeqUp);
    std::tuple<bool, int, int> pDown1 = checkAlign(contigSeq, scSeqDown);
    std::tuple<bool, int, int> pDown2 = checkAlign(revCompContigSeq, scSeqDown);

    // found up and down sc on one contig
    if (std::get<0>(pUp1) && std::get<0>(pDown1))
    {
        int st = std::get<1>(pUp1);
        int en = std::get<2>(pDown1);
        if (st >= contigSeq.size())
            return std::make_tuple(-1, std::string());
        std::string insSeq = contigSeq.substr(st, en - st);
        return std::make_tuple(0, insSeq);
    }
    else if (std::get<0>(pUp2) && std::get<0>(pDown2))
    {
        int st = std::get<1>(pUp2);
        int en = std::get<2>(pDown2);
        if (st >= revCompContigSeq.size())
            return std::make_tuple(-1, std::string());
        std::string insSeq = revCompContigSeq.substr(st, en - st);
        return std::make_tuple(0, insSeq);
    }
    // found up or down sc on one contig
    else if (std::get<0>(pUp1))
    {
        int st = std::get<1>(pUp1);
        if (st >= contigSeq.size())
            return std::make_tuple(-1, std::string());
        std::string insSeq = contigSeq.substr(st);
        return std::make_tuple(1, insSeq);
    }
    else if (std::get<0>(pUp2))
    {
        int st = std::get<1>(pUp2);
        if (st >= revCompContigSeq.size())
            return std::make_tuple(-1, std::string());
        std::string insSeq = revCompContigSeq.substr(st);
        return std::make_tuple(1, insSeq);
    }
    else if (std::get<0>(pDown1))
    {
        int en = std::get<2>(pDown1);
        std::string insSeq = contigSeq.substr(0, en);
        return std::make_tuple(2, insSeq);
    }
    else if (std::get<0>(pDown2))
    {
        int en = std::get<2>(pDown2);
        std::string insSeq = revCompContigSeq.substr(0, en);
        return std::make_tuple(2, insSeq);
    }
    return std::make_tuple(-1, std::string());
}

void readFile(std::ofstream &ofs, std::string &outFolderPath, std::string fName)
{
    std::size_t foundIdx = fName.find_last_of(".");
    std::string fPrefixName = fName.substr(0, foundIdx);
    std::string fContigsName = outFolderPath + "tmp/ins/" + fPrefixName + ".contigs.fa";
    std::string fSeqName = outFolderPath + "tmp/ins/" + fPrefixName + ".seq";

    std::ifstream ifs;
    ifs.open(fSeqName, std::ifstream::in);
    std::string scSeqUp, scSeqDown, refName;
    int pos, sup1, sup2;
    ifs >> scSeqUp >> scSeqDown >> refName >> pos >> sup1 >> sup2;
    for (int i = 0; i < scSeqUp.size(); i++)
        scSeqUp.at(i) = toupper(scSeqUp.at(i));
    for (int i = 0; i < scSeqDown.size(); i++)
        scSeqDown.at(i) = toupper(scSeqDown.at(i));
    ifs.close();

    ifs.open(fContigsName, std::ifstream::in);
    std::string contigHead, contigSeq;
    std::string insSeq, insUpSeq, insDownSeq;
    int found = -1;
    while (std::getline(ifs, contigHead))
    {
        std::getline(ifs, contigSeq);
        std::tuple<int, std::string> ret = align(contigSeq, scSeqUp, scSeqDown);
        int ty = std::get<0>(ret);
        std::string hereSeq = std::get<1>(ret);
        if (ty == 0)
        {
            found = 0;
            insSeq = hereSeq;
            break;
        }
        // imprecise
        else if (ty == 1)
        {
            found = 1;
            insUpSeq = hereSeq;
        }
        else if (ty == 2)
        {
            found = 1;
            insDownSeq = hereSeq;
        }
    }
    ifs.close();

    if (found == 0)
    {
        ofs << refName << '\t' << pos << '\t' << sup1 << '\t' << sup2 << '\t' << insSeq << '\n';
    }
    else if (found == 1)
    {
        ofs << refName << '\t' << pos << '\t' << sup1 << '\t' << sup2 << '\t' << insUpSeq << '\t' << insDownSeq << '\n';
    }
}

void traverseFiles(std::ofstream &ofs, std::string &outFolderPath)
{
    DIR *dir;
    struct dirent *ent;
    std::string tmpInsFolder = outFolderPath + "tmp/ins/";
    dir = opendir(tmpInsFolder.c_str());
    if (dir != NULL)
    {
        while ((ent = readdir(dir)) != NULL)
        {
            std::string fName = ent->d_name;
            if (isValidExtension(fName))
            {
                readFile(ofs, outFolderPath, fName);
            }
        }
        closedir(dir);
    }
}

int main(int argc, char const *argv[])
{
    std::string outFolderPath;
    if (!parse(argc, argv, outFolderPath))
        return EXIT_FAILURE;
    std::cerr << "[Step5] Align Insertions start\n";
    std::cerr << "Output folder: " << outFolderPath << '\n';

    std::ofstream ofs;
    std::string outFile = outFolderPath + "insertions.bed";
    ofs.open(outFile, std::ofstream::out);

    traverseFiles(ofs, outFolderPath);

    ofs.close();

    std::cerr << "[Step5] Align Insertions end\n";

    return EXIT_SUCCESS;
}
