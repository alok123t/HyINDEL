#include "Parser.hpp"

const std::string COMMA_DELIM = ",";

bool isValidExtension(const std::string &fileName)
{
    std::size_t foundIdx = fileName.find_last_of(".");
    if (foundIdx != std::string::npos)
    {
        return fileName.substr(foundIdx) == ".bam";
    }
    else
    {
        return false;
    }
}

void checkFolder(const std::string &folderPath, std::string &outFolderPath)
{
    outFolderPath = folderPath;
    if (outFolderPath.back() != '/')
    {
        outFolderPath += "/";
    }

    DIR *dir;
    dir = opendir(outFolderPath.c_str());

    if (dir != NULL)
    {
        closedir(dir);
    }
}

bool parseArgs(int argc, char const *argv[], struct ArgsParams &ap)
{
    args::ArgumentParser parser("This program identifies insertions and deletions in NGS data");
    args::HelpFlag help(parser, "help", "Help menu", {'h', "help"});

    args::Group groupInputFileParams(parser, "Input file parameters", args::Group::Validators::All);
    args::ValueFlag<int> insSz(groupInputFileParams, "insSz", "Insert Size", {'s', "insSz"});
    args::ValueFlag<int> stdDev(groupInputFileParams, "stdDev", "Standard Deviation", {'d', "stdDev"});
    args::ValueFlag<int> readLen(groupInputFileParams, "readLen", "Read length", {'l', "readLen"});
    args::ValueFlag<double> cov(groupInputFileParams, "cov", "Coverage", {'c', "cov"});

    args::Group groupInputFiles(parser, "Input File", args::Group::Validators::All);
    args::ValueFlag<std::string> inpFilePath(groupInputFiles, "inpFile", "Input File", {'i', "inp"});

    args::Group groupOutputFolder(parser, "Output Folder", args::Group::Validators::All);
    args::ValueFlag<std::string> outFolderName(groupOutputFolder, "outFolder", "Output Folder", {'o', "out"});

    args::Group groupThreads(parser, "Threads", args::Group::Validators::AllOrNone);
    args::ValueFlag<unsigned int> threads(groupThreads, "threads", "Threads", {'t', "threads"});

    args::Group groupVerboseFlag(parser, "Verbose", args::Group::Validators::AllOrNone);
    args::ValueFlag<int> verbose(groupVerboseFlag, "verbose", "Verbose", {'v', "verbose"});

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

    ap.insSz = args::get(insSz);
    ap.stdDev = args::get(stdDev);
    ap.readLen = args::get(readLen);
    ap.cov = args::get(cov);

    if (inpFilePath)
    {
        ap.inpFilePath = args::get(inpFilePath);
    }

    if (outFolderName)
    {
        checkFolder(args::get(outFolderName), ap.outFolderPath);
    }

    ap.threads = std::thread::hardware_concurrency();
    if (threads)
    {
        ap.threads = std::max((unsigned int)1, std::min(args::get(threads), ap.threads));
    }

    if (verbose)
    {
        ap.verbose = args::get(verbose);
    }

    return true;
}
