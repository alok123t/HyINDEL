#include "ArgsParser.hpp"

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

void parseFolder(std::string folderPath, std::vector<std::string> &filePaths)
{
    if (folderPath.back() != '/')
    {
        folderPath += "/";
    }

    DIR *dir;
    struct dirent *ent;
    dir = opendir(folderPath.c_str());

    if (dir != NULL)
    {
        ent = readdir(dir);
        while (ent != NULL)
        {
            if (isValidExtension(ent->d_name))
            {
                filePaths.push_back(folderPath + ent->d_name);
            }
            ent = readdir(dir);
        }
        closedir(dir);
    }
    else
    {
        std::cerr << "Input folder path does not exist\n";
    }
}

void splitFilePaths(const std::string mergedFileNames, std::vector<std::string> &filePaths)
{
    std::size_t idx = mergedFileNames.find_first_not_of(COMMA_DELIM, 0);
    std::size_t commaIdx = mergedFileNames.find_first_of(COMMA_DELIM, idx);

    while (idx != std::string::npos || commaIdx != std::string::npos)
    {
        filePaths.push_back(mergedFileNames.substr(idx, commaIdx - idx));

        idx = mergedFileNames.find_first_not_of(COMMA_DELIM, commaIdx);
        commaIdx = mergedFileNames.find_first_of(COMMA_DELIM, idx);
    }
    return;
}

void checkFolder(const std::string &folderPath, std::string &outFolderPath)
{
    outFolderPath = folderPath;
    if (outFolderPath.back() != '/')
    {
        outFolderPath += "/";
    }

    DIR *dir;
    struct dirent *ent;
    dir = opendir(outFolderPath.c_str());

    if (dir != NULL)
    {
        closedir(dir);
    }
    else
    {
        outFolderPath = "./";
        std::cerr << "Output folder path does not exist\n";
        std::cerr << "Changing output folder path to \".\"";
    }
}

bool parseArgs(int argc, char const *argv[], struct ArgsParams &ap)
{
    args::ArgumentParser parser("This program identifies insertions and deletions in NGS data");
    args::HelpFlag help(parser, "help", "Help menu", {'h', "help"});

    args::Group groupInsertSize(parser, "Insert size parameters", args::Group::Validators::All);
    args::ValueFlag<int> insSz(groupInsertSize, "insSz", "Insert Size", {'s', "insSz"});
    args::ValueFlag<int> stdDev(groupInsertSize, "stdDev", "Standard Deviation", {'d', "stdDev"});

    args::Group groupFilteredInput(parser, "Discordant and Soft clip file paths", args::Group::Validators::AllOrNone);
    args::ValueFlag<std::string> inpFilterDiscPaths(groupFilteredInput, "discFiles", "Discordant File Paths", {"disc"});
    args::ValueFlag<std::string> inpFilterSCPaths(groupFilteredInput, "softFiles", "Softclip File Paths", {"soft"});

    args::Group groupOutputFolder(parser, "Output folder name", args::Group::Validators::AllOrNone);
    args::ValueFlag<std::string> outFolderName(groupOutputFolder, "outFolder", "Output Folder", {'o', "out"});

    args::Group groupThreads(parser, "Threads", args::Group::Validators::AllOrNone);
    args::ValueFlag<int> threads(groupThreads, "threads", "Threads", {'t', "threads"});

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

    if (inpFilterDiscPaths && inpFilterSCPaths)
    {
        splitFilePaths(args::get(inpFilterDiscPaths), ap.discFilePaths);
        splitFilePaths(args::get(inpFilterSCPaths), ap.softFilePaths);
    }

    if (outFolderName)
    {
        checkFolder(args::get(outFolderName), ap.outFolderPath);
    }

    if (threads)
    {
        ap.threads = std::max(1, std::min(args::get(threads), (int)(std::thread::hardware_concurrency())));
    }

    if (verbose)
    {
        ap.verbose = args::get(verbose);
    }

    return true;
}
