#include "ArgsParser.hpp"

const std::string COMMA_DELIM = ",";

bool isValidExtension(const std::string & fileName) {
    std::size_t foundIdx = fileName.find_last_of(".");
    if (foundIdx != std::string::npos) {
        return fileName.substr(foundIdx) == ".bam";
    }
    else return false;
}

void parseFolder(std::string folderPath, std::vector<std::string> &filePaths) {
    if (folderPath.back() != '/') folderPath += "/";

    DIR *dir;
    struct dirent *ent;
    dir = opendir(folderPath.c_str());

    if (dir != NULL) {
        ent = readdir(dir);
        while (ent != NULL) {
            if (isValidExtension(ent->d_name)) {
                filePaths.push_back(folderPath + ent->d_name);
            }
            ent = readdir(dir);
        }
        closedir(dir);
    }
    else {
        std::cerr << "Input folder path does not exist\n";
    }
}

void splitFilePaths(const std::string mergedFileNames, std::vector<std::string> &filePaths) {
    std::size_t idx = mergedFileNames.find_first_not_of(COMMA_DELIM, 0);
    std::size_t commaIdx = mergedFileNames.find_first_of(COMMA_DELIM, idx);

    while (idx != std::string::npos || commaIdx != std::string::npos) {
        filePaths.push_back(mergedFileNames.substr(idx, commaIdx - idx));

        idx = mergedFileNames.find_first_not_of(COMMA_DELIM, commaIdx);
        commaIdx = mergedFileNames.find_first_of(COMMA_DELIM, idx);
    }
    return;
}

bool parseInput(int argc, char const *argv[], struct InputParams *ip) {
    args::ArgumentParser parser("This program identifies insertions and deletions in NGS data");
    args::HelpFlag help(parser, "help", "Help menu", {'h', "help"});

    args::Group groupInsertSize(parser, "Insert size parameters", args::Group::Validators::All);
    args::ValueFlag<int> insSz(groupInsertSize, "insSz", "Insert Size", {'s', "insSz"});
    args::ValueFlag<int> stdDev(groupInsertSize, "stdDev", "Standard Deviation", {'d', "stdDev"});

    args::Group groupInputFile(parser, "Input file name", args::Group::Validators::Xor);
    args::ValueFlag<std::string> inpFileName(groupInputFile, "inpFile", "Input FileName", {'i', "inFile"});
    args::ValueFlag<std::string> inpFolderPath(groupInputFile, "inpFolder", "Input FileName", {'f', "inFolder"});
    args::ValueFlag<std::string> inpFilesList(groupInputFile, "inpFiles", "Input FileName", {'m', "inFiles"});

    args::Group groupOutputFile(parser, "Output file name", args::Group::Validators::AllOrNone);
    args::ValueFlag<std::string> outFileName(groupOutputFile, "outFile", "Output FileName", {'o', "outFile"});

    args::Group groupVerboseFlag(parser, "Verbose", args::Group::Validators::AllOrNone);
    args::ValueFlag<int> verbose(groupVerboseFlag, "verbose", "Verbose", {'v', "verbose"});
    
    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
        std::cout << parser;
        return false;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return false;
    }
    catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return false;
    }

    (*ip).insSz = args::get(insSz);
    (*ip).stdDev = args::get(stdDev);
    
    if (inpFileName) (*ip).filePaths.push_back(args::get(inpFileName));
    else if (inpFolderPath) parseFolder(args::get(inpFolderPath), (*ip).filePaths);
    else if (inpFilesList) splitFilePaths(args::get(inpFilesList), (*ip).filePaths);

    if (verbose) (*ip).verbose = args::get(verbose);

    return true;
}
