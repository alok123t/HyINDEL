#include "OutputWriter.hpp"

const std::string suffixDels = "_dels.txt";
const std::string suffixDelsSmall = "_dels_small.txt";

void parseOutput(const std::string outFilePrefix,
				 const std::vector<std::vector<std::string>> &info_dels_large,
				 const std::vector<std::vector<std::string>> &info_dels_small)
{
	std::ofstream ofsDels;
	std::string outFileDels = outFilePrefix + suffixDels;
	ofsDels.open(outFileDels);

	ofsDels << "Chromosome" << '\t' << "Start" << '\t' << "End" << '\t' << "Size" << '\t' << "Support (PE)" << '\t' << "Support (SR)" << '\n';
	for (int i = 0; i < info_dels_large.size(); i++)
	{
		for (int j = 0; j < info_dels_large[i].size(); j++)
		{
			ofsDels << info_dels_large[i][j] << '\t';
		}
		ofsDels << '\n';
	}

	ofsDels.close();

	std::ofstream ofsDelsSmall;
	std::string outFileDelsSmall = outFilePrefix + suffixDelsSmall;
	ofsDelsSmall.open(outFileDelsSmall);

	ofsDelsSmall << "Chromosome" << '\t' << "Start" << '\t' << "End" << '\t' << "Size" << '\t' << "Support (PE)" << '\t' << "Support (SR)" << '\n';
	for (int i = 0; i < info_dels_small.size(); i++)
	{
		for (int j = 0; j < info_dels_small[i].size(); j++)
		{
			ofsDelsSmall << info_dels_small[i][j] << '\t';
		}
		ofsDelsSmall << '\n';
	}

	ofsDelsSmall.close();
}
