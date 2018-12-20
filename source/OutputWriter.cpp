#include "OutputWriter.hpp"

const std::string suffixDels = "_dels_large.txt";
const std::string suffixDelsSmall = "_dels_small.txt";
const std::string suffixDelsLargeSplit = "_dels_split.txt";

void parseOutput(const std::string outFilePrefix,
				 const std::vector<std::vector<std::string>> &output,
				 int outputType)
{
	std::ofstream ofs;
	auto suffixFn = [](int outputType) {
		if (outputType == 0)
			return suffixDels;
		else if (outputType == 1)
			return suffixDelsSmall;
		else if (outputType == 2)
			return suffixDelsLargeSplit;
		else
			return std::string();
	};
	std::string outFileDels = outFilePrefix + suffixFn(outputType);
	ofs.open(outFileDels);

	ofs << "#Chromosome" << '\t' << "Start" << '\t' << "End" << '\t' << "Size" << '\t' << "Support (PE)" << '\t' << "Support (SC)" << '\n';
	for (int i = 0; i < output.size(); ++i)
	{
		for (int j = 0; j < output.at(i).size(); ++j)
		{
			ofs << output.at(i).at(j);
			if (j != output.at(i).size() - 1)
				ofs << '\t';
		}
		ofs << '\n';
	}

	ofs.close();
}

void parseOutputN(const std::string outFilePrefix, const std::vector<OutNode> &output)
{
	std::ofstream ofs;
	std::string outFileDels = outFilePrefix + "_deletions.txt";
	ofs.open(outFileDels);

	ofs << "Chr" << '\t' << "Start" << '\t' << "End" << '\t' << "Size" << '\t' << "Support (PE)" << '\t' << "Support (SR)" << '\t' << "Support (SC)" << '\n';
	for (OutNode out : output)
	{
		ofs << out.chr << '\t' << out.st << '\t' << out.en << '\t' << out.sz << '\t' << out.supDisc << '\t' << out.supSR << '\t' << out.supSC << '\n';
	}

	ofs.close();
}
