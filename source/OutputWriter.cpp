#include "OutputWriter.hpp"

const std::string suffixDels = "_dels.txt";
const std::string suffixDelsSmall = "_dels_small.txt";

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
		else
			return std::string();
	};
	std::string outFileDels = outFilePrefix + suffixFn(outputType);
	ofs.open(outFileDels);

	ofs << "Chromosome" << '\t' << "Start" << '\t' << "End" << '\t' << "Size" << '\t' << "Support (PE)" << '\t' << "Support (SR)" << '\n';
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
