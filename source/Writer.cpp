#include "Writer.hpp"

const std::string suffixDelsLarge = "tmp/dels/large.txt";
const std::string suffixDelsSmall = "tmp/dels/small.txt";
const std::string suffixDelsLargeImprecise = "tmp/dels/large_imprecise.txt";
const std::string suffixIns = "tmp/ins/ins.txt";

std::string suffixFn(int outputType)
{
	if (outputType == 0)
		return suffixDelsSmall;
	else if (outputType == 1)
		return suffixDelsLarge;
	else if (outputType == 2)
		return suffixDelsLargeImprecise;
	else
		return std::string();
}

void headerOutput(const std::string &outFilePrefix, const int type)
{
	std::ofstream ofs;
	std::string outFile = outFilePrefix + suffixFn(type);
	ofs.open(outFile, std::ofstream::out);

	ofs << "#Chr" << '\t' << "Start" << '\t' << "End" << '\t' << "Size" << '\t' << "Support(PE)" << '\t' << "Support(SR)" << '\t' << "Support(SC)" << '\n';
	ofs.close();
}

void parseOutput(const std::string &outFilePrefix, const std::vector<OutNode> &output, const int type)
{
	std::ofstream ofs;
	std::string outFile = outFilePrefix + suffixFn(type);
	ofs.open(outFile, std::ofstream::out | std::ofstream::app);

	for (OutNode out : output)
	{
		ofs << out.chr << '\t' << out.st << '\t' << out.en << '\t' << out.en - out.st + 1 << '\t' << out.supDisc << '\t' << out.supSR << '\t' << out.supSC << '\n';
	}

	ofs.close();
}
