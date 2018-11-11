#include "OutputWriter.hpp"

const std::string suffixDels = "_dels.txt";

void parseOutput(const std::string outFilePrefix,
				 const std::vector<std::vector<std::string>> &info_dels_large,
				 const std::vector<std::vector<std::string>> &info_dels_small,
				 const std::vector<std::vector<std::string>> &info_ins,
				 const std::vector<std::vector<int>> &pred_dels_large,
				 const std::vector<std::vector<int>> &pred_dels_small,
				 const std::vector<std::vector<int>> &pred_ins)
{
	std::ofstream ofsDels;
	std::string outFileDels = outFilePrefix + suffixDels;
	ofsDels.open(outFileDels);

	ofsDels << "Chromosome" << '\t' << "Start" << '\t' << "End" << '\t' << "Size" << '\t' << "Support (PE)" << '\t' << "Support (SR)" << '\t' << "Breakpoint-Sequence" << '\n';
	for (int i = 0; i < pred_dels_large.size(); i++)
	{
		ofsDels << info_dels_large[i][0] << '\t';
		for (int j = 0; j < pred_dels_large[i].size(); j++)
		{
			ofsDels << pred_dels_large[i][j] << '\t';
		}
		ofsDels << info_dels_large[i][1] << info_dels_large[i][2] << '\n';
	}

	ofsDels.close();

	// //cout << "Small Deletions" << '\n';
	// for(int i = 0; i < pred_dels_small.size(); i++) {
	// 	std::cout << info_dels_small[i][0] << '\t';
	// 	for(int j = 0; j < pred_dels_small[i].size(); j++) {
	// 		std::cout << pred_dels_small[i][j] << '\t';
	// 	}
	// 	std::cout << info_dels_small[i][1] << info_dels_small[i][2];
	// 	std::cout << '\n';
	// }

	// //cerr << pred_ins.size() << ' ' << info_ins.size() << '\n';

	// //cout << "Insertions" << '\n';
	// for(int i = 0; i < pred_ins.size(); i++) {
	// 	std::cout << info_ins[i][0] << '\t';
	// 	for(int j = 0; j < pred_ins[i].size(); j++) {
	// 		std::cout << pred_ins[i][j] << '\t';
	// 	}
	// 	std::cout << info_ins[i][1] << info_ins[i][2];
	// 	std::cout << '\n';
	// }
}
