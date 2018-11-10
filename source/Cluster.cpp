// TODO: Part of sub hpp
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <set>
#include <vector>
#include <chrono>

#include "transwarp.h"

#include "ArgsParser.hpp"
#include "InputReader.hpp"

// TODO: Remove namespaces
using namespace std;
using namespace BamTools;

set < int > dbg_refid;

/*




void DelsParse(BamReader &br) {
	set < int > sl_dels;
	auto prv = std::chrono::steady_clock::now();
	cerr << "Size SV Ranges : " << dels_large.size() << endl;
	int co = 0, edges_co = 0;
	for(int i = 0; i < dels_large.size(); i++) {
		// if ( i > 10 ) break;
		Range cur = dels_large[i];
		if(IsValid(cur)) {
			co++;
			// cerr << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << ' ' << cur.isSC1 << ' ' << cur.isSC2 << '\t' << cur.support << endl;
			co5 = co3 = 0;isIns = true;
			nodes.clear();clusters.clear();pred.clear();info.clear();

			// cerr << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << ' ';

			CreateNodes(cur, br);

			if(CreateLinks()) {
				// cerr << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << '\n';
				edges_co++;
			}

			// for (int cl = 0; cl < clusters.size(); cl++) cerr << "{" << clusters[cl].n_up.size() << ',' << clusters[cl].n_down.size() << "} ";
			// cerr << '\n';

			ClustersParse();

			if(isIns && pred.size()) {
				cerr << "large \n";
				pred[3] = cur.support;
				info[0] = br.GetReferenceData()[cur.refID1].RefName;
				info_ins.push_back(info);
				pred_ins.push_back(pred);
				continue;
			}

			if(pred.size()) {
				pred[3] = cur.support;
				// if(dbg_refid.find(cur.refID1) == dbg_refid.end()) {
				// 	cout << "cur refid : " << cur.refID1 << endl;
				// 	dbg_refid.insert(cur.refID1);
				// }
				info[0] = br.GetReferenceData()[cur.refID1].RefName;
				info_dels_large.push_back(info);
				pred_dels_large.push_back(pred);
			}
		}
		// else {
		// 	//cout << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << ' ' << cur.isSC1 << ' ' << cur.isSC2 << '\t' << cur.support << endl;
		// }
		int percentDone = (((i+1)*100/dels_large.size())/5)*5;
		if(percentDone != 0 && percentDone%5 == 0 && sl_dels.find(percentDone) == sl_dels.end()) {
			sl_dels.insert(percentDone);
			auto cur = std::chrono::steady_clock::now();
			// cerr << percentDone << "%\r" << std::flush;// << ' ' << setprecision(1) << (std::chrono::duration_cast<std::chrono::minutes>(cur - prv).count()) << endl;
			prv = cur;
		}
	}
	cerr << "Valid: " << co << '\n';
	cerr << "Edges built : " << edges_co << '\n';
	// cerr << "SVOutput count : " << svout_co << '\n';
	// cerr << "Edges failed count : " << fedg_co << '\n';
	// cerr << "Vertices failed count : " << fver_co << '\n';
	// cerr << "Final Deletions Count : " << out_co << '\n';

	// for(int i = 0; i < dels_small.size(); i++) {
	// 	Range cur = dels_small[i];
	// 	if(IsValid(cur)) {
	// 		co5 = co3 = 0;isIns = true;
	// 		nodes.clear();clusters.clear();pred.clear();info.clear();

	// 		CreateNodes(cur, br);

	// 		CreateLinks();

	// 		ClustersParse();

	// 		if(isIns && pred.size()) {
	// 			cerr << "small \n";
	// 			pred[3] = cur.support;
	// 			info[0] = br.GetReferenceData()[cur.refID1].RefName;
	// 			info_ins.push_back(info);
	// 			pred_ins.push_back(pred);
	// 			continue;
	// 		}

	// 		if(pred.size()) {
	// 			pred[3] = cur.support;
	// 			info[0] = br.GetReferenceData()[cur.refID1].RefName;
	// 			info_dels_small.push_back(info);
	// 			pred_dels_small.push_back(pred);
	// 		}
	// 	}
	// 	else {
	// 		//cout << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << ' ' << cur.isSC1 << ' ' << cur.isSC2 << '\t' << cur.support << endl;
	// 	}
	// }
}

void Output() {
	cout << "Chromosome" << '\t' << "Start" << '\t' << "End" << '\t' << "Size" << '\t' << "Support (PE)" << '\t' << "Support (SR)" << '\t' << "Breakpoint-Sequence" << endl;
	//cout << "Large Deletions" << endl;
	for(int i = 0; i < pred_dels_large.size(); i++) {
		cout << info_dels_large[i][0] << '\t';
		for(int j = 0; j < pred_dels_large[i].size(); j++) {
			cout << pred_dels_large[i][j] << '\t';
		}
		cout << info_dels_large[i][1] << info_dels_large[i][2];
		cout << endl;
		//cout << pred_dels_large[i][0] << '\t' << pred_dels_large[i][1] << '\t' << pred_dels_large[i][2] << endl;
	}

	//cout << "Small Deletions" << endl;
	for(int i = 0; i < pred_dels_small.size(); i++) {
		cout << info_dels_small[i][0] << '\t';
		for(int j = 0; j < pred_dels_small[i].size(); j++) {
			cout << pred_dels_small[i][j] << '\t';
		}
		cout << info_dels_small[i][1] << info_dels_small[i][2];
		cout << endl;
	}

	//cerr << pred_ins.size() << ' ' << info_ins.size() << endl;

	//cout << "Insertions" << endl;
	for(int i = 0; i < pred_ins.size(); i++) {
		cout << info_ins[i][0] << '\t';
		for(int j = 0; j < pred_ins[i].size(); j++) {
			cout << pred_ins[i][j] << '\t';
		}
		cout << info_ins[i][1] << info_ins[i][2];
		cout << endl;
	}
}
*/

int main(int argc, char const *argv[]) {

	std::chrono::steady_clock::time_point timeStart = std::chrono::steady_clock::now();

	ArgsParams ap;
	
	if (!parseArgs(argc, argv, &ap)) return EXIT_FAILURE;

	if (ap.verbose) {
		std::cerr << "Insert size: " << ap.insSz << '\n';
		std::cerr << "Standard deviation: " << ap.stdDev << '\n';

		std::cerr << "Input files: " << '\n';
		for (std::string fp: ap.filePaths) {
			std::cerr << '\t' << fp << '\n';
		}

		std::cerr << "Threads: " <<  ap.threads << '\n';
	}

	transwarp::parallel executor{ap.threads};

	std::vector<std::shared_ptr<transwarp::task<void>>> tasks;

	for (std::string fp: ap.filePaths) {
		auto fpTask = transwarp::make_value_task(fp);
		auto meanTask = transwarp::make_value_task(ap.insSz);
		auto stdDevTask = transwarp::make_value_task(ap.stdDev);
		auto inputTask = transwarp::make_task(transwarp::consume, processInput, fpTask, meanTask, stdDevTask);
		
		tasks.emplace_back(inputTask);
	}

	for (auto task: tasks) {
		task->schedule(executor);
	}

	for (auto task: tasks) {
		task->get();
	}

	std::chrono::steady_clock::time_point timeEnd = std::chrono::steady_clock::now();

	std::cerr << "Time taken: " << setprecision(1) << (std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeStart).count()) << " secs" << '\n';

	// DelsParse(br);

	// Output();

	return EXIT_SUCCESS;
}
