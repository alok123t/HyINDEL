#include <cmath>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <set>
#include <vector>
#include <chrono>
// args api
#include "args.hxx"
// bamtools api
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

using namespace std;
using namespace BamTools;

int Mean, StdDev;

int minSCLen = 10, bpTol = 5, maxSc = 20;
int maxCluster = 10;

// TODO: Change to 1 and check, do more analysis on what overlap criteria should be
int matchScore = 2, mismatchScore = -2, gapPenalty = -2;
double overlap = 0.8;

int qualOffset = 33, minScQual = 10;
map < char, int > qual2phred;

struct Range {
	int refID1, refID2, start1, start2, end1, end2, support;
	bool isPlus1, isPlus2, isMinus1, isMinus2, isSC1, isSC2, isShort;
};

// TODO: Remove readName
struct Node {
	string readName, seq, scSeq, nscSeq;
	int start, end, len, scLen, scPos;
	bool at5, scAt5, ins;
};

// TODO: Refactor struct into class, first have a constructor
struct Cluster {
	vector < Node > n_up, n_down;
};

vector < Range > dels_large, dels_small;

int co5, co3;

vector < Node > nodes;

vector < Cluster > clusters;

vector < string > info;
vector < int > pred;
vector < vector < string > > info_dels_large, info_dels_small, info_ins;
vector < vector < int > > pred_dels_large, pred_dels_small, pred_ins;

bool isIns;

void Init() {
	qual2phred['!'] = 33;
	qual2phred['\"'] = 34;
	qual2phred['#'] = 35;
	qual2phred['$'] = 36;
	qual2phred['%'] = 37;
	qual2phred['&'] = 38;
	qual2phred['\''] = 39;
	qual2phred['('] = 40;
	qual2phred[')'] = 41;
	qual2phred['*'] = 42;
	qual2phred['+'] = 43;
	qual2phred[','] = 44;
	qual2phred['-'] = 45;
	qual2phred['.'] = 46;
	qual2phred['/'] = 47;
	qual2phred['0'] = 48;
	qual2phred['1'] = 49;
	qual2phred['2'] = 50;
	qual2phred['3'] = 51;
	qual2phred['4'] = 52;
	qual2phred['5'] = 53;
	qual2phred['6'] = 54;
	qual2phred['7'] = 55;
	qual2phred['8'] = 56;
	qual2phred['9'] = 57;
	qual2phred[':'] = 58;
	qual2phred[';'] = 59;
	qual2phred['<'] = 60;
	qual2phred['='] = 61;
	qual2phred['>'] = 62;
	qual2phred['?'] = 63;
	qual2phred['@'] = 64;
	qual2phred['A'] = 65;
	qual2phred['B'] = 66;
	qual2phred['C'] = 67;
	qual2phred['D'] = 68;
	qual2phred['E'] = 69;
	qual2phred['F'] = 70;
	qual2phred['G'] = 71;
	qual2phred['H'] = 72;
	qual2phred['I'] = 73;
	qual2phred['J'] = 74;
	qual2phred['K'] = 75;
	qual2phred['L'] = 76;
	qual2phred['M'] = 77;
	qual2phred['N'] = 78;
	qual2phred['O'] = 79;
	qual2phred['P'] = 80;
	qual2phred['Q'] = 81;
	qual2phred['R'] = 82;
	qual2phred['S'] = 83;
	qual2phred['T'] = 84;
	qual2phred['U'] = 85;
	qual2phred['V'] = 86;
	qual2phred['W'] = 87;
	qual2phred['X'] = 88;
	qual2phred['Y'] = 89;
	qual2phred['Z'] = 90;
	qual2phred['['] = 91;
	qual2phred['\\'] = 92;
	qual2phred[']'] = 93;
	qual2phred['^'] = 94;
	qual2phred['_'] = 95;
	qual2phred['`'] = 96;
	qual2phred['a'] = 97;
	qual2phred['b'] = 98;
	qual2phred['c'] = 99;
	qual2phred['d'] = 100;
	qual2phred['e'] = 101;
	qual2phred['f'] = 102;
	qual2phred['g'] = 103;
	qual2phred['h'] = 104;
	qual2phred['i'] = 105;
	qual2phred['j'] = 106;
	qual2phred['k'] = 107;
	qual2phred['l'] = 108;
	qual2phred['m'] = 109;
	qual2phred['n'] = 110;
	qual2phred['o'] = 111;
	qual2phred['p'] = 112;
	qual2phred['q'] = 113;
	qual2phred['r'] = 114;
	qual2phred['s'] = 115;
	qual2phred['t'] = 116;
	qual2phred['u'] = 117;
	qual2phred['v'] = 118;
	qual2phred['w'] = 119;
	qual2phred['x'] = 120;
	qual2phred['y'] = 121;
	qual2phred['z'] = 122;
	qual2phred['{'] = 123;
	qual2phred['|'] = 124;
	qual2phred['}'] = 125;
	qual2phred['~'] = 126;
}

set < int > dbg_refid;

Range CreateRange(int refID1, int start1, int end1, bool isRev1, bool isSC1, int refID2, int start2, int end2, bool isRev2, bool isSC2, int support, bool isShort) {
	Range ret;

	// if(dbg_refid.find(refID1) == dbg_refid.end()) {
	// 	cout << "refid : " << refID1 << endl;
	// 	// dbg_refid.insert(refID1);
	// }

	ret.refID1 = refID1;ret.start1 = start1;ret.end1 = end1;ret.isPlus1 = !isRev1;ret.isMinus1 = isRev1;ret.isSC1 = isSC1;
	ret.refID2 = refID2;ret.start2 = start2;ret.end2 = end2;ret.isPlus2 = !isRev2;ret.isMinus2 = isRev2;ret.isSC2 = isSC2;
	ret.support = support;ret.isShort = isShort;
	return ret;
}

bool OverlapRange(Range r1, Range r2) {
	if(r1.refID1 == r2.refID1 && r1.refID2 == r2.refID2
		&& r1.start1 <= r2.end1 && r1.end1 >= r2.start1
		&& r1.start2 <= r2.end2 && r1.end2 >= r2.start2
		) return true;
	else return false;
}

bool OverlapOrientation(Range r1, Range r2) {
	if(r1.isPlus1 == r2.isPlus1 && r1.isPlus2 == r2.isPlus2
		&& r1.isMinus1 == r2.isMinus1 && r1.isMinus2 == r2.isMinus2
		) return true;
	else return false;
}

bool IsValid(Range r) {
	//Large SV
	if(!r.isShort) {
		//bp region shouldn't be less than 10bp
		if(r.isSC1 || r.isSC2 || abs(r.end1 - r.start1 + 1) < 10 || abs(r.end2 - r.start2 + 1) < 10)
			return false;
		else
			return true;
	}
	else {	//Small SV
		//Both SC need to be present
		if(r.isSC1 && r.isSC2)
			return true;
		else
			return false;
	}
}

void AddRange(vector < Range > & v, Range r) {
	if(v.size() == 0) {
		v.push_back(r);
		return;
	}

	vector < Range >::reverse_iterator rit = v.rbegin();
	while(rit != v.rend()) {
		//Diff chr
		if(rit->refID1 != r.refID1) break;
		//Doesn't overlap
		if((rit->end2 < r.start1) || (rit->start1 > r.end2)) break;

		if(OverlapRange(*rit, r) && OverlapOrientation(*rit, r)) {
			//Update
			if(r.start1 > rit->start1) rit->start1 = r.start1;
			if(r.start2 > rit->start2) rit->start2 = r.start2;
			if(r.end1 < rit->end1) rit->end1 = r.end1;
			if(r.end2 < rit->end2) rit->end2 = r.end2;

			rit->isPlus1 |= r.isPlus1;
			rit->isPlus2 |= r.isPlus2;
			rit->isMinus1 |= r.isMinus1;
			rit->isMinus2 |= r.isMinus2;
			rit->support += r.support;
			//Large
			if(!rit->isShort) {
				rit->isSC1 &= r.isSC1;
				rit->isSC2 &= r.isSC2;
			}//Small
			else {
				rit->isSC1 |= r.isSC1;
				rit->isSC2 |= r.isSC2;
			}
			return;
		}
		rit++;
	}
	v.push_back(r);
	return;
}

void Input(string arg, BamReader &br) {
	string inpBamFile = arg;
	string inpBamIdxFile = inpBamFile + ".bai";

	br.Open(inpBamFile);
	br.OpenIndex(inpBamIdxFile);
}

void ProcessBam(BamReader &br) {
	BamAlignment aln;
	bool aln_exists;

	int refID1, refID2, start1, start2, end1, end2, insSize, support, readLength, centre, bound1, bound2, Threshold, bpRegion;
	bool isRev1, isRev2, isSC1, isSC2, isShort, isUp, at5;
	vector <int> clipSizes, readPositions, genomePositions;
	vector <int> down(5), up(5);
	Range tmpRange;


	aln_exists = br.GetNextAlignmentCore(aln);
	while(aln_exists) {
		//Parse
		if(aln.IsMapped() && aln.IsMateMapped()) {
			clipSizes.clear();readPositions.clear();genomePositions.clear();
			//Large SV
			if(!aln.IsProperPair()) {
				//First Read
				readLength = aln.Length;
				isShort = false;
				Threshold = bpRegion = Mean + (3*StdDev);
				if(aln.Position < aln.MatePosition) {
					refID1 = aln.RefID;
					start1 = aln.Position;
					end1 = aln.GetEndPosition();
					isRev1 = aln.IsReverseStrand();
					isSC1 = aln.GetSoftClips(clipSizes, readPositions, genomePositions);

					refID2 = aln.MateRefID;
					start2 = aln.MatePosition;
					end2 = start2 + readLength;
					isRev2 = aln.IsMateReverseStrand();
					isSC2 = true;

					insSize = aln.InsertSize;
					support = 1;
				}//Second Read
				else {
					refID1 = aln.MateRefID;
					start1 = aln.MatePosition;
					end1 = start1 + readLength;
					isRev1 = aln.IsMateReverseStrand();
					isSC1 = true;

					refID2 = aln.RefID;
					start2 = aln.Position;
					end2 = aln.GetEndPosition();
					isRev2 = aln.IsReverseStrand();
					isSC2 = aln.GetSoftClips(clipSizes, readPositions, genomePositions);

					insSize = aln.InsertSize * -1;
					support = 0;
				}
				//Not a SC read
				if(!(isSC1 && isSC2)) {
					if(refID1 == refID2 && insSize >= Threshold && !isRev1 && isRev2) {
						centre = ceil((start2 + end1)/2);
						bound1 = min(centre, end1 + bpRegion);
						bound2 = max(centre, start2 - bpRegion);
						bound1 = end1 + bpRegion;
						bound2 = start2 - bpRegion;

						tmpRange = CreateRange(refID1, start1, bound1, isRev1, isSC1, refID2, bound2, end2, isRev2, isSC2, support, isShort);
						AddRange(dels_large, tmpRange);
					}
				}
			}
			//Small SV
			if(aln.IsProperPair()) {
				if(aln.RefID == aln.MateRefID && aln.GetSoftClips(clipSizes, readPositions, genomePositions) && clipSizes.size() == 1 && clipSizes[0] >= minSCLen) {
					bpRegion = 4*StdDev;
					isShort = true;
					support = 0;

					at5 = (aln.CigarData[0].Type == 'S');
					isUp = (aln.Position < aln.MatePosition);
					isRev1 = aln.IsReverseStrand();
					isRev2 = aln.IsMateReverseStrand();
					refID1 = refID2 = aln.RefID;

					//Going Downstream
					if(at5) {
						start1 = aln.Position - bpTol;
						end1 = aln.Position + bpTol;
					}
					else {
						start1 = aln.GetEndPosition() - bpTol;
						end1 = aln.GetEndPosition() + bpTol;
					}
					start2 = end1 + 1;
					end2 = start2 + bpRegion;
					down = {start1, end1, start2, end2};

					//Going Upstream
					if(at5) {
						start2 = aln.Position - bpTol;
						end2 = aln.Position + bpTol;
					}
					else {
						start2 = aln.GetEndPosition() - bpTol;
						end2 = aln.GetEndPosition() + bpTol;
					}
					end1 = start2 - 1;
					start1 = end1 - bpRegion;
					up = {start1, end1, start2, end2};

					//Upstream
					if(at5) {
						//Down bp1 ... sc/bp2
						//Up sc/bp1 ... bp2
						if((isUp && !isRev1) || (!isUp && isRev1)) {
							isSC1 = false;
							isSC2 = true;

							tmpRange = CreateRange(refID1, up[0], up[1], false, isSC1, refID2, up[2], up[3], true, isSC2, support, isShort);
							AddRange(dels_small, tmpRange);
						}
					}
					else {	//Downstream
						//Down bp1 ... bp2/sc
						//Up bp1/sc ... bp2
						if((isUp && !isRev1) || (!isUp && isRev1)) {
							isSC1 = true;
							isSC2 = false;

							tmpRange = CreateRange(refID1, down[0], down[1], false, isSC1, refID2, down[2], down[3], true, isSC2, support, isShort);
							AddRange(dels_small, tmpRange);
						}
					}
				}
			}
		}
		aln_exists = br.GetNextAlignmentCore(aln);
	}
}

bool Align(string s1, string s2) {
	int rows = s1.size() + 1;
    int cols = s2.size() + 1;

    vector < vector < int > > H(rows, vector < int > (cols, 0));

    int sim, max_score = 0;

	for(int i = 1; i < rows; i++) {
		for(int j = 1; j < cols; j++) {
			if(s1[i-1] == s2[j-1]) sim = matchScore;
			else sim = mismatchScore;
			H[i][j] = max(H[i-1][j-1] + sim,
				max(H[i-1][j] + gapPenalty,
					max(H[i][j-1] + gapPenalty, 0)));
			if(H[i][j] > max_score) max_score = H[i][j];
		}
	}

	int min_len = min(s1.size(), s2.size());
	if(max_score/2.0 >= overlap*min_len) {
		return true;
	}
	else return false;
}

bool Match(Node n1, Node n2) {
	bool here, here1, here2, here3;
	int scDist;
	string interSC1, interSC2;

	if (n1.readName != n2.readName && n1.start == n2.start) return false;

	//Same SC pos req
	//if(n1.start != n2.start) return false;

	here = here1 = here2 = here3 = false;

	//both reads at same end, sc at 3' end for 5' reads and vice versa
	if(n1.at5 == n2.at5 && n1.at5 != n1.scAt5 && n2.at5 != n2.scAt5 && abs(n1.scPos - n2.scPos) <= bpTol
		&& n1.at5 != n1.scAt5 && n2.at5 != n2.scAt5
		) {
		here = Align(n1.scSeq, n2.scSeq);
		return here;
	}
	//n1 at 5', n2 at 3'
	if(n1.at5 && !n2.at5 && !n1.scAt5 && n2.scAt5
		&& n1.scPos < n2.scPos
		) {
		here1 = Align(n1.scSeq, n2.nscSeq);
		here2 = Align(n1.nscSeq, n2.scSeq);
		//Insertion
		if(n1.scPos > (n2.scPos - n2.scLen)
			&& n2.scPos < (n1.scPos + n1.scLen)
			) {
			scDist = abs(n2.scPos - n1.scPos) + 1;
			if(scDist >= 5) {
				//cout << "st3" << endl;
				interSC1 = n1.scSeq.substr(0, scDist);
				interSC2 = n2.scSeq.substr(n2.scLen - scDist, scDist);
				//cout << "en3" << endl;
				here3 = Align(interSC1, interSC2);
				if(here3) {
					n1.ins = true;n2.ins = true;
				}
			}
		}
		return here1 & here2;
	}
	//n1 at 3', n2 at 5'
	if(!n1.at5 && n2.at5 && n1.scAt5 && !n2.scAt5
		&& n1.scPos > n2.scPos
		) {
		here1 = Align(n1.scSeq, n2.nscSeq);
		here2 = Align(n1.nscSeq, n2.scSeq);
		return here1 & here2;
	}

	return false;
}

bool Check(Node x, Cluster c) {
	for(int i = 0; i < c.n_up.size(); i++) {
		if(!Match(x, c.n_up[i])) {
			return false;
		}
	}
	for(int i = 0; i < c.n_down.size(); i++) {
		if(!Match(x, c.n_down[i])) {
			return false;
		}
	}
	return true;
}

bool CreateLinks() {
	if(co5 == 0 || co3 == 0) return false;
	if (nodes.size() == 0) return false;
	int add_co = 0;
	
	bool ret = false;
	// for (int i = 0; i < nodes.size(); i++) {
	// 	for (int j = i+1; j < nodes.size(); j++) {
	// 		if (Match(nodes[i], nodes[j])) {
	// 			add_co++;
	// 			ret = true;
	// 		}
	// 	}
	// }
	// // cerr << "Vertices: " << nodes.size() << ' ' << add_co << '\n';
	// return ret;

	
	int found;
	for(int i = 0; i < nodes.size(); i++) {
		found = -1;
		for(int j = 0; j < clusters.size(); j++) {
			if(Check(nodes[i], clusters[j])) {
				found = j;
				if(nodes[i].at5) clusters[j].n_up.push_back(nodes[i]);
				else clusters[j].n_down.push_back(nodes[i]);

				if (clusters[j].n_up.size() && clusters[j].n_down.size()) {
					ret = true;
				}

				break;
			}
		}
		if(found == -1) {
			Cluster tmp;
			if(nodes[i].at5) tmp.n_up.push_back(nodes[i]);
			else tmp.n_down.push_back(nodes[i]);
			clusters.push_back(tmp);
		}
	}
	return ret;
}

int getDelSize(BamAlignment aln) {

	int ret = 0;
	//contains del and sc
	//10M2D10M8S
	if(aln.CigarData.size() >= 4) {
		for(int i = 0; i < aln.CigarData.size(); i++) {
			if(aln.CigarData[i].Type == 'D') {
				ret += aln.CigarData[i].Length;
			}
		}
	}
	return ret;
}

int trimSeq(string &scSeq, string &scQual, bool scAt5) {
	int trimLen = 0;
	//3'
	if(!scAt5) {
		for(int i = scQual.size()-1; i >= 0 && (qual2phred[scQual[i]] - qualOffset < minScQual); i--) {
			trimLen++;
		}
		scSeq = scSeq.substr(0, scSeq.size() - trimLen);
	}	//5'
	else {
		for(int i = 0; i < scQual.size() && (qual2phred[scQual[i]] - qualOffset < minScQual); i++) {
			trimLen++;
		}
		scSeq = scSeq.substr(trimLen);
	}
	return trimLen;
}

void CreateNodes(Range r, BamReader &br) {
	int prvNodesSz = 0;
	int refID, start, end;
	int scIdx, scLen, seqLen;
	bool at5, hasSC, scAt5;
	string seq, seqQual, scSeq, nscSeq;
	BamAlignment aln;
	vector <int> clipSizes, readPositions, genomePositions;
	// refactor, remove for loop, have up_nodes, down_nodes
	for(int bp = 0; bp <= 1; bp++) {
		int totNewSC = 0;
		if(bp == 0) {
			refID = r.refID1;start = r.start1 - 1;end = r.end1 + 1;
			at5 = true;
		}
		else {
			refID = r.refID2;start = r.start2 - 1;end = r.end2 + 1;
			at5 = false;
		}

		br.SetRegion(refID, start, refID, end);
		// cout << "Nodes range : " << start << ' ' << end << endl;
		while(br.GetNextAlignmentCore(aln)) {
			clipSizes.clear();readPositions.clear();genomePositions.clear();

			hasSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);

			if(!hasSC) continue;

			scIdx = 0;
			if(clipSizes.size() == 2 && clipSizes[0] < clipSizes[1]){
				scIdx = 1;
			}

			scLen = clipSizes[scIdx];
			// remove if to continue
			if(scLen >= minSCLen && genomePositions[scIdx] >= start && genomePositions[scIdx] <= end) {

				aln.BuildCharData();

				seq = aln.QueryBases;
				seqLen = seq.size();
				seqQual = aln.Qualities;

				scAt5 = (aln.CigarData[scIdx].Type == 'S');

				int trimLen = trimSeq(seq, seqQual, scAt5);
				seqLen -= trimLen;
				scLen -= trimLen;

				if(scLen < minSCLen) continue;
				//sssnnn
				if(scAt5) {
					//cout << "st1" << endl;
					scSeq = seq.substr(0, scLen);
					nscSeq = seq.substr(scLen, seqLen - scLen);
					//cout << "en1" << endl;
				}
				else {	//nnnsss
					int scPos = readPositions[scIdx];
					scPos -= getDelSize(aln);
					//cout << "st2" << endl;
					//cout << scPos << ' ' << scLen << ' ' << seq << endl;
					scSeq = seq.substr(scPos, scLen);
					//cout << "mi2" << endl;
					nscSeq = seq.substr(0, seqLen - scLen);
					//cout << "en2" << endl;
				}

				int addNewSC = 1;
				int scCo = 0;
				for(int i = 0; i < nodes.size() && scCo <= maxSc; i++) {
					Node x = nodes[i];
					if(scAt5 == x.scAt5 && abs(genomePositions[scIdx] - x.scPos) <= bpTol) {
						scCo++;
						addNewSC = 0;
					}
				}
				totNewSC += addNewSC;
				if(totNewSC > maxCluster) {
					// cerr << co5 << ' ' << co3 << " reset\n";
					nodes.clear();
					return;
				}
				// cerr << scCo << '\n';
				if(scCo <= maxSc) {
					if(at5) co5++;
					else co3++;
					Node n;
					n.readName = aln.Name;
					n.start = aln.Position;n.end = aln.GetEndPosition();
					n.len = aln.Length;n.scLen = scLen;n.scPos = genomePositions[scIdx];
					n.seq = aln.QueryBases;n.scSeq = scSeq;n.nscSeq = nscSeq;
					n.at5 = at5;n.scAt5 = scAt5;n.ins = false;
					//cout << n.at5 << ' ' << n.scAt5 << ' ' << n.start << ' ' << n.end << ' ' << n.readName << endl;

					nodes.push_back(n);
				}
			}
		}
		//cout << "End" << endl;
	}
	// cerr << co5 << ' ' << co3 << '\n';
}

pair < int, int > Median(vector < int > v) {
	if(v.size() > 0) {
		vector < int > v_org(v);
		std::sort(v.begin(), v.end());
		int mid = v.size()/2;
		int median = v[mid];

		auto median_it = std::find(v_org.begin(), v_org.end(), median);
		int median_index = median_it - v_org.begin();
		return make_pair(median, median_index);
	}
	return {-1, -1};
}

void ComputeBreakpoints(Cluster c) {
	//ins parsing
	vector < int > bp1, bp2;
	vector < string > bp1_seq, bp2_seq;
	int start, end, size, nscSeqSize;
	string nscSeq;

	for(int i = 0; i < c.n_up.size(); i++) {
		bp1.push_back(c.n_up[i].scPos);
		nscSeq = c.n_up[i].nscSeq;
		nscSeqSize = min((int)floor(c.n_up[i].len/2), (int)nscSeq.size());
		bp1_seq.push_back(nscSeq.substr(nscSeq.size() - nscSeqSize, nscSeqSize));
		isIns &= c.n_up[i].ins;
	}
	for(int i = 0; i < c.n_down.size(); i++) {
		bp2.push_back(c.n_down[i].scPos);
		nscSeq = c.n_down[i].nscSeq;
		nscSeqSize = min((int)floor(c.n_down[i].len/2), (int)nscSeq.size());
		bp2_seq.push_back(nscSeq.substr(0, nscSeqSize));
		isIns &= c.n_down[i].ins;
	}
	pair < int, int > bp1_med = Median(bp1), bp2_med = Median(bp2);
	start = bp1_med.first;
	end = bp2_med.first;
	if(start != -1 && end != -1) {
		size = abs(end - (start+1)) + 1;
		//cout << start << '\t' << end << '\t' << size << endl;
		pred.resize(5);pred[0] = start+1;pred[1] = end;pred[2] = size;pred[4] = c.n_up.size() + c.n_down.size();
		info.resize(3);info[1] = bp1_seq[bp1_med.second];info[2] = bp2_seq[bp2_med.second];
	}
}

void ClustersParse() {
	int idx = -1, score = 0, here;
	for(int i = 0; i < clusters.size(); i++) {
		here = clusters[i].n_up.size() * clusters[i].n_down.size();
		if(here > score) {
			score = here;
			idx = i;
		}
	}
	if(idx != -1) {
		ComputeBreakpoints(clusters[idx]);
	}
}

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

	/*for(int i = 0; i < dels_small.size(); i++) {
		Range cur = dels_small[i];
		if(IsValid(cur)) {
			co5 = co3 = 0;isIns = true;
			nodes.clear();clusters.clear();pred.clear();info.clear();

			CreateNodes(cur, br);

			CreateLinks();

			ClustersParse();

			if(isIns && pred.size()) {
				cerr << "small \n";
				pred[3] = cur.support;
				info[0] = br.GetReferenceData()[cur.refID1].RefName;
				info_ins.push_back(info);
				pred_ins.push_back(pred);
				continue;
			}

			if(pred.size()) {
				pred[3] = cur.support;
				info[0] = br.GetReferenceData()[cur.refID1].RefName;
				info_dels_small.push_back(info);
				pred_dels_small.push_back(pred);
			}
		}
		else {
			//cout << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << ' ' << cur.isSC1 << ' ' << cur.isSC2 << '\t' << cur.support << endl;
		}
	}*/
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

int main(int argc, char const *argv[]) {

	args::ArgumentParser parser("This program identifies insertions and deletions in NGS data");
    args::HelpFlag help(parser, "help", "Help menu", {'h', "help"});

    args::Group groupInsertSize(parser, "Insert size parameters", args::Group::Validators::AllOrNone);
    args::ValueFlag<int> insSz(groupInsertSize, "insSz", "Insert Size", {'s', "insSz"});
    args::ValueFlag<int> stdDev(groupInsertSize, "stdDev", "Standard Deviation", {'d', "stdDev"});

    args::Group groupInputFile(parser, "Input file name", args::Group::Validators::All);
    args::ValueFlag<std::string> inpFileName(groupInputFile, "inpFile", "Input FileName", {'i', "inFile"});

    args::Group groupOutputFile(parser, "Output file name", args::Group::Validators::AllOrNone);
    args::ValueFlag<std::string> outFileName(groupOutputFile, "outFile", "Output FileName", {'o', "outFile"});
    
    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
        std::cout << parser;
        return EXIT_SUCCESS;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return EXIT_FAILURE;
    }
    catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return EXIT_FAILURE;
    }

    if (insSz) std::cout << "Insert size: " << args::get(insSz) << std::endl;
    if (stdDev) std::cout << "Standard deviation: " << args::get(stdDev) << std::endl;
    
	if (inpFileName) std::cout << "Input file: " << args::get(inpFileName) << std::endl;

	BamReader br;

	Mean = args::get(insSz);
	StdDev = args::get(stdDev);
	Input(args::get(inpFileName), br);

	auto start_t1 = std::chrono::steady_clock::now();

	Init();

	ProcessBam(br);

	auto end_t1 = std::chrono::steady_clock::now();
	cerr << "Step 1 : Done (" << setprecision(1) << (std::chrono::duration_cast<std::chrono::minutes>(end_t1 - start_t1).count()) << " mins)" << endl;

	auto start_t2 = std::chrono::steady_clock::now();

	DelsParse(br);

	auto end_t2 = std::chrono::steady_clock::now();
	cerr << "Step 2 : Done (" << setprecision(1) << (std::chrono::duration_cast<std::chrono::minutes>(end_t2 - start_t2).count()) << " mins)" << endl;

	br.Close();

	Output();

	return EXIT_SUCCESS;
}
