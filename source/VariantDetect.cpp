#include "VariantDetect.hpp"

bool IsValid(const Range &r)
{
	//Large SV
	if (!r.isShort)
	{
		//bp region shouldn't be less than 10bp
		if (r.isSC1 || r.isSC2 || abs(r.end1 - r.start1 + 1) < 10 || abs(r.end2 - r.start2 + 1) < 10)
			return false;
		else
			return true;
	}
	else
	{ //Small SV
		//Both SC need to be present
		if (r.isSC1 && r.isSC2)
			return true;
		else
			return false;
	}
}

int getDelSize(BamTools::BamAlignment &aln)
{

	int ret = 0;
	//contains del and sc
	//10M2D10M8S
	if (aln.CigarData.size() >= 4)
	{
		for (int i = 0; i < aln.CigarData.size(); i++)
		{
			if (aln.CigarData[i].Type == 'D')
			{
				ret += aln.CigarData[i].Length;
			}
		}
	}
	return ret;
}

int trimSeq(std::string &scSeq, std::string &scQual, bool scAt5)
{
	int trimLen = 0;
	//3'
	if (!scAt5)
	{
		for (int i = scQual.size() - 1; i >= 0 && (qual2phred.at(scQual[i]) - qualOffset < minScQual); i--)
		{
			trimLen++;
		}
		scSeq = scSeq.substr(0, scSeq.size() - trimLen);
	} //5'
	else
	{
		for (int i = 0; i < scQual.size() && (qual2phred.at(scQual[i]) - qualOffset < minScQual); i++)
		{
			trimLen++;
		}
		scSeq = scSeq.substr(trimLen);
	}
	return trimLen;
}

bool Align(const std::string &s1, const std::string &s2)
{
	int rows = s1.size() + 1;
	int cols = s2.size() + 1;

	std::vector<std::vector<int>> H(rows, std::vector<int>(cols, 0));

	int sim, max_score = 0;

	for (int i = 1; i < rows; i++)
	{
		for (int j = 1; j < cols; j++)
		{
			if (s1[i - 1] == s2[j - 1])
				sim = matchScore;
			else
				sim = mismatchScore;
			H[i][j] = std::max(H[i - 1][j - 1] + sim,
							   std::max(H[i - 1][j] + gapPenalty,
										std::max(H[i][j - 1] + gapPenalty, 0)));
			if (H[i][j] > max_score)
				max_score = H[i][j];
		}
	}

	int min_len = std::min(s1.size(), s2.size());
	if (max_score / 2.0 >= overlap * min_len)
	{
		return true;
	}
	else
		return false;
}

bool Match(Node &n1, Node &n2)
{
	bool here, here1, here2, here3;
	int scDist;
	std::string interSC1, interSC2;

	if (n1.readName != n2.readName && n1.start == n2.start)
		return false;

	//Same SC pos req
	//if(n1.start != n2.start) return false;

	here = here1 = here2 = here3 = false;

	//both reads at same end, sc at 3' end for 5' reads and vice versa
	if (n1.at5 == n2.at5 && n1.at5 != n1.scAt5 && n2.at5 != n2.scAt5 && abs(n1.scPos - n2.scPos) <= bpTol && n1.at5 != n1.scAt5 && n2.at5 != n2.scAt5)
	{
		here = Align(n1.scSeq, n2.scSeq);
		return here;
	}
	//n1 at 5', n2 at 3'
	if (n1.at5 && !n2.at5 && !n1.scAt5 && n2.scAt5 && n1.scPos < n2.scPos)
	{
		here1 = Align(n1.scSeq, n2.nscSeq);
		here2 = Align(n1.nscSeq, n2.scSeq);
		//Insertion
		if (n1.scPos > (n2.scPos - n2.scLen) && n2.scPos < (n1.scPos + n1.scLen))
		{
			scDist = abs(n2.scPos - n1.scPos) + 1;
			if (scDist >= 5)
			{
				//cout << "st3" << endl;
				interSC1 = n1.scSeq.substr(0, scDist);
				interSC2 = n2.scSeq.substr(n2.scLen - scDist, scDist);
				//cout << "en3" << endl;
				here3 = Align(interSC1, interSC2);
				if (here3)
				{
					n1.ins = true;
					n2.ins = true;
				}
			}
		}
		return here1 & here2;
	}
	//n1 at 3', n2 at 5'
	if (!n1.at5 && n2.at5 && n1.scAt5 && !n2.scAt5 && n1.scPos > n2.scPos)
	{
		here1 = Align(n1.scSeq, n2.nscSeq);
		here2 = Align(n1.nscSeq, n2.scSeq);
		return here1 & here2;
	}

	return false;
}

bool Check(Node &x, Cluster &c)
{
	for (int i = 0; i < c.n_up.size(); i++)
	{
		if (!Match(x, c.n_up[i]))
		{
			return false;
		}
	}
	for (int i = 0; i < c.n_down.size(); i++)
	{
		if (!Match(x, c.n_down[i]))
		{
			return false;
		}
	}
	return true;
}

bool CreateLinks(std::vector<Node> &nodes, std::vector<Cluster> &clusters, const int &co5, const int &co3)
{
	if (co5 == 0 || co3 == 0)
		return false;
	if (nodes.size() == 0)
		return false;
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
	for (int i = 0; i < nodes.size(); i++)
	{
		found = -1;
		for (int j = 0; j < clusters.size(); j++)
		{
			if (Check(nodes[i], clusters[j]))
			{
				found = j;
				if (nodes[i].at5)
					clusters[j].n_up.push_back(nodes[i]);
				else
					clusters[j].n_down.push_back(nodes[i]);

				if (clusters[j].n_up.size() && clusters[j].n_down.size())
				{
					ret = true;
				}

				break;
			}
		}
		if (found == -1)
		{
			Cluster tmp;
			if (nodes[i].at5)
				tmp.n_up.push_back(nodes[i]);
			else
				tmp.n_down.push_back(nodes[i]);
			clusters.push_back(tmp);
		}
	}
	return ret;
}

void CreateNodes(const Range &r, BamTools::BamReader &br, std::vector<Node> &nodes, int &co5, int &co3)
{
	int prvNodesSz = 0;
	int refID, start, end;
	int scIdx, scLen, seqLen;
	bool at5, hasSC, scAt5;
	std::string seq, seqQual, scSeq, nscSeq;
	BamTools::BamAlignment aln;
	std::vector<int> clipSizes, readPositions, genomePositions;
	// refactor, remove for loop, have up_nodes, down_nodes
	for (int bp = 0; bp <= 1; bp++)
	{
		int totNewSC = 0;
		if (bp == 0)
		{
			refID = r.refID1;
			start = r.start1 - 1;
			end = r.end1 + 1;
			at5 = true;
		}
		else
		{
			refID = r.refID2;
			start = r.start2 - 1;
			end = r.end2 + 1;
			at5 = false;
		}

		br.SetRegion(refID, start, refID, end);
		// cout << "Nodes range : " << start << ' ' << end << endl;
		while (br.GetNextAlignmentCore(aln))
		{
			clipSizes.clear();
			readPositions.clear();
			genomePositions.clear();

			hasSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);

			if (!hasSC)
				continue;

			scIdx = 0;
			if (clipSizes.size() == 2 && clipSizes[0] < clipSizes[1])
			{
				scIdx = 1;
			}

			scLen = clipSizes[scIdx];
			// remove if to continue
			if (scLen >= minSCLen && genomePositions[scIdx] >= start && genomePositions[scIdx] <= end)
			{

				aln.BuildCharData();

				seq = aln.QueryBases;
				seqLen = seq.size();
				seqQual = aln.Qualities;

				scAt5 = (aln.CigarData[scIdx].Type == 'S');

				int trimLen = trimSeq(seq, seqQual, scAt5);
				seqLen -= trimLen;
				scLen -= trimLen;

				if (scLen < minSCLen)
					continue;
				//sssnnn
				if (scAt5)
				{
					//cout << "st1" << endl;
					scSeq = seq.substr(0, scLen);
					nscSeq = seq.substr(scLen, seqLen - scLen);
					//cout << "en1" << endl;
				}
				else
				{ //nnnsss
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
				for (int i = 0; i < nodes.size() && scCo <= maxSc; i++)
				{
					Node x = nodes[i];
					if (scAt5 == x.scAt5 && abs(genomePositions[scIdx] - x.scPos) <= bpTol)
					{
						scCo++;
						addNewSC = 0;
					}
				}
				totNewSC += addNewSC;
				if (totNewSC > maxCluster)
				{
					// cerr << co5 << ' ' << co3 << " reset\n";
					nodes.clear();
					return;
				}
				// cerr << scCo << '\n';
				if (scCo <= maxSc)
				{
					if (at5)
						co5++;
					else
						co3++;
					Node n;
					n.readName = aln.Name;
					n.start = aln.Position;
					n.end = aln.GetEndPosition();
					n.len = aln.Length;
					n.scLen = scLen;
					n.scPos = genomePositions[scIdx];
					n.seq = aln.QueryBases;
					n.scSeq = scSeq;
					n.nscSeq = nscSeq;
					n.at5 = at5;
					n.scAt5 = scAt5;
					n.ins = false;
					//cout << n.at5 << ' ' << n.scAt5 << ' ' << n.start << ' ' << n.end << ' ' << n.readName << endl;

					nodes.push_back(n);
				}
			}
		}
		//cout << "End" << endl;
	}
	// cerr << co5 << ' ' << co3 << '\n';
}

std::pair<int, int> Median(std::vector<int> &v)
{
	if (v.size() > 0)
	{
		std::vector<int> v_org(v);
		std::sort(v.begin(), v.end());
		int mid = v.size() / 2;
		int median = v[mid];

		auto median_it = std::find(v_org.begin(), v_org.end(), median);
		int median_index = median_it - v_org.begin();
		return std::make_pair(median, median_index);
	}
	return {-1, -1};
}

void ComputeBreakpoints(const Cluster &c, bool &isIns, std::vector<int> &pred, std::vector<std::string> &info)
{
	//ins parsing
	std::vector<int> bp1, bp2;
	std::vector<std::string> bp1_seq, bp2_seq;
	int start, end, size, nscSeqSize;
	std::string nscSeq;

	for (int i = 0; i < c.n_up.size(); i++)
	{
		bp1.push_back(c.n_up[i].scPos);
		nscSeq = c.n_up[i].nscSeq;
		nscSeqSize = std::min((int)floor(c.n_up[i].len / 2), (int)nscSeq.size());
		bp1_seq.push_back(nscSeq.substr(nscSeq.size() - nscSeqSize, nscSeqSize));
		isIns &= c.n_up[i].ins;
	}
	for (int i = 0; i < c.n_down.size(); i++)
	{
		bp2.push_back(c.n_down[i].scPos);
		nscSeq = c.n_down[i].nscSeq;
		nscSeqSize = std::min((int)floor(c.n_down[i].len / 2), (int)nscSeq.size());
		bp2_seq.push_back(nscSeq.substr(0, nscSeqSize));
		isIns &= c.n_down[i].ins;
	}
	std::pair<int, int> bp1_med = Median(bp1), bp2_med = Median(bp2);
	start = bp1_med.first;
	end = bp2_med.first;
	if (start != -1 && end != -1)
	{
		size = abs(end - (start + 1)) + 1;
		//cout << start << '\t' << end << '\t' << size << endl;
		pred.resize(5);
		pred[0] = start + 1;
		pred[1] = end;
		pred[2] = size;
		pred[4] = c.n_up.size() + c.n_down.size();
		info.resize(3);
		info[1] = bp1_seq[bp1_med.second];
		info[2] = bp2_seq[bp2_med.second];
	}
}

void ClustersParse(const std::vector<Cluster> &clusters, bool &isIns, std::vector<int> &pred, std::vector<std::string> &info)
{
	int idx = -1, score = 0, here;
	for (int i = 0; i < clusters.size(); i++)
	{
		here = clusters[i].n_up.size() * clusters[i].n_down.size();
		if (here > score)
		{
			score = here;
			idx = i;
		}
	}
	if (idx != -1)
	{
		ComputeBreakpoints(clusters[idx], isIns, pred, info);
	}
}

void DelsParse(BamTools::BamReader &br, const std::vector<Range> &dels_large,
			   std::vector<std::vector<std::string>> &info_dels_large,
			   std::vector<std::vector<std::string>> &info_dels_small,
			   std::vector<std::vector<std::string>> &info_ins,
			   std::vector<std::vector<int>> &pred_dels_large,
			   std::vector<std::vector<int>> &pred_dels_small,
			   std::vector<std::vector<int>> &pred_ins)
{

	int co5 = 0, co3 = 0;
	bool isIns = false;

	std::vector<Node> nodes;
	std::vector<Cluster> clusters;

	std::vector<std::string> info;
	std::vector<int> pred;

	std::set<int> sl_dels;
	auto prv = std::chrono::steady_clock::now();
	std::cerr << "Size SV Ranges : " << dels_large.size() << '\n';
	int co = 0, edges_co = 0;
	for (int i = 0; i < dels_large.size(); i++)
	{
		// if ( i > 10 ) break;
		Range cur = dels_large[i];
		if (IsValid(cur))
		{
			co++;
			// cerr << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << ' ' << cur.isSC1 << ' ' << cur.isSC2 << '\t' << cur.support << endl;
			co5 = co3 = 0;
			isIns = true;
			nodes.clear();
			clusters.clear();
			pred.clear();
			info.clear();

			// cerr << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << ' ';

			CreateNodes(cur, br, nodes, co5, co3);

			if (CreateLinks(nodes, clusters, co5, co3))
			{
				// cerr << cur.start1 << ' ' << cur.end1 << ' ' << cur.start2 << ' ' << cur.end2 << '\n';
				edges_co++;
			}

			// for (int cl = 0; cl < clusters.size(); cl++) cerr << "{" << clusters[cl].n_up.size() << ',' << clusters[cl].n_down.size() << "} ";
			// cerr << '\n';

			ClustersParse(clusters, isIns, pred, info);

			if (isIns && pred.size())
			{
				std::cerr << "large \n";
				pred[3] = cur.support;
				info[0] = br.GetReferenceData()[cur.refID1].RefName;
				info_ins.push_back(info);
				pred_ins.push_back(pred);
				continue;
			}

			if (pred.size())
			{
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
		int percentDone = (((i + 1) * 100 / dels_large.size()) / 5) * 5;
		if (percentDone != 0 && percentDone % 5 == 0 && sl_dels.find(percentDone) == sl_dels.end())
		{
			sl_dels.insert(percentDone);
			auto cur = std::chrono::steady_clock::now();
			// cerr << percentDone << "%\r" << std::flush;// << ' ' << setprecision(1) << (std::chrono::duration_cast<std::chrono::minutes>(cur - prv).count()) << endl;
			prv = cur;
		}
	}
	std::cerr << "Valid: " << co << '\n';
	std::cerr << "Edges built : " << edges_co << '\n';
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
