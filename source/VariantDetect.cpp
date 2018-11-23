#include "VariantDetect.hpp"

inline bool inBetween(const int &a, const int &b, const int &x)
{
	return a <= x && x <= b;
}

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

bool MatchInCluster(const Node &n1, const Node &n2)
{
	// both should be at either start or end
	if (n1.at5 != n2.at5)
		return false;
	if (n1.at5 == n1.scAt5 || n2.at5 == n2.scAt5)
		return false;
	if (abs(n1.scPos - n2.scPos) <= bpTol)
	{
		return Align(n1.scSeq, n2.scSeq);
	}
	return false;
}

bool FindCluster(const Node &n, const std::vector<Node> &cluster)
{
	for (int i = 0; i < cluster.size(); i++)
	{
		if (!MatchInCluster(n, cluster[i]))
			return false;
	}
	return true;
}

void ClusterSC(std::vector<Node> &nodes, std::vector<std::vector<Node>> &clusters)
{
	for (int i = 0; i < nodes.size(); ++i)
	{
		bool found = false;
		for (int j = 0; j < clusters.size(); j++)
		{
			if (FindCluster(nodes[i], clusters[j]))
			{
				found = true;
				clusters[j].push_back(nodes[i]);
				if (clusters[j].size() > 100)
				{
					clusters[j].clear();
				}
				break;
			}
		}
		if (!found)
		{
			clusters.push_back({nodes[i]});
			if (clusters.size() > 175)
			{
				clusters.clear();
				return;
			}
		}
	}
}

bool MatchOutNodes(const Node &nUp, const Node &nDown)
{
	// both are at same end
	if (nUp.at5 == nDown.at5)
		return false;
	if (nUp.scAt5 == nDown.scAt5)
		return false;
	if (nUp.at5 && !nDown.at5 && nUp.scPos < nDown.scPos)
	{
		bool match1 = Align(nUp.scSeq, nDown.nscSeq);
		bool match2 = Align(nUp.nscSeq, nDown.scSeq);
		return match1 | match2;
	}
	return false;
}

bool ClustersMatch(const std::vector<Node> &up, const std::vector<Node> &down)
{
	for (int i = 0; i < up.size(); i++)
	{
		Node nUp = up[i];
		for (int j = 0; j < down.size(); j++)
		{
			Node nDown = down[j];
			if (!MatchOutNodes(nUp, nDown))
				return false;
		}
	}
	return true;
}

void ExtractSC(BamTools::BamReader &br, const int &refID, const int &start, const int &end, std::vector<std::vector<Node>> &clusters, const bool &at5)
{
	std::vector<Node> nodes;

	BamTools::BamAlignment aln;
	std::vector<int> clipSizes, readPositions, genomePositions;

	br.SetRegion(refID, start, refID, end);
	while (br.GetNextAlignmentCore(aln))
	{
		clipSizes.clear();
		readPositions.clear();
		genomePositions.clear();

		bool hasSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
		if (!hasSC)
			continue;

		int scIdx = 0;
		if (clipSizes.size() == 2 && clipSizes[0] < clipSizes[1])
		{
			scIdx = 1;
		}

		int scLen = clipSizes[scIdx];
		if (scLen < minSCLen)
			continue;

		if (!inBetween(start, end, genomePositions[scIdx]))
			continue;

		aln.BuildCharData();
		std::string seq = aln.QueryBases;
		int seqLen = seq.size();
		std::string seqQual = aln.Qualities;

		bool scAtLeft = (aln.CigarData[scIdx].Type == 'S');

		int trimLen = trimSeq(seq, seqQual, scAtLeft);
		seqLen -= trimLen;
		scLen -= trimLen;

		if (scLen < minSCLen)
			continue;

		std::string scSeq, nscSeq;
		if (scAtLeft)
		{
			scSeq = seq.substr(0, scLen);
			nscSeq = seq.substr(scLen, seqLen - scLen);
		}
		else
		{
			int scPos = readPositions[scIdx];
			scPos -= getDelSize(aln);
			scSeq = seq.substr(scPos, scLen);
			nscSeq = seq.substr(0, seqLen - scLen);
		}

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
		n.scAt5 = scAtLeft;
		n.at5 = at5;
		n.ins = false;

		nodes.push_back(n);
		if (nodes.size() > 500)
		{
			nodes.clear();
			return;
		}
	}

	ClusterSC(nodes, clusters);
}

void MedianBP(const std::vector<Node> &nodes, int &bpMedian)
{
	std::vector<int> bp;
	for (int i = 0; i < nodes.size(); i++)
	{
		bp.push_back(nodes[i].scPos);
	}

	int mid = bp.size() / 2;
	std::nth_element(bp.begin(), bp.begin() + mid, bp.end());
	bpMedian = bp[mid];
}

void ParseSC(BamTools::BamReader &br, const DiscNode &r, std::vector<std::string> &out)
{
	std::vector<std::vector<Node>> clusterUp, clusterDown;
	ExtractSC(br, r.refID, r.upStart - 1, r.upEnd + 1, clusterUp, true);
	if (clusterUp.size() == 0)
		return;
	
	ExtractSC(br, r.refID, r.downStart - 1, r.downEnd + 1, clusterDown, false);
	if (clusterDown.size() == 0)
		return;

	int maxScore = 0, maxIdxUp = -1, maxIdxDown = -1;

	for (int i = 0; i < clusterUp.size(); i++)
	{
		for (int j = 0; j < clusterDown.size(); j++)
		{
			if (ClustersMatch(clusterUp[i], clusterDown[j]))
			{
				int hereScore = clusterUp[i].size() * clusterDown[j].size();
				if (hereScore > maxScore)
				{
					maxScore = hereScore;
					maxIdxUp = i;
					maxIdxDown = j;
				}
			}
		}
	}

	if (maxScore != 0)
	{
		int bpUp, bpDown;
		MedianBP(clusterUp[maxIdxUp], bpUp);
		MedianBP(clusterDown[maxIdxDown], bpDown);

		out.resize(6);
		out[0] = br.GetReferenceData()[r.refID].RefName;
		out[1] = std::to_string(bpUp + 1);
		out[2] = std::to_string(bpDown);
		out[3] = std::to_string(bpDown - (bpUp + 1) + 1);
		out[4] = std::to_string(r.support);
		out[5] = std::to_string(clusterUp[maxIdxUp].size() + clusterDown[maxIdxDown].size());
	}

	/*
	std::cout << "Upstream\n";
	for (int i = 0; i < clusterUp.size(); i++)
	{
		std::cout << clusterUp[i].size() << '\n';
		Node hc;
		for (int j = 0; j < clusterUp[i].size(); j++)
		{
			hc = clusterUp[i][j];
			std::cout << hc.scPos << ' ' << hc.at5 << ' ' << hc.scAt5 << ' ' << hc.scSeq << ' ' << hc.nscSeq << '\n';
		}
	}
	std::cout << "Downstream\n";
	for (int i = 0; i < clusterDown.size(); i++)
	{
		std::cout << clusterDown[i].size() << '\n';
		Node hc;
		for (int j = 0; j < clusterDown[i].size(); j++)
		{
			hc = clusterDown[i][j];
			std::cout << hc.scPos << ' ' << hc.at5 << ' ' << hc.scAt5 << ' ' << hc.nscSeq << ' ' << hc.scSeq << '\n';
		}
	}
	*/
}

std::vector<std::vector<std::string>> DelsParseParallel(std::string fPath, int st, int en, const std::vector<DiscCluster> dels)
{
	BamTools::BamReader br;
	br.Open(fPath);
	br.OpenIndex(fPath + ".bai");

	std::vector<std::vector<std::string>> ret;
	for (int i = st; i < en; i++)
	{
		DiscNode cur = dels[i].info;

		std::vector<std::string> out;
		ParseSC(br, cur, out);

		if (out.size())
		{
			ret.push_back(out);
		}
	}
	br.Close();
	return ret;
}

void DelsParse(std::string fPath, const std::vector<DiscCluster> &dels, std::vector<std::vector<std::string>> &output)
{
	unsigned int th = 4;
	int till = dels.size();

	transwarp::parallel ex{th};
	std::vector<std::shared_ptr<transwarp::task<std::vector<std::vector<std::string>>>>> tasks;

	for (int i = 0; i < th; i++)
	{
		int st = (till / th) * i;
		int en = (i == th - 1) ? dels.size() : (till / th) * (i + 1);
		auto fpathTask = transwarp::make_value_task(fPath);
		auto stTask = transwarp::make_value_task(st);
		auto enTask = transwarp::make_value_task(en);
		auto delsTask = transwarp::make_value_task(dels);
		auto fnTask = transwarp::make_task(transwarp::consume, DelsParseParallel, fpathTask, stTask, enTask, delsTask);
		tasks.emplace_back(fnTask);
	}

	for (auto task : tasks)
	{
		task->schedule(ex);
	}
	for (auto task : tasks)
	{
		std::vector<std::vector<std::string>> outputPL = task->get();
		for (std::vector<std::string> vs : outputPL)
		{
			output.push_back(vs);
		}
	}
}
