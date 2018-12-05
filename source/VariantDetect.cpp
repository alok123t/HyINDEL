#include "VariantDetect.hpp"

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

		bool down = false;
		if (aln.Position > aln.MatePosition)
		{
			down = true;
		}

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

		int scPos = genomePositions[scIdx];

		std::string scSeq, nscSeq;
		if (scAtLeft)
		{
			scSeq = seq.substr(0, scLen);
			nscSeq = seq.substr(scLen, seqLen - scLen);
		}
		else
		{
			scPos = readPositions[scIdx];
			scPos -= getDelSize(aln);
			scSeq = seq.substr(scPos, scLen);
			nscSeq = seq.substr(0, seqLen - scLen);
		}

		std::string refName = br.GetReferenceData()[aln.RefID].RefName;
		SoftNode sn(scPos, aln.RefID, refName, scSeq, nscSeq, down, at5);

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

std::vector<std::vector<std::string>> DelsLargeParallel(std::string fPath, int st, int en, const std::vector<DiscCluster> dels)
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

void DelsLarge(std::string fPath, const std::vector<DiscCluster> &dels, std::vector<std::vector<std::string>> &output)
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
		auto fnTask = transwarp::make_task(transwarp::consume, DelsLargeParallel, fpathTask, stTask, enTask, delsTask);
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

void SoftMedianBP(const std::vector<SoftNode> &nodes, int &bpMedian)
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

bool MatchOutClusters(const SoftCluster &c1, const SoftCluster &c2)
{
	int diffLen = c2.info.scPos - c1.info.scPos;
	if (diffLen < MIN_DEL_LEN || diffLen > MAX_DEL_LEN)
		return false;

	std::string nLeft, nRight, curLeft, curRight;
	for (int i = 0; i < c1.nodes.size(); ++i)
	{
		if (c1.nodes.at(i).scAtRight)
		{
			nLeft = c1.nodes.at(i).nscSeq;
			nRight = c1.nodes.at(i).scSeq;
		}
		else
		{
			nLeft = c1.nodes.at(i).scSeq;
			nRight = c1.nodes.at(i).nscSeq;
		}
		for (int j = 0; j < c2.nodes.size(); ++j)
		{
			if (c2.nodes.at(j).scAtRight)
			{
				curLeft = c2.nodes.at(j).nscSeq;
				curRight = c2.nodes.at(j).scSeq;
			}
			else
			{
				curLeft = c2.nodes.at(j).scSeq;
				curRight = c2.nodes.at(j).nscSeq;
			}
			if (!Align(nLeft, curLeft))
				return false;
			if (!Align(nRight, curRight))
				return false;
		}
	}
	return true;
}

void FindInCluster(const SoftCluster &sc1, const std::vector<SoftCluster> &c, std::vector<std::string> &out)
{
	int maxScore = 0, maxIdx = -1;
	for (int i = 0; i < c.size(); ++i)
	{
		SoftCluster sc2 = c[i];
		if (MatchOutClusters(sc1, sc2))
		{
			int hereScore = sc1.nodes.size() * sc2.nodes.size();
			if (hereScore > maxScore)
			{
				maxScore = hereScore;
				maxIdx = i;
			}
		}
	}
	if (maxScore != 0)
	{
		int bpUp, bpDown;
		SoftMedianBP(sc1.nodes, bpUp);
		SoftMedianBP(c[maxIdx].nodes, bpDown);

		out.resize(6);
		out[0] = sc1.info.refName;
		out[1] = std::to_string(bpUp + 1);
		out[2] = std::to_string(bpDown);
		out[3] = std::to_string(bpDown - (bpUp + 1) + 1);
		out[4] = std::string("0");
		int scSup = sc1.nodes.size() + c[maxIdx].nodes.size();
		if (scSup < MIN_SC_SUP)
		{
			out.clear();
			return;
		}
		out[5] = std::to_string(scSup);
	}
}

void DelsSmall(const std::vector<std::vector<SoftCluster>> &cUp, const std::vector<std::vector<SoftCluster>> &cDown, std::vector<std::vector<std::string>> &output)
{
	for (int i = 0; i < cUp.size(); ++i)
	{
		for (int j = 0; j < cUp.at(i).size(); ++j)
		{
			std::vector<std::string> out;
			FindInCluster(cUp.at(i).at(j), cDown.at(i), out);
			if (out.size())
			{
				output.push_back(out);
				continue;
			}

			// search in next bin
			if (j + 1 >= cDown.size())
				continue;

			FindInCluster(cUp.at(i).at(j), cDown.at(i + 1), out);
			if (out.size())
			{
				output.push_back(out);
				continue;
			}
		}
	}
}
