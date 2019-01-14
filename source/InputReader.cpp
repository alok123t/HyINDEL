#include "InputReader.hpp"

bool overlapDiscRange(const DiscNode &r, const std::vector<DiscNode> &nodes)
{
	const double RO = 0.75;

	for (int i = 0; i < nodes.size(); i++)
	{
		double len1 = (r.downEnd - r.upStart + 1) * 1.0;
		double len2 = (nodes[i].downEnd - nodes[i].upStart + 1) * 1.0;

		int st = std::max(r.upStart, nodes[i].upStart);
		int en = std::min(r.downEnd, nodes[i].downEnd);
		double len = (en - st + 1) * 1.0;
		double overlap = (len > 0.0) ? len : 0.0;

		double ro1 = overlap / len1;
		double ro2 = overlap / len2;
		// std::cerr << ro1 << ' ' << ro2 << '\n';
		bool ok = ro1 >= RO && ro2 >= RO;
		if (!ok)
			return false;
	}
	return true;
}

void clusterDiscRange(const DiscNode &r, std::vector<DiscCluster> &c)
{
	for (int i = c.size() - 1; i >= 0; i--)
	{
		// diff chr
		if (c[i].info.refID != r.refID)
			break;
		// not overlapping
		if ((c[i].info.downEnd < r.upStart) || (r.downEnd < c[i].info.upStart))
			continue;
		if (overlapDiscRange(r, c[i].nodes))
		{
			// update info
			c[i].info.upStart = std::min(c[i].info.upStart, r.upStart);
			c[i].info.upEnd = std::max(c[i].info.upEnd, r.upEnd);
			c[i].info.downStart = std::min(c[i].info.downStart, r.downStart);
			c[i].info.downEnd = std::max(c[i].info.downEnd, r.downEnd);
			c[i].info.support += r.support;
			// push new node
			c[i].nodes.emplace_back(r);
			return;
		}
	}
	// add new cluster
	DiscCluster dc(r);
	c.emplace_back(dc);
}

void DelsLargeSplit(const BamTools::RefVector &ref, const std::vector<DiscCluster> &discClusters, const std::vector<SplitCluster> &splitClusters, std::vector<int> &usedSR, std::vector<OutNode> &output)
{
	for (int i = 0; i < splitClusters.size(); ++i)
	{
		if (usedSR.at(i))
			continue;
		SplitCluster sc = splitClusters.at(i);
		if (sc.nodes.size() < 2)
			continue;
		OutNode on(ref.at(sc.info.refID).RefName, sc.info.st, sc.info.en, sc.info.en - sc.info.st + 1, 0, sc.nodes.size(), 0);
		output.emplace_back(on);
	}
}

void deletionsSR(const BamTools::RefVector &ref, const std::vector<DiscCluster> &discClusters, const std::vector<SplitCluster> &splitClusters, const std::string inpFilePath, const std::string folderPath)
{
	std::vector<int> usedSR(splitClusters.size(), 0);
	std::vector<OutNode> outputDelsSplitLarge;
	DelsLargeSplit(ref, discClusters, splitClusters, usedSR, outputDelsSplitLarge);

	// print deletion to file
	std::string outputFile = folderPath;
	parseOutputN(outputFile, outputDelsSplitLarge, 2);
}

inline void splitMultipleTag(const std::string &tag, std::vector<std::string> &tags)
{
	std::size_t idx = tag.find_first_not_of(SEMICOLON, 0);
	std::size_t semicolonIdx = tag.find_first_of(SEMICOLON, idx);

	while (idx != std::string::npos || semicolonIdx != std::string::npos)
	{
		tags.emplace_back(tag.substr(idx, semicolonIdx - idx));

		idx = tag.find_first_not_of(SEMICOLON, semicolonIdx);
		semicolonIdx = tag.find_first_of(SEMICOLON, idx);
	}
}

inline void splitTag(const std::vector<std::string> &tags, std::vector<std::vector<std::string>> &infoTags)
{
	for (std::string tag : tags)
	{
		std::vector<std::string> infoTag;
		std::size_t idx = tag.find_first_not_of(COMMA, 0);
		std::size_t commaIdx = tag.find_first_of(COMMA, idx);

		while (idx != std::string::npos || commaIdx != std::string::npos)
		{
			infoTag.emplace_back(tag.substr(idx, commaIdx - idx));

			idx = tag.find_first_not_of(COMMA, commaIdx);
			commaIdx = tag.find_first_of(COMMA, idx);
		}
		infoTags.emplace_back(infoTag);
	}
}

int parseCigar(std::string cig)
{
	int ret = 0;
	std::string s = "";
	for (char c : cig)
	{
		if (isalpha(c))
		{
			int val = std::stoi(s);
			s = "";
			if (c == 'M')
				ret += val;
			else if (c == 'S')
				break;
			else if (c == 'I')
				continue;
			else if (c == 'D')
				ret -= val;
		}
		else
			s += c;
	}
	return ret;
}

void clusterSR(const std::vector<SplitNode> &nodes, std::vector<SplitCluster> &clusters)
{
	for (int i = 0; i < nodes.size(); ++i)
	{
		bool found = false;
		for (int j = clusters.size() - 1; j >= 0; j--)
		{
			int d1 = std::abs(clusters.at(j).info.st - nodes.at(i).st);
			int d2 = std::abs(clusters.at(j).info.en - nodes.at(i).en);
			if (d1 <= 5 && d2 <= 5)
			{
				found = true;
				clusters.at(j).nodes.emplace_back(nodes.at(i));
				break;
			}
			if (d1 > MAX_SPLIT_LEN || d2 > MAX_SPLIT_LEN)
				break;
		}
		if (!found)
		{
			SplitCluster sc(nodes.at(i));
			clusters.emplace_back(sc);
		}
	}

	std::sort(clusters.begin(), clusters.end(), SplitCmp);
}

void clusterDisc(const std::vector<DiscNode> &nodes, std::vector<DiscCluster> &clusters)
{
	std::vector<DiscCluster> discClusters;
	for (DiscNode dn : nodes)
	{
		if (dn.downEnd - dn.upStart > MAX_LARGE_DEL_LEN)
			continue;
		clusterDiscRange(dn, discClusters);
	}
	for (DiscCluster dc : discClusters)
	{
		if (dc.nodes.size() > MIN_DISC_CLUSTER_SUPPORT)
			clusters.emplace_back(dc);
	}
}

int alignSeq(const std::string &a, const std::string &b)
{
	int rows = a.size() + 1;
	int cols = b.size() + 1;

	std::vector<std::vector<int>> m(rows, std::vector<int>(cols));
	int maxScore = 0;
	for (int i = 0; i < rows; ++i)
		m.at(i).at(0) = 0;
	for (int j = 0; j < cols; ++j)
		m.at(0).at(j) = j * a_gap;
	for (int i = 1; i < rows; ++i)
	{
		for (int j = 1; j < cols; ++j)
		{
			int sim = a.at(i - 1) == b.at(j - 1) ? a_matchScore : a_misMatchScore;
			m.at(i).at(j) = std::max(sim + m.at(i - 1).at(j - 1),
									 std::max(a_gap + m.at(i - 1).at(j), a_gap + m.at(i).at(j - 1)));
			maxScore = std::max(maxScore, m.at(i).at(j));
		}
	}
	return maxScore;
}

bool overlapSC(const SoftNode &n, const SoftCluster &c)
{
	int sum = 0;
	for (SoftNode cn : c.nodes)
	{
		std::string a = n.seq, b = cn.seq;
		if (n.scAtRight)
		{
			if (n.start > cn.start)
			{
				a = cn.seq;
				b = n.seq;
			}
		}
		else
		{
			if (n.end > cn.end)
			{
				a = cn.seq;
				b = n.seq;
			}
		}
		int here = alignSeq(a, b);
		if (here < interMinAlignScore)
			return false;
		sum += here;
	}
	double avgScore = 1.0 * sum / c.nodes.size();
	return avgScore >= interAlignScore;
}

void clusterSC(const std::vector<SoftNode> &nodes, std::vector<SoftCluster> &clusters)
{
	for (SoftNode n : nodes)
	{
		bool found = false;
		for (int i = clusters.size() - 1; i >= 0; i--)
		{
			// diff chr
			if (clusters.at(i).info.refID != n.refID)
				break;
			// too far apart, will not find a cluster
			if (n.scPos > clusters.at(i).info.scPos + 100)
				break;
			if (abs(n.scPos - clusters.at(i).info.scPos) > MIN_SC_DISTANCE)
				continue;
			if (overlapSC(n, clusters.at(i)))
			{
				found = true;
				clusters.at(i).nodes.emplace_back(n);
				break;
			}
		}
		// add new cluster
		if (!found)
		{
			SoftCluster sc(n);
			clusters.emplace_back(n);
		}
	}
}

void addDisc(BamTools::BamAlignment &aln, const int &bpRegion, std::vector<DiscNode> &nodes)
{
	if (aln.IsProperPair())
		return;

	int upSt = aln.Position;
	int upEn = aln.GetEndPosition();
	bool upIsRev = aln.IsReverseStrand();
	int downSt = aln.MatePosition;
	int downEn = downSt + aln.Length;
	bool downIsRev = aln.IsMateReverseStrand();
	int insSz = aln.InsertSize;
	int support = 1;

	if (aln.Position > aln.MatePosition)
	{
		std::swap(upSt, downSt);
		std::swap(upEn, downEn);
		std::swap(upIsRev, downIsRev);
		insSz *= -1;
		// do not add the same read twice
		support = 0;
	}

	// up should be forward and down should be reverse
	if (upIsRev || !downIsRev)
		return;

	std::vector<int> clipSizes, readPositions, genomePositions;
	bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
	// read should not have SC
	if (isSC)
		return;

	int upEnBound = upEn + bpRegion;
	int downStBound = downSt - bpRegion;
	DiscNode dn(upSt, upEnBound, downStBound, downEn, aln.RefID, support);
	nodes.emplace_back(dn);
}

void addSR(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<SplitNode> &nodes)
{
	if (!aln.HasTag("SA"))
		return;
	std::string tag;
	aln.GetTag("SA", tag);

	std::vector<std::string> mTags;
	std::vector<std::vector<std::string>> tags;
	splitMultipleTag(tag, mTags);
	splitTag(mTags, tags);

	std::string alnChr = ref.at(aln.RefID).RefName;
	for (std::vector<std::string> mt : tags)
	{
		std::string sChr = mt.at(0);
		int sPos = std::stoi(mt.at(1));
		std::string sCig = mt.at(3);

		if (alnChr != sChr)
			continue;

		int stAprx = aln.Position;
		int enAprx = sPos;
		int stBp, enBp;
		if (stAprx < enAprx)
		{
			stBp = aln.GetEndPosition();
			enBp = sPos;
		}
		else
		{
			stBp = sPos + parseCigar(sCig) - 1;
			enBp = aln.Position + 1;
		}
		int sz = enBp - stBp;
		if (sz > MAX_SPLIT_LEN || sz < MIN_SPLIT_LEN)
			continue;
		SplitNode sn(stBp, enBp, aln.RefID);
		nodes.emplace_back(sn);
	}
}

void addSC(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<SoftNode> &nodesUp, std::vector<SoftNode> &nodesDown)
{
	std::vector<int> clipSizes, readPositions, genomePositions;
	bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
	// read should have SC
	if (!isSC)
		return;

	bool down = false;
	if (aln.Position > aln.MatePosition)
	{
		down = true;
	}

	int scIdx;
	bool scAtRight = false;
	// softclips are either at beginning or end or both
	// sssnnn cigar:[X S .*] down bp
	// nnnsss cigar:[.* X S] up bp
	if (clipSizes.size() == 1)
	{
		scIdx = 0;
		if (aln.CigarData.front().Type == 'S')
			scAtRight = false;
		else
			scAtRight = true;
	}
	else if (clipSizes.size() == 2)
	{
		if (clipSizes.front() > clipSizes.back())
		{
			scAtRight = false;
			scIdx = 0;
		}
		else
		{
			scAtRight = true;
			scIdx = 1;
		}
	}

	int scPos = genomePositions[scIdx];

	std::string seq = aln.QueryBases;
	int seqLen = seq.size();

	std::string refName = ref.at(aln.RefID).RefName;
	std::string scSeq, nscSeq;
	int start = aln.Position, end = aln.Position + aln.AlignedBases.size();
	SoftNode sn(scPos, start, end, aln.RefID, refName, seq, aln.Name, nscSeq, down, scAtRight);

	if (sn.scAtRight)
		nodesUp.emplace_back(sn);
	else
		nodesDown.emplace_back(sn);
}

bool inExclude(long pos, std::vector<std::pair<long, long>> v)
{
	int low = 0, high = v.size() - 1;
	int it = -1;
	while (low <= high)
	{
		int mid = (low + high) >> 1;
		if (pos < v[mid].first)
			high = mid - 1;
		else
		{
			it = mid;
			low = mid + 1;
		}
	}
	if (it != -1 && v[it].first <= pos && pos <= v[it].second)
	{
		return true;
	}
	return false;
}

int getIdx(int pos, const std::vector<SoftCluster> &clusters)
{
	int ret = -1;
	int low = 0, high = clusters.size() - 1;
	while (low <= high)
	{
		int mid = (low + high) >> 1;
		if (clusters.at(mid).info.scPos > pos)
		{
			ret = mid;
			high = mid - 1;
		}
		else
			low = mid + 1;
	}
	return ret;
}

bool overlapCC(const SoftCluster &c1, const SoftCluster &c2)
{
	int sum = 0;
	int maxScore = 0;
	for (SoftNode n1 : c1.nodes)
	{
		for (SoftNode n2 : c2.nodes)
		{
			int here = alignSeq(n1.seq, n2.seq);
			sum += here;
		}
	}
	double avgScore = 1.0 * sum / (c1.nodes.size() * c2.nodes.size());
	return avgScore >= intraAlignScore;
}

void mergeSC(const std::vector<SoftCluster> &c, std::vector<SoftCluster> &mergedC)
{
	for (int i = c.size() - 1; i >= 0; i--)
	{
		bool found = false;
		if (!mergedC.empty() && mergedC.back().info.scPos - c.at(i).info.scPos <= MIN_SC_DISTANCE)
		{
			if (overlapCC(mergedC.back(), c.at(i)))
			{
				found = true;
				for (SoftNode n : c.at(i).nodes)
					mergedC.back().nodes.emplace_back(n);
			}
		}
		if (!found)
		{
			mergedC.emplace_back(c.at(i));
		}
	}
}

void readInput(BamTools::BamReader &br, const BamTools::RefVector &ref, const int &bpRegion, const bool verbose, const std::map<std::string, std::vector<std::pair<long, long>>> &exRegions, std::vector<DiscCluster> &discClusters, std::vector<SplitCluster> &srClusters, std::vector<SoftCluster> &scClustersUp, std::vector<SoftCluster> &scClustersDown)
{
	std::vector<DiscNode> discNodes;
	std::vector<SplitNode> srNodes;
	std::vector<SoftNode> scNodesUp, scNodesDown;

	BamTools::BamAlignment aln;
	while (br.GetNextAlignment(aln))
	{
		if (aln.MapQuality < MIN_MAP_QUAL)
			continue;
		// both reads should be mapped
		if (!(aln.IsMapped() && aln.IsMateMapped()))
			continue;
		// both pairs should be on same chr
		if (aln.RefID != aln.MateRefID)
			continue;

		if (exRegions.find(ref.at(aln.RefID).RefName) != exRegions.end())
		{
			if (inExclude(aln.Position, exRegions.at(ref.at(aln.RefID).RefName)))
				continue;
		}

		addDisc(aln, bpRegion, discNodes);
		addSR(aln, ref, srNodes);
		addSC(aln, ref, scNodesUp, scNodesDown);
	}

	if (debug)
	{
		std::cerr << "Disc nodes: " << discNodes.size() << '\n';
		std::cerr << "SR nodes: " << srNodes.size() << '\n';
		std::cerr << "SC nodes: " << scNodesUp.size() + scNodesDown.size() << '\n';
	}

	std::vector<SoftCluster> clustersUp, clustersDown, mergedUp, mergedDown;

	clusterDisc(discNodes, discClusters);
	clusterSR(srNodes, srClusters);
	clusterSC(scNodesUp, clustersUp);
	clusterSC(scNodesDown, clustersDown);

	std::sort(clustersUp.begin(), clustersUp.end(), SoftClusterCmp);
	std::sort(clustersDown.begin(), clustersDown.end(), SoftClusterCmp);

	mergeSC(clustersUp, mergedUp);
	mergeSC(clustersDown, mergedDown);

	for (auto c : mergedUp)
	{
		if (c.nodes.size() >= MIN_SC_CLUSTER_SUPPORT)
			scClustersUp.emplace_back(c);
	}
	for (auto c : mergedDown)
	{
		if (c.nodes.size() >= MIN_SC_CLUSTER_SUPPORT)
			scClustersDown.emplace_back(c);
	}
	std::sort(scClustersUp.begin(), scClustersUp.end(), SoftClusterCmp);
	std::sort(scClustersDown.begin(), scClustersDown.end(), SoftClusterCmp);

	if (debug)
	{
		std::cerr << "Disc clusters: " << discClusters.size() << '\n';
		std::cerr << "SR clusters: " << srClusters.size() << '\n';
		std::cerr << "SC clusters: " << scClustersUp.size() + scClustersDown.size() << '\n';
	}
}

void smallDeletionsSC(const BamTools::RefVector &ref, const std::vector<SoftCluster> &scClustersUp, const std::vector<SoftCluster> &scClustersDown, const std::string &folderPath)
{
	std::vector<OutNode> output;
	for (SoftCluster scUp : scClustersUp)
	{
		int stBound = getIdx(scUp.info.scPos + MIN_DEL_LEN, scClustersDown);

		if (stBound == -1)
			continue;

		int bestIdx = -1, bestSup = 0;
		for (int i = stBound; i < scClustersDown.size() - 1; ++i)
		{
			if (scClustersDown.at(i).info.scPos - scUp.info.scPos > MAX_DEL_LEN)
				break;
			if (scClustersDown.at(i).info.scPos - scUp.info.scPos < MIN_DEL_LEN)
				continue;
			if (overlapCC(scUp, scClustersDown.at(i)))
			{
				int curSup = scUp.nodes.size() + scClustersDown.at(i).nodes.size();
				if (curSup > bestSup)
				{
					bestSup = curSup;
					bestIdx = i;
				}
			}
		}
		if (bestIdx != -1)
		{
			OutNode on(ref.at(scUp.info.refID).RefName, scUp.info.scPos, scClustersDown.at(bestIdx).info.scPos, scClustersDown.at(bestIdx).info.scPos - scUp.info.scPos + 1, 0, 0, scUp.nodes.size() + scClustersDown.at(bestIdx).nodes.size());
			output.emplace_back(on);
		}
	}
	parseOutputN(folderPath, output, 1);
}

void largeDeletionsSC(const BamTools::RefVector &ref, const std::vector<DiscCluster> &discClusters, std::vector<SoftCluster> &scClustersUp, const std::vector<SoftCluster> &scClustersDown, const std::string &folderPath)
{
	std::vector<OutNode> output, impreciseOutput;
	for (DiscCluster dc : discClusters)
	{
		std::vector<SoftCluster> up, down;
		int upStBound = getIdx(dc.info.upStart, scClustersUp);
		for (int i = upStBound; i < scClustersUp.size(); ++i)
		{
			if (scClustersUp.at(i).info.scPos > dc.info.upEnd)
				break;
			up.emplace_back(scClustersUp.at(i));
		}
		int downStBound = getIdx(dc.info.downStart, scClustersDown);
		for (int i = downStBound; i < scClustersDown.size(); ++i)
		{
			if (scClustersDown.at(i).info.scPos > dc.info.downEnd)
				break;
			down.emplace_back(scClustersDown.at(i));
		}

		int bestSup = 0;
		SoftCluster bestUp, bestDown;
		for (SoftCluster sc1 : up)
		{
			for (SoftCluster sc2 : down)
			{
				int p1 = sc1.info.scPos, p2 = sc2.info.scPos;
				if (p1 > p2 || p2 - p1 < MIN_LARGE_DEL)
					continue;
				if (overlapCC(sc1, sc2))
				{
					int curSup = sc1.nodes.size() + sc2.nodes.size();
					if (curSup > bestSup)
					{
						bestSup = curSup;
						bestUp = sc1;
						bestDown = sc2;
					}
				}
			}
		}
		if (bestSup != 0)
		{
			OutNode on(ref.at(bestUp.info.refID).RefName, bestUp.info.scPos, bestDown.info.scPos, bestDown.info.scPos - bestUp.info.scPos + 1, dc.nodes.size(), 0, bestUp.nodes.size() + bestDown.nodes.size());
			output.emplace_back(on);
		}
		else
		{
			OutNode on(ref.at(dc.info.refID).RefName, dc.info.upStart, dc.info.downEnd, dc.info.downEnd - dc.info.upStart + 1, dc.nodes.size(), 0, 0);
			impreciseOutput.emplace_back(on);
		}
	}
	parseOutputN(folderPath, output, 0);
	parseOutputN(folderPath, impreciseOutput, 3);
}

void openExcludeRegions(const std::string fPath, std::map<std::string, std::vector<std::pair<long, long>>> &exRegions)
{
	std::string exChr;
	long exSt, exEn, exCount;
	double exCov;

	std::fstream fs;
	fs.open(fPath, std::fstream::in);
	if (!fs.is_open())
	{
		std::cerr << "Could not find exclude regions file at: " << fPath << '\n';
	}

	while (fs >> exChr >> exSt >> exEn >> exCount >> exCov)
	{
		exRegions[exChr].emplace_back(std::make_pair(exSt, exEn));
	}

	fs.close();

	// sort vector for each chr
	for (auto reg : exRegions)
	{
		std::sort(exRegions.at(reg.first).begin(), exRegions.at(reg.first).end());
	}
}

void openInputIntervals(const std::string fPath, std::vector<std::tuple<std::string, int, int>> &inp)
{
	std::fstream fs;
	fs.open(fPath, std::fstream::in);
	if (!fs.is_open())
	{
		std::cerr << "Could not find input regions file at: " << fPath << '\n';
		return;
	}

	std::string chr;
	int st, en;

	while (fs >> chr >> st >> en)
	{
		inp.emplace_back(std::make_tuple(chr, st, en));
	}
	fs.close();
}

void parallelProcess(const std::tuple<std::string, int, int> &region, const BamTools::RefVector &ref, const std::map<std::string, int> &revRefMap, const int bpRegion, const bool verbose, const std::map<std::string, std::vector<std::pair<long, long>>> &exRegions, const std::string &inpFilePath, const std::string &outFolderPath)
{
	BamTools::BamReader br;
	if (!openInput(inpFilePath, br))
	{
		return;
	}

	std::string chr = std::get<0>(region);
	int regSt = std::get<1>(region);
	int regEn = std::get<2>(region);
	int regRefId = revRefMap.at(chr);

	br.SetRegion(regRefId, regSt, regRefId, regEn);

	std::vector<DiscCluster> discClusters;
	std::vector<SplitCluster> srClusters;
	std::vector<SoftCluster> scClustersUp, scClustersDown;

	readInput(br, ref, bpRegion, verbose, exRegions, discClusters, srClusters, scClustersUp, scClustersDown);

	deletionsSR(ref, discClusters, srClusters, inpFilePath, outFolderPath);

	smallDeletionsSC(ref, scClustersUp, scClustersDown, outFolderPath);

	largeDeletionsSC(ref, discClusters, scClustersUp, scClustersDown, outFolderPath);

	br.Close();
}

void processInput(const std::string inpFilePath, const int mean, const int stdDev, const std::string outFolderPath, const bool verbose, const unsigned int threads)
{
	std::map<std::string, std::vector<std::pair<long, long>>> exRegions;
	std::string excludeFilePath = outFolderPath + excludeFileName;
	openExcludeRegions(excludeFilePath, exRegions);

	std::vector<std::tuple<std::string, int, int>> inpRegions;
	std::string inpIntervalsFilePath = outFolderPath + inpIntervalFileName;
	openInputIntervals(inpIntervalsFilePath, inpRegions);

	BamTools::BamReader br;
	if (!openInput(inpFilePath, br))
	{
		return;
	}

	BamTools::RefVector ref = br.GetReferenceData();

	br.Close();

	std::map<std::string, int> revRefMap;

	int refCo = 0;
	for (auto d : ref)
	{
		revRefMap[d.RefName] = refCo;
		refCo++;
	}

	for (int i = 0; i < 4; i++)
		headerOutput(outFolderPath, i);

	int bpRegion = mean + (3 * stdDev);

	transwarp::parallel executor{threads};

	std::vector<std::shared_ptr<transwarp::task<void>>> tasks;

	for (auto inpRegion : inpRegions)
	{
		auto regionTask = transwarp::make_value_task(inpRegion);
		auto refTask = transwarp::make_value_task(ref);
		auto mapTask = transwarp::make_value_task(revRefMap);
		auto bpRegionTask = transwarp::make_value_task(bpRegion);
		auto verboseTask = transwarp::make_value_task(verbose);
		auto exRegionsTask = transwarp::make_value_task(exRegions);
		auto inpFileTask = transwarp::make_value_task(inpFilePath);
		auto outFolderTask = transwarp::make_value_task(outFolderPath);

		auto processTask = transwarp::make_task(transwarp::consume, parallelProcess, regionTask, refTask, mapTask, bpRegionTask, verboseTask, exRegionsTask, inpFileTask, outFolderTask);

		tasks.emplace_back(processTask);
	}

	for (auto task : tasks)
	{
		task->schedule(executor);
	}

	for (auto task : tasks)
	{
		task->get();
	}
}
