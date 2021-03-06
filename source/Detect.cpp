#include "Detect.hpp"

int INTRA_ALIGN_SCORE;
int MIN_SC_CLUSTER_SUPPORT;
int MIN_SPLIT_SUPPORT;
int MIN_DISC_CLUSTER_SUPPORT;

bool overlapDiscRange(const DiscNode &r, const std::vector<DiscNode> &nodes)
{
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

		bool ok = ro1 >= DISC_RO && ro2 >= DISC_RO;
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

void clusterDisc(const std::vector<DiscNode> &nodes, std::vector<DiscCluster> &clusters)
{
	std::vector<DiscCluster> discClusters;
	for (DiscNode dn : nodes)
	{
		clusterDiscRange(dn, discClusters);
	}
	for (DiscCluster dc : discClusters)
	{
		if (dc.info.support > MIN_DISC_CLUSTER_SUPPORT)
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
		m.at(0).at(j) = j * GAP_PENALTY;
	for (int i = 1; i < rows; ++i)
	{
		for (int j = 1; j < cols; ++j)
		{
			int sim = a.at(i - 1) == b.at(j - 1) ? MATCH_SCORE : MIS_MATCH_SCORE;
			m.at(i).at(j) = std::max(sim + m.at(i - 1).at(j - 1),
									 std::max(GAP_PENALTY + m.at(i - 1).at(j), GAP_PENALTY + m.at(i).at(j - 1)));
			maxScore = std::max(maxScore, m.at(i).at(j));
		}
	}
	return maxScore;
}

bool overlapSC(const SoftCluster &c, const SoftNode &n)
{
	for (SoftNode cn : c.nodes)
	{
		if (abs(n.scPos - cn.scPos) > MAX_SC_OVERLAP_DISTANCE)
			return false;
		// compare reads in sorted order, first seq should be more towards 5` than second seq
		double alnScore = 1.0 * alignSeq(cn.seq, n.seq);
		double o_len = 1.0 * (std::min(n.end, cn.end) - std::max(n.start, cn.start));
		if (alnScore < ALN_THRESHOLD * o_len)
			return false;
	}
	return true;
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
			if (overlapSC(clusters.at(i), n))
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
			clusters.emplace_back(sc);
		}
	}
}

void addDisc(BamTools::BamAlignment &aln, const int &bpRegion, std::vector<DiscNode> &nodes)
{
	// both reads should be mapped
	if (!aln.IsMapped() || !aln.IsMateMapped())
		return;
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

	// current read is 3` and its mate at 5`
	if (aln.Position > aln.MatePosition)
	{
		std::swap(upSt, downSt);
		std::swap(upEn, downEn);
		std::swap(upIsRev, downIsRev);
		insSz *= -1;
		// do not add the same read twice
		support = 0;
	}

	if (insSz >= MAX_LARGE_DEL_LEN)
		return;

	// up should be forward and down should be reverse
	if (upIsRev || !downIsRev)
		return;

	std::vector<int> clipSizes, readPositions, genomePositions;
	bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
	// read should not have SC
	if (isSC)
		return;

	DiscNode dn(upSt, upEn, downSt, downEn, aln.RefID, support);
	nodes.emplace_back(dn);
}

void addHC(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<SoftNode> &nodesUp, std::vector<SoftNode> &nodesDown)
{
	bool isSC = false;
	bool bothMapped = aln.IsMapped() & aln.IsMateMapped();
	if (aln.CigarData.back().Type == 'H')
	{
		int start = aln.Position, end = aln.GetEndPosition();
		SoftNode sn(end, start, end, aln.RefID, ref.at(aln.RefID).RefName, aln.QueryBases, aln.Name, true, isSC);
		nodesUp.emplace_back(sn);
	}
	else if (aln.CigarData.front().Type == 'H')
	{
		int start = aln.Position, end = aln.GetEndPosition();
		SoftNode sn(start, start, end, aln.RefID, ref.at(aln.RefID).RefName, aln.QueryBases, aln.Name, false, isSC);
		nodesDown.emplace_back(sn);
	}
}

void addSC(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<SoftNode> &nodesUp, std::vector<SoftNode> &nodesDown)
{
	std::vector<int> clipSizes, readPositions, genomePositions;
	bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
	// read should have SC
	if (!isSC)
		return;

	int scIdx = 0;
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
	std::string refName = ref.at(aln.RefID).RefName;
	bool bothMapped = aln.IsMapped() & aln.IsMateMapped();
	SoftNode sn(scPos, -1, -1, aln.RefID, refName, seq, aln.Name, scAtRight, isSC);

	if (sn.scAtRight)
	{
		sn.start = aln.Position;
		sn.end = aln.Position + aln.QueryBases.size();
		nodesUp.emplace_back(sn);
	}
	else
	{
		sn.start = aln.GetEndPosition() - aln.QueryBases.size();
		sn.end = aln.GetEndPosition();
		nodesDown.emplace_back(sn);
	}
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
		if (pos <= clusters.at(mid).info.scPos)
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
	std::vector<SoftNode> up(c1.nodes.begin(), c1.nodes.end());
	std::vector<SoftNode> down(c2.nodes.begin(), c2.nodes.end());

	// check if split read already present
	std::set<std::string> stUp, stDown;
	for (SoftNode sn : up)
	{
		if (!sn.isSC)
			stUp.insert(sn.readName);
	}
	for (SoftNode sn : down)
	{
		if (!sn.isSC)
			stDown.insert(sn.readName);
	}
	int srCo = 0;
	for (SoftNode sn : up)
	{
		if (sn.isSC)
		{
			if (stDown.find(sn.readName) != stDown.end())
				srCo += 1;
		}
	}
	for (SoftNode sn : down)
	{
		if (sn.isSC)
		{
			if (stUp.find(sn.readName) != stUp.end())
				srCo += 1;
		}
	}
	if (srCo >= MIN_SPLIT_SUPPORT)
		return true;

	std::vector<SoftNode> up3, down3;
	for (int i = up.size() - 1; i >= 0; i--)
	{
		if (up.at(i).isSC)
			up3.emplace_back(up.at(i));
		if (up3.size() >= 3)
			break;
	}
	for (int i = down.size() - 1; i >= 0; i--)
	{
		if (down.at(i).isSC)
			down3.emplace_back(down.at(i));
		if (down3.size() >= 3)
			break;
	}

	for (SoftNode n1 : up3)
	{
		for (SoftNode n2 : down3)
		{
			int here = alignSeq(n1.seq, n2.seq);
			sum += here;
		}
	}
	double avgScore = 1.0 * sum / (up3.size() + down3.size());
	return avgScore >= INTRA_ALIGN_SCORE;
}

void addDDel(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<DirectDelNode> &nodes)
{
	int pos = aln.Position;
	int delLen = 0;
	for (int i = 0; i < aln.CigarData.size(); ++i)
	{
		if (aln.CigarData.at(i).Type == 'I')
			continue;
		else if (aln.CigarData.at(i).Type == 'M')
			pos += aln.CigarData.at(i).Length;
		else if (aln.CigarData.at(i).Type == 'S')
			continue;
		else if (aln.CigarData.at(i).Type == 'D')
		{
			if (aln.CigarData.at(i).Length >= MIN_SMALL_DEL_LEN)
			{
				DirectDelNode dn(pos, aln.CigarData.at(i).Length, ref.at(aln.RefID).RefName);
				nodes.emplace_back(dn);
			}
		}
	}
}

bool checkRO(const DirectDelNode &n1, const DirectDelNode &n2)
{
	int st1 = n1.pos, en1 = n1.pos + n1.delLen;
	int st2 = n2.pos, en2 = n2.pos + n2.delLen;
	double overlap = 1.0 * (std::min(en1, en2) - std::max(st1, st2));
	double len1 = 1.0 * (en1 - st1);
	double len2 = 1.0 * (en2 - st2);
	double ro1 = overlap / len1;
	double ro2 = overlap / len2;
	return ro1 >= DDEL_RO && ro2 >= DDEL_RO;
}

void clusterDDel(const std::vector<DirectDelNode> &nodes, std::vector<DirectDelCluster> &clusters)
{
	for (DirectDelNode n : nodes)
	{
		bool found = false;
		for (int i = clusters.size() - 1; i >= 0; i--)
		{
			// diff chr
			if (clusters.at(i).info.refName != n.refName)
				break;
			if (checkRO(clusters.at(i).info, n))
			{
				// push new node
				clusters.at(i).nodes.emplace_back(n);
				found = true;
				break;
			}
		}
		// add new cluster
		if (!found)
		{
			DirectDelCluster dc(n);
			clusters.emplace_back(dc);
		}
	}
}

void directDeletions(const std::vector<DirectDelCluster> &clusters, const std::string &folderPath)
{
	std::vector<OutNode> output;

	for (auto dc : clusters)
	{
		if (dc.nodes.size() >= MIN_DDEL_SUPPORT)
		{
			OutNode on(dc.info.refName, dc.info.pos, dc.info.pos + dc.info.delLen - 1, 0, dc.nodes.size(), 0);
			output.emplace_back(on);
		}
	}
	parseOutput(folderPath, output, 1);
}

void outputFastq(BamTools::BamAlignment &aln, std::ofstream &ofs)
{
	std::string name = aln.Name;
	if (aln.IsPaired())
		name.append((aln.IsFirstMate() ? "/1" : "/2"));

	std::string qualities = aln.Qualities;
	std::string sequence = aln.QueryBases;
	if (aln.IsReverseStrand())
	{
		reverse(qualities.begin(), qualities.end());

		std::size_t seqLength = sequence.length();
		for (std::size_t i = 0; i < seqLength; ++i)
			sequence.replace(i, 1, 1, REVCOMP_LOOKUP[(int)sequence.at(i) - 65]);
		reverse(sequence.begin(), sequence.end());
	}

	ofs << '@' << name << '\n'
		<< sequence << '\n'
		<< '+' << '\n'
		<< qualities << '\n';
}

void writeReads(BamTools::BamReader &br, const BamTools::RefVector &ref, int refID, const std::string &refName, int st, int en, int bpRegion, const std::string &folderPath)
{
	std::ofstream ofs;
	std::string outFile = folderPath + "tmp/ins/reads.fastq";
	ofs.open(outFile, std::ofstream::out | std::ofstream::app);

	std::vector<SoftNode> insNodesUp, insNodesDown;
	std::vector<SoftCluster> insClustersUp, insClustersDown;

	br.SetRegion(refID, st - bpRegion, refID, en + bpRegion);
	BamTools::BamAlignment aln;
	while (br.GetNextAlignment(aln))
	{
		// OEA unmapped reads have mapping qual 0
		// current aln is unmapped, but it's mate lies in region, hence we can even parse cur aln
		if (!aln.IsMapped())
		{
			outputFastq(aln, ofs);
			continue;
		}

		if (aln.MapQuality < MIN_MAP_QUAL)
			continue;

		std::vector<int> clipSizes, readPositions, genomePositions;
		bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
		if (isSC)
		{
			outputFastq(aln, ofs);
			addSC(aln, ref, insNodesUp, insNodesDown);
		}
	}
	ofs.close();
}

void insertions(BamTools::BamReader &br, const BamTools::RefVector &ref, const int bpRegion, const std::vector<SoftCluster> &scClustersUp, const std::vector<SoftCluster> &scClustersDown, const std::string &folderPath)
{
	for (SoftCluster scUp : scClustersUp)
	{
		int stBound = getIdx(scUp.info.scPos + MIN_INS_DIFF_CLUSTERS, scClustersDown);
		if (stBound == -1)
			continue;

		int bestIdx = -1, bestSup = 0;
		for (int i = stBound; i < scClustersDown.size(); ++i)
		{
			if (scClustersDown.at(i).info.scPos - scUp.info.scPos > MAX_INS_DIFF_CLUSTERS)
				break;
			if (scClustersDown.at(i).info.scPos - scUp.info.scPos < MIN_INS_DIFF_CLUSTERS)
				continue;
			int curSup = scUp.nodes.size() + scClustersDown.at(i).nodes.size();
			if (curSup > bestSup)
			{
				bestSup = curSup;
				bestIdx = i;
			}
		}
		if (bestIdx != -1)
		{
			int st = scUp.info.scPos;
			int en = scClustersDown.at(bestIdx).info.scPos;
			writeReads(br, ref, scUp.info.refID, ref.at(scUp.info.refID).RefName, st, en, bpRegion, folderPath);
		}
	}
}

void readInput(BamTools::BamReader &br, const BamTools::RefVector &ref, const int &bpRegion, const bool verbose, const std::map<std::string, std::vector<std::pair<long, long>>> &exRegions, std::vector<DiscCluster> &discClusters, std::vector<SoftCluster> &scClustersUp, std::vector<SoftCluster> &scClustersDown, std::vector<DirectDelCluster> &dDelClusters, const std::string &outFolderPath)
{
	std::vector<DiscNode> discNodes;
	std::vector<SoftNode> scNodesUp, scNodesDown;
	std::vector<DirectDelNode> dDelNodes;

	BamTools::BamAlignment aln;
	while (br.GetNextAlignment(aln))
	{
		// alignment quality check
		if (aln.MapQuality < MIN_MAP_QUAL)
			continue;
		// read should be mapped
		if (!(aln.IsMapped()))
			continue;
		// both pairs should be on same chr
		if (aln.RefID != aln.MateRefID)
			continue;
		// not part of excluded region
		if (exRegions.find(ref.at(aln.RefID).RefName) != exRegions.end())
		{
			if (inExclude(aln.Position, exRegions.at(ref.at(aln.RefID).RefName)))
				continue;
		}

		addDisc(aln, bpRegion, discNodes);
		addHC(aln, ref, scNodesUp, scNodesDown);
		addSC(aln, ref, scNodesUp, scNodesDown);
		addDDel(aln, ref, dDelNodes);
	}

	std::sort(scNodesUp.begin(), scNodesUp.end(), SoftCmp);
	std::sort(scNodesDown.begin(), scNodesDown.end(), SoftCmp);

	if (debug)
	{
		std::cerr << "Disc nodes: " << discNodes.size() << '\n';
		std::cerr << "SC nodes: " << scNodesUp.size() + scNodesDown.size() << '\n';
		std::cerr << "Direct Del nodes: " << dDelNodes.size() << '\n';
	}

	std::vector<SoftCluster> clustersUp, clustersDown;

	clusterDDel(dDelNodes, dDelClusters);

	clusterDisc(discNodes, discClusters);

	clusterSC(scNodesUp, clustersUp);
	clusterSC(scNodesDown, clustersDown);

	std::sort(clustersUp.begin(), clustersUp.end(), SoftClusterCmp);
	std::sort(clustersDown.begin(), clustersDown.end(), SoftClusterCmp);

	for (auto c : clustersUp)
	{
		if (c.nodes.size() > MIN_SC_CLUSTER_SUPPORT)
			scClustersUp.emplace_back(c);
	}
	for (auto c : clustersDown)
	{
		if (c.nodes.size() > MIN_SC_CLUSTER_SUPPORT)
			scClustersDown.emplace_back(c);
	}

	std::sort(scClustersUp.begin(), scClustersUp.end(), SoftClusterCmp);
	std::sort(scClustersDown.begin(), scClustersDown.end(), SoftClusterCmp);

	if (debug)
	{
		std::cerr << "Disc clusters: " << discClusters.size() << '\n';
		std::cerr << "SC clusters: " << scClustersUp.size() + scClustersDown.size() << '\n';
		std::cerr << "Direct Del clusters: " << dDelClusters.size() << '\n';
	}
}

bool checkMergeDels(const std::string &chr1, int st1, int en1, const std::string &chr2, int st2, int en2, double RO)
{
	if (chr1 != chr2)
		return false;
	bool closeEnough = (abs(st1 - st2) <= DEL_MERGE_DIS) & (abs(en1 - en2) <= DEL_MERGE_DIS);
	double overlap = 1.0 * (std::min(en1, en2) - std::max(st1, st2));
	double len1 = 1.0 * (en1 - st1);
	double len2 = 1.0 * (en2 - st2);
	double ro1 = overlap / len1;
	double ro2 = overlap / len2;
	return ((ro1 >= RO && ro2 >= RO) | closeEnough);
}

void removeCommonPreds(std::vector<OutNode> &small, std::vector<OutNode> &large)
{
	std::vector<OutNode> filterSmall;
	for (OutNode onS : small)
	{
		bool found = false;
		for (OutNode onL : large)
		{
			if (checkMergeDels(onS.chr, onS.st, onS.en, onL.chr, onL.st, onL.en, DEL_MERGE_RO))
			{
				found = true;
				break;
			}
		}
		if (!found)
			filterSmall.emplace_back(onS);
	}
	small.clear();
	for (OutNode on : filterSmall)
		small.emplace_back(on);
}

void removeMultiplePreds(std::vector<OutNode> &events)
{
	std::vector<OutNode> filterEvents(events);
	std::vector<bool> toRemove(events.size(), false);
	for (int i = 0; i < events.size(); i++)
	{
		OutNode onI = events.at(i);
		int supI = onI.supSC + onI.supSR;
		for (int j = i + 1; j < events.size(); j++)
		{
			if (toRemove.at(j))
				continue;
			OutNode onJ = events.at(j);
			int supJ = onJ.supSC + onJ.supSR;
			if (checkMergeDels(onI.chr, onI.st, onI.en, onJ.chr, onJ.st, onJ.en, DEL_MERGE_RO))
			{
				if (supI >= supJ)
					toRemove.at(j) = true;
				else
					toRemove.at(i) = true;
			}
		}
	}
	events.clear();
	for (int i = 0; i < filterEvents.size(); ++i)
	{
		if (!toRemove.at(i))
			events.emplace_back(filterEvents.at(i));
	}
}

std::pair<int, int> calcSup(const std::vector<SoftNode> &sc1, const std::vector<SoftNode> &sc2)
{
	std::set<std::string> st1, st2;
	for (SoftNode sn : sc1)
	{
		if (!sn.isSC)
			st1.insert(sn.readName);
	}
	for (SoftNode sn : sc2)
	{
		if (!sn.isSC)
			st2.insert(sn.readName);
	}
	int srCo = 0, scCo = 0;
	for (SoftNode sn : sc1)
	{
		if (sn.isSC)
		{
			if (st2.find(sn.readName) != st2.end())
				srCo += 1;
			else
				scCo += 1;
		}
	}
	for (SoftNode sn : sc2)
	{
		if (sn.isSC)
		{
			if (st1.find(sn.readName) != st1.end())
				srCo += 1;
			else
				scCo += 1;
		}
	}
	return std::make_pair(srCo, scCo);
}

void smallDeletions(const BamTools::RefVector &ref, const std::vector<SoftCluster> &scClustersUp, const std::vector<SoftCluster> &scClustersDown, std::vector<OutNode> &smallOutput)
{
	for (SoftCluster scUp : scClustersUp)
	{
		int stBound = getIdx(scUp.info.scPos + MIN_SMALL_DEL_LEN, scClustersDown);

		if (stBound == -1)
			continue;

		int bestIdx = -1, bestSup = 0;
		for (int i = stBound; i < scClustersDown.size(); ++i)
		{
			if (scClustersDown.at(i).info.scPos - scUp.info.scPos > MAX_SMALL_DEL_LEN)
				break;
			if (scClustersDown.at(i).info.scPos - scUp.info.scPos < MIN_SMALL_DEL_LEN)
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
			std::pair<int, int> sup = calcSup(scUp.nodes, scClustersDown.at(bestIdx).nodes);
			OutNode on(ref.at(scUp.info.refID).RefName, scUp.info.scPos, scClustersDown.at(bestIdx).info.scPos, 0, sup.first, sup.second);
			smallOutput.emplace_back(on);
		}
	}
}

void largeDeletions(const BamTools::RefVector &ref, const std::vector<DiscCluster> &discClusters, std::vector<SoftCluster> &scClustersUp, const std::vector<SoftCluster> &scClustersDown, const int bpRegion, std::vector<OutNode> &largeOutput, std::vector<OutNode> &impreciseOutput)
{
	for (DiscCluster dc : discClusters)
	{
		std::vector<SoftCluster> up, down;
		int upStBound = getIdx(dc.info.upStart, scClustersUp);
		for (int i = upStBound; i < scClustersUp.size(); ++i)
		{
			if (scClustersUp.at(i).info.scPos > dc.info.upEnd + bpRegion)
				break;
			up.emplace_back(scClustersUp.at(i));
		}

		int downStBound = getIdx(dc.info.downStart - bpRegion, scClustersDown);
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
				if (p1 > p2)
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
			std::pair<int, int> sup = calcSup(bestUp.nodes, bestDown.nodes);
			OutNode on(ref.at(bestUp.info.refID).RefName, bestUp.info.scPos, bestDown.info.scPos, dc.info.support, sup.first, sup.second);
			largeOutput.emplace_back(on);
		}
		else
		{
			if (dc.info.upEnd < dc.info.downStart)
			{
				OutNode on(ref.at(dc.info.refID).RefName, dc.info.upEnd, dc.info.downStart, dc.info.support, 0, 0);
				impreciseOutput.emplace_back(on);
			}
		}
	}
}

void openExcludeRegions(const std::string &fPath, std::map<std::string, std::vector<std::pair<long, long>>> &exRegions)
{
	std::string exChr;
	long exSt, exEn;

	std::fstream fs;
	fs.open(fPath, std::fstream::in);
	if (!fs.is_open())
	{
		std::cerr << "Could not find exclude regions file at: " << fPath << '\n';
		return;
	}
	else
	{
		std::cerr << "        Exclude regions file: " << fPath << '\n';
	}

	while (fs >> exChr >> exSt >> exEn)
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

void openInputIntervals(const std::string &fPath, std::vector<std::tuple<std::string, int, int>> &inp)
{
	std::fstream fs;
	fs.open(fPath, std::fstream::in);
	if (!fs.is_open())
	{
		std::cerr << "Could not find input regions file at: " << fPath << '\n';
		return;
	}
	else
	{
		std::cerr << "        Input regions file: " << fPath << '\n';
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
	std::vector<SoftCluster> scClustersUp, scClustersDown;
	std::vector<DirectDelCluster> dDelClusters;

	readInput(br, ref, bpRegion, verbose, exRegions, discClusters, scClustersUp, scClustersDown, dDelClusters, outFolderPath);

	std::vector<OutNode> smallOutput, largeOutput, impreciseOutput;

	insertions(br, ref, bpRegion, scClustersUp, scClustersDown, outFolderPath);

	directDeletions(dDelClusters, outFolderPath);

	smallDeletions(ref, scClustersUp, scClustersDown, smallOutput);

	largeDeletions(ref, discClusters, scClustersUp, scClustersDown, bpRegion, largeOutput, impreciseOutput);

	removeCommonPreds(smallOutput, largeOutput);
	removeMultiplePreds(smallOutput);
	removeMultiplePreds(largeOutput);

	parseOutput(outFolderPath, smallOutput, 0);
	parseOutput(outFolderPath, largeOutput, 1);
	parseOutput(outFolderPath, impreciseOutput, 2);

	br.Close();
}

void processInput(const std::string &inpFilePath, const int insSz, const int stdDev, const int readLen, const double coverage, const std::string &outFolderPath, const bool verbose, const unsigned int threads)
{
	INTRA_ALIGN_SCORE = readLen / 2;
	MIN_SC_CLUSTER_SUPPORT = MIN_SPLIT_SUPPORT = MIN_DISC_CLUSTER_SUPPORT = std::max(1, (int)(coverage / 10));

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

	for (int i = 0; i < 3; i++)
		headerOutput(outFolderPath, i);

	int bpRegion = insSz + (3 * stdDev);

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
