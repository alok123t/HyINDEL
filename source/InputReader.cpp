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

void readDiscordant(BamTools::BamReader &br, const int &Mean, const int &StdDev, std::vector<DiscCluster> &clusters)
{
	int Threshold = Mean + (3 * StdDev);
	int bpRegion = Mean + (3 * StdDev);

	BamTools::BamAlignment aln;
	std::vector<int> clipSizes, readPositions, genomePositions;

	while (br.GetNextAlignmentCore(aln))
	{
		// both reads should be mapped
		if (!(aln.IsMapped() && aln.IsMateMapped()))
			continue;

		if (!aln.IsProperPair())
		{
			// both pairs should be on same chr
			if (aln.RefID != aln.MateRefID)
				continue;

			int upSt = aln.Position;
			int upEn = aln.GetEndPosition();
			bool upIsRev = aln.IsReverseStrand();
			int downSt = aln.MatePosition;
			int downEn = downSt + aln.Length;
			bool downIsRev = aln.IsMateReverseStrand();
			int insSz = aln.InsertSize;
			int support = 1;
			// Second read
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
				continue;
			// only large dels
			// if (insSz < Threshold)
			// 	continue;

			clipSizes.clear();
			readPositions.clear();
			genomePositions.clear();
			bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
			// read should not have SC
			if (isSC)
				continue;
			int upEnBound = upEn + bpRegion;
			int downStBound = downSt - bpRegion;
			DiscNode dn(upSt, upEnBound, downStBound, downEn, aln.RefID, support);
			clusterDiscRange(dn, clusters);
		}
	}
}

bool overlapIntervals(int st1, int en1, int st2, int en2, double RO)
{
	int oSt = std::max(st1, st2);
	int oEn = std::min(en1, en2);

	double len1 = (en1 - st1 + 1) * 1.0;
	double len2 = (en2 - st2 + 1) * 1.0;
	double oLen = (oEn - oSt + 1) * 1.0;

	double ro1 = oLen / len1;
	double ro2 = oLen / len2;

	bool ret = ro1 >= RO && ro2 >= RO;
	// std::cout << st1 << ' ' << en1 << ' ' << st2 << ' ' << en2 << ' ' << ro1 << ' ' << ro2 << '\n';

	return ret;
}

void DelsLargeSplit(BamTools::RefVector &ref, const std::vector<DiscCluster> &discClusters, const std::vector<SplitCluster> &splitClusters, std::vector<int> &usedSR, std::vector<OutNode> &output)
{
	/*
	for (int i = 0; i < discClusters.size(); ++i)
	{
		DiscCluster dc = discClusters.at(i);
		if (dc.nodes.size() < 2)
			continue;
		for (int j = 0; j < splitClusters.size(); ++j)
		{
			if (usedSR.at(j))
				continue;
			SplitCluster sc = splitClusters.at(j);
			if (sc.nodes.size() < 2)
				continue;
			if (dc.info.refID != sc.info.refID)
				continue;
			if (dc.info.upStart <= sc.info.st && sc.info.en < dc.info.downEnd)
			{
				usedSR.at(j) = 1;
				OutNode on(ref.at(sc.info.refID).RefName, sc.info.st, sc.info.en, sc.info.en - sc.info.st + 1, dc.nodes.size(), sc.nodes.size(), 0);
				output.emplace_back(on);
				break;
			}
		}
	}
	*/
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

void largeDeletions(const std::vector<SplitCluster> &splitClusters, const std::string inpFilePath, const int mean, const int stdDev, const std::string folderPath)
{
	BamTools::BamReader br;
	if (!openInput(inpFilePath, br))
	{
		return;
	}

	BamTools::RefVector ref = br.GetReferenceData();

	// candidate large deletions
	std::vector<DiscCluster> clusters, filteredClusters;
	readDiscordant(br, mean, stdDev, clusters);

	std::string inpFname;
	getFileName(br.GetFilename(), inpFname);

	for (int i = 0; i < clusters.size(); ++i)
	{
		DiscNode dr = clusters[i].info;
		if (dr.support > MIN_DISC_CLUSTER_SUPPORT && dr.downEnd - dr.upStart <= MAX_DISC_CLUSTER_DEL_LEN)
		{
			filteredClusters.emplace_back(clusters[i]);
		}
	}
	std::cerr << inpFname << " Large deletions cluster done\n";

	br.Close();

	// parse softclip around candidate regions
	std::vector<std::vector<std::string>> outputDelsLarge;
	std::vector<OutNode> outputDelsSplitLarge;
	DelsLarge(inpFilePath, filteredClusters, outputDelsLarge);
	std::cerr << inpFname << " Large deletions overlap done\n";

	// DelsLargeSplit(ref, clusters, splitClusters, outputDelsSplitLarge);

	// print deletion to file
	std::string outputFile = folderPath;
	getFileName(inpFilePath, outputFile);
	parseOutput(outputFile, outputDelsLarge, 0);
	// parseOutputN(outputFile, outputDelsSplitLarge);
}

void large(BamTools::RefVector &ref, const std::vector<DiscCluster> &discClusters, const std::vector<SplitCluster> &splitClusters, const std::string inpFilePath, const std::string folderPath)
{
	std::vector<int> usedSR(splitClusters.size(), 0);
	std::vector<OutNode> outputDelsSplitLarge;
	DelsLargeSplit(ref, discClusters, splitClusters, usedSR, outputDelsSplitLarge);

	// print deletion to file
	std::string outputFile = folderPath;
	getFileName(inpFilePath, outputFile);
	parseOutputN(outputFile, outputDelsSplitLarge);
}

bool overlapSoft(const SoftNode &n, const SoftCluster c)
{
	if (n.scPos > c.info.scPos + 20)
		return false;

	std::string nLeft, nRight, curLeft, curRight;

	if (n.scAtRight)
	{
		nLeft = n.nscSeq;
		nRight = n.scSeq;
	}
	else
	{
		nLeft = n.scSeq;
		nRight = n.nscSeq;
	}

	for (int i = 0; i < c.nodes.size(); i++)
	{
		int diff = std::abs(n.scPos - c.nodes[i].scPos);
		if (diff > bpTol)
			return false;

		if (c.nodes.at(i).scAtRight)
		{
			curLeft = c.nodes.at(i).nscSeq;
			curRight = c.nodes.at(i).scSeq;
		}
		else
		{
			curLeft = c.nodes.at(i).scSeq;
			curRight = c.nodes.at(i).nscSeq;
		}
		if (!Align(nLeft, curLeft))
			return false;
		if (!Align(nRight, curRight))
			return false;
	}
	return true;
}

void clusterSoft(const SoftNode &n, std::vector<SoftCluster> &c)
{
	for (int i = c.size() - 1; i >= 0; i--)
	{
		// diff chr
		if (c[i].info.refID != n.refID)
			break;
		// too far apart, will not find a cluster
		if (n.scPos > c[i].info.scPos + 100)
			break;
		if (overlapSoft(n, c[i]))
		{
			c[i].nodes.emplace_back(n);
			if (c[i].nodes.size() > MAX_NODES_IN_CLUSTER)
			{
				c[i].nodes.clear();
			}
			return;
		}
	}
	// add new cluster
	SoftCluster sc(n);
	c.emplace_back(n);
}

void readSoft(BamTools::BamReader &br, std::vector<SoftNode> &nodes)
{
	int co = 0;
	// br.SetRegion(18, 10000000, 18, 20000000);
	BamTools::BamAlignment aln;
	std::vector<int> clipSizes, readPositions, genomePositions;
	std::unordered_set<int> ust;

	while (br.GetNextAlignmentCore(aln))
	{
		if (!(aln.IsMapped() && aln.IsMateMapped()))
			continue;

		// if (!aln.IsProperPair())
		// 	continue;
		clipSizes.clear();
		readPositions.clear();
		genomePositions.clear();
		bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
		// read should have SC
		if (!isSC)
			continue;

		if (aln.RefID != aln.MateRefID)
			continue;

		bool down = false;
		bool upIsRev = aln.IsReverseStrand();
		bool downIsRev = aln.IsMateReverseStrand();
		if (aln.Position > aln.MatePosition)
		{
			down = true;
			std::swap(upIsRev, downIsRev);
		}
		// up should be forward and down should be reverse
		if (upIsRev || !downIsRev)
			continue;

		bool scAtRight = false;
		int scIdx = 0;
		// softclips are either at beginning or end or both
		// sssnnn cigar:[X S .*] down bp
		// nnnsss cigar:[.* X S] up bp
		if (clipSizes.size() == 1)
		{
			scIdx = 0;
			if (clipSizes.front() < minSCLen)
				continue;
			if (aln.CigarData.front().Type == 'S')
				scAtRight = false;
			else
				scAtRight = true;
		}
		else if (clipSizes.size() == 2)
		{
			if (clipSizes.front() < 10 && clipSizes.back() < 10)
				continue;

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

		int scLen = clipSizes[scIdx];
		int scPos = genomePositions[scIdx];

		aln.BuildCharData();
		std::string seq = aln.QueryBases;
		int seqLen = seq.size();
		std::string seqQual = aln.Qualities;

		int trimLen = trimSeq(seq, seqQual, down);
		seqLen -= trimLen;
		scLen -= trimLen;

		if (scLen < minSCLen)
			continue;

		std::string refName = br.GetReferenceData()[aln.RefID].RefName;
		std::string scSeq, nscSeq;
		SoftNode sn(scPos, aln.RefID, refName, scSeq, nscSeq, down, scAtRight);

		int bId = scPos / 1000;

		if (ust.find(bId) != ust.end())
			continue;

		if (!scAtRight)
		{
			sn.scSeq = seq.substr(0, scLen);
			sn.nscSeq = seq.substr(scLen, seqLen - scLen);

			// bktDown[bId].push_back(sn);
		}
		else
		{
			// if del exists in cigar
			int rPos = readPositions[scIdx];
			rPos -= getDelSize(aln);

			sn.nscSeq = seq.substr(0, seqLen - scLen);
			sn.scSeq = seq.substr(rPos, scLen);

			// bktUp[bId].push_back(sn);
		}

		nodes.emplace_back(sn);

		// if (bktUp[bId].size() > 1000 || bktDown[bId].size() > 1000)
		// {
		// 	ust.insert(bId);
		// 	bktUp[bId].clear();
		// 	bktDown[bId].clear();
		// }
	}
	std::sort(nodes.begin(), nodes.end(), SoftCmp);
}

void localizeSoft(const std::vector<SoftNode> &nodes, std::vector<std::vector<SoftCluster>> &clustersUp, std::vector<std::vector<SoftCluster>> &clustersDown)
{
	std::vector<std::vector<SoftNode>> bktNodes;
	bktNodes.resize(MAX_SZ);
	for (int i = 0; i < nodes.size(); i++)
	{
		int bId = nodes.at(i).scPos / BUCKET_SIZE;
		if (bktNodes.at(bId).size() < MAX_NODES_IN_BUCKET)
		{
			bktNodes[bId].emplace_back(nodes.at(i));
		}
	}

	for (int i = 0; i < bktNodes.size(); i++)
	{
		std::chrono::steady_clock::time_point timeStart = std::chrono::steady_clock::now();
		for (int j = 0; j < bktNodes.at(i).size(); j++)
		{
			if (bktNodes.at(i).at(j).scAtRight)
			{
				clusterSoft(bktNodes.at(i).at(j), clustersUp.at(i));
			}
			else
			{
				clusterSoft(bktNodes.at(i).at(j), clustersDown.at(i));
			}
		}
		std::chrono::steady_clock::time_point timeEnd = std::chrono::steady_clock::now();
		// std::cout << std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeStart).count() << '\n';
	}
}

void smallDeletions(const std::string softPath, const std::string folderPath)
{
	BamTools::BamReader softBr;
	if (!openInput(softPath, softBr))
	{
		return;
	}

	std::string outputFile = folderPath;
	getFileName(softPath, outputFile);

	std::vector<SoftNode> nodes;
	readSoft(softBr, nodes);

	softBr.Close();

	std::vector<std::vector<SoftCluster>> clustersUp(MAX_SZ), clustersDown(MAX_SZ);
	localizeSoft(nodes, clustersUp, clustersDown);
	std::cerr << "Small deletions cluster done\n";

	std::vector<std::vector<std::string>> outputDelsSmall;
	DelsSmall(clustersUp, clustersDown, outputDelsSmall);
	std::cerr << "Small deletions overlap done\n";

	// print deletion to file
	parseOutput(outputFile, outputDelsSmall, 1);
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

void readSplit(BamTools::BamReader &br, BamTools::RefVector &ref, std::vector<SplitNode> &nodes)
{
	BamTools::BamAlignment aln;

	int co = 0;
	while (br.GetNextAlignmentCore(aln))
	{
		// both reads should be mapped
		if (!(aln.IsMapped() && aln.IsMateMapped()))
			continue;
		// both pairs should be on same chr
		if (aln.RefID != aln.MateRefID)
			continue;
		aln.BuildCharData();
		// aln should have supplementary bit set
		if (!aln.HasTag("SA"))
			continue;
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

void splitReads(BamTools::BamReader &br, BamTools::RefVector &ref, std::vector<SplitCluster> &clusters)
{
	std::vector<SplitNode> nodes;
	readSplit(br, ref, nodes);

	clusterSR(nodes, clusters);

	std::sort(clusters.begin(), clusters.end(), SplitCmp);

	std::cerr << "Split reads cluster done\n";
}

void clusterDisc(const std::vector<DiscNode> &nodes, std::vector<DiscCluster> &clusters)
{
	for (DiscNode dn : nodes)
	{
		clusterDiscRange(dn, clusters);
	}
}

void clusterSC(const std::vector<SoftNode> &nodes, std::vector<std::vector<SoftCluster>> &clustersUp, std::vector<std::vector<SoftCluster>> &clustersDown)
{
	localizeSoft(nodes, clustersUp, clustersDown);
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

void addSC(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<SoftNode> &nodes)
{
	std::vector<int> clipSizes, readPositions, genomePositions;
	bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
	// read should have SC
	if (!isSC)
		return;

	bool down = false;
	bool upIsRev = aln.IsReverseStrand();
	bool downIsRev = aln.IsMateReverseStrand();
	if (aln.Position > aln.MatePosition)
	{
		down = true;
		std::swap(upIsRev, downIsRev);
	}
	// up should be forward and down should be reverse
	if (upIsRev || !downIsRev)
		return;

	bool scAtRight = false;
	int scIdx = 0;
	// softclips are either at beginning or end or both
	// sssnnn cigar:[X S .*] down bp
	// nnnsss cigar:[.* X S] up bp
	if (clipSizes.size() == 1)
	{
		scIdx = 0;
		if (clipSizes.front() < minSCLen)
			return;
		if (aln.CigarData.front().Type == 'S')
			scAtRight = false;
		else
			scAtRight = true;
	}
	else if (clipSizes.size() == 2)
	{
		if (clipSizes.front() < 10 && clipSizes.back() < 10)
			return;

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

	int scLen = clipSizes[scIdx];
	int scPos = genomePositions[scIdx];

	std::string seq = aln.QueryBases;
	int seqLen = seq.size();
	std::string seqQual = aln.Qualities;

	if (aln.MapQuality <= MIN_ALN_QUAL)
		return;

	int trimLen = trimSeq(seq, seqQual, down);
	seqLen -= trimLen;
	scLen -= trimLen;

	if (scLen < minSCLen)
		return;

	std::string refName = ref.at(aln.RefID).RefName;
	std::string scSeq, nscSeq;
	SoftNode sn(scPos, aln.RefID, refName, scSeq, nscSeq, down, scAtRight);

	if (!scAtRight)
	{
		sn.scSeq = seq.substr(0, scLen);
		sn.nscSeq = seq.substr(scLen, seqLen - scLen);
	}
	else
	{
		// if del exists in cigar
		int rPos = readPositions[scIdx];
		rPos -= getDelSize(aln);

		sn.nscSeq = seq.substr(0, seqLen - scLen);
		sn.scSeq = seq.substr(rPos, scLen);
	}

	nodes.emplace_back(sn);
}

void readInput(BamTools::BamReader &br, const BamTools::RefVector &ref, const int &bpRegion, const bool verbose, std::vector<DiscCluster> &discClusters, std::vector<SplitCluster> &srClusters, std::vector<std::vector<SoftCluster>> &scClustersUp, std::vector<std::vector<SoftCluster>> &scClustersDown)
{
	std::vector<DiscNode> discNodes;
	std::vector<SplitNode> srNodes;
	std::vector<SoftNode> scNodes;

	BamTools::BamAlignment aln;
	while (br.GetNextAlignment(aln))
	{
		// both reads should be mapped
		if (!(aln.IsMapped() && aln.IsMateMapped()))
			continue;
		// both pairs should be on same chr
		if (aln.RefID != aln.MateRefID)
			continue;

		addDisc(aln, bpRegion, discNodes);
		addSR(aln, ref, srNodes);
		addSC(aln, ref, scNodes);
	}

	if (debug)
	{
		std::cerr << "Disc nodes: " << discNodes.size() << '\n';
		std::cerr << "SR nodes: " << srNodes.size() << '\n';
		std::cerr << "SC nodes: " << scNodes.size() << '\n';
	}

	clusterDisc(discNodes, discClusters);
	clusterSR(srNodes, srClusters);
	// clusterSC(scNodes, scClustersUp, scClustersDown);

	if (debug)
	{
		std::cerr << "Disc clusters: " << discClusters.size() << '\n';
		std::cerr << "SR clusters: " << srClusters.size() << '\n';
	}
	if (verbose)
	{
		std::string inpFname;
		getFileName(br.GetFilename(), inpFname);
		std::cerr << inpFname << " Cluster done\n";
	}
}

void processInput(const std::string inpFilePath, const int mean, const int stdDev, const std::string outFolderPath, const bool verbose)
{
	BamTools::BamReader br;
	if (!openInput(inpFilePath, br))
	{
		return;
	}

	BamTools::RefVector ref = br.GetReferenceData();

	int bpRegion = mean + (3 * stdDev);

	std::vector<DiscCluster> discClusters;
	std::vector<SplitCluster> srClusters;
	std::vector<std::vector<SoftCluster>> scClustersUp(MAX_SZ), scClustersDown(MAX_SZ);

	// readInput(br, ref, bpRegion, verbose, discClusters, srClusters, scClustersUp, scClustersDown);

	br.Close();

	// large(ref, discClusters, srClusters, inpFilePath, outFolderPath);
	largeDeletions(srClusters, inpFilePath, mean, stdDev, outFolderPath);

	// smallDeletions(inpFilePath, outFolderPath);
}
