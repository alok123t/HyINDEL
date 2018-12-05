#include "InputReader.hpp"

bool openInput(const std::string filePath, BamTools::BamReader &br)
{
	if (!br.Open(filePath))
	{
		std::cerr << "ERROR:BamTools Could not open " << filePath << '\n';
		return false;
	}
	std::string indexPath = filePath + indexExt;
	if (!br.OpenIndex(indexPath))
	{
		std::cerr << "ERROR:BamTools Could not open index " << indexPath << '\n';
		return false;
	}
	return true;
}

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

void largeDeletions(const std::string discPath, const std::string softPath, const int mean, const int stdDev, const std::string folderPath)
{
	BamTools::BamReader discBr;
	if (!openInput(discPath, discBr))
	{
		return;
	}

	// candidate large deletions
	std::vector<DiscCluster> clusters, filteredClusters;
	readDiscordant(discBr, mean, stdDev, clusters);

	for (int i = 0; i < clusters.size(); ++i)
	{
		DiscNode dr = clusters[i].info;
		if (dr.support > MIN_DISC_CLUSTER_SUPPORT && dr.downEnd - dr.upStart <= MAX_DISC_CLUSTER_DEL_LEN)
		{
			filteredClusters.emplace_back(clusters[i]);
		}
	}
	clusters.clear();
	std::cerr << "Large deletions cluster done\n";

	discBr.Close();

	// parse softclip around candidate regions
	std::vector<std::vector<std::string>> outputDelsLarge;
	DelsLarge(softPath, filteredClusters, outputDelsLarge);
	std::cerr << "Large deletions overlap done\n";

	// print deletion to file
	std::string outputFile = folderPath;
	getFileName(discPath, outputFile);
	parseOutput(outputFile, outputDelsLarge, 0);
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

void processInput(const std::string discPath, const std::string softPath, const int mean, const int stdDev, const std::string folderPath)
{
	largeDeletions(discPath, softPath, mean, stdDev, folderPath);

	smallDeletions(softPath, folderPath);
}
