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

void getFileName(const std::string &filePath, std::string &outFilePath)
{
	std::string fileName = filePath.substr(filePath.find_last_of('/') + 1);
	std::size_t dotLen = fileName.find_last_of('.');
	std::string fileNameNoExt = fileName.substr(0, dotLen);
	outFilePath += fileNameNoExt;
}

Range CreateRange(int refID1, int start1, int end1, bool isRev1, bool isSC1, int refID2, int start2, int end2, bool isRev2, bool isSC2, int support, bool isShort)
{
	Range ret;
	ret.refID1 = refID1;
	ret.start1 = start1;
	ret.end1 = end1;
	ret.isPlus1 = !isRev1;
	ret.isMinus1 = isRev1;
	ret.isSC1 = isSC1;
	ret.refID2 = refID2;
	ret.start2 = start2;
	ret.end2 = end2;
	ret.isPlus2 = !isRev2;
	ret.isMinus2 = isRev2;
	ret.isSC2 = isSC2;
	ret.support = support;
	ret.isShort = isShort;
	return ret;
}

bool OverlapRange(Range r1, Range r2)
{
	if (r1.refID1 == r2.refID1 && r1.refID2 == r2.refID2 && r1.start1 <= r2.end1 && r1.end1 >= r2.start1 && r1.start2 <= r2.end2 && r1.end2 >= r2.start2)
		return true;
	else
		return false;
}

bool OverlapOrientation(Range r1, Range r2)
{
	if (r1.isPlus1 == r2.isPlus1 && r1.isPlus2 == r2.isPlus2 && r1.isMinus1 == r2.isMinus1 && r1.isMinus2 == r2.isMinus2)
		return true;
	else
		return false;
}

void AddRange(std::vector<Range> &v, Range r)
{
	if (v.size() == 0)
	{
		v.push_back(r);
		return;
	}

	std::vector<Range>::reverse_iterator rit = v.rbegin();
	while (rit != v.rend())
	{
		//Diff chr
		if (rit->refID1 != r.refID1)
			break;
		//Doesn't overlap
		if ((rit->end2 < r.start1) || (rit->start1 > r.end2))
			break;

		if (OverlapRange(*rit, r) && OverlapOrientation(*rit, r))
		{
			//Update
			if (r.start1 > rit->start1)
				rit->start1 = r.start1;
			if (r.start2 > rit->start2)
				rit->start2 = r.start2;
			if (r.end1 < rit->end1)
				rit->end1 = r.end1;
			if (r.end2 < rit->end2)
				rit->end2 = r.end2;

			rit->isPlus1 |= r.isPlus1;
			rit->isPlus2 |= r.isPlus2;
			rit->isMinus1 |= r.isMinus1;
			rit->isMinus2 |= r.isMinus2;
			rit->support += r.support;
			//Large
			if (!rit->isShort)
			{
				rit->isSC1 &= r.isSC1;
				rit->isSC2 &= r.isSC2;
			} //Small
			else
			{
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

void readInput(std::string filePath, BamTools::BamReader &br, const int mean, const int stdDev, const std::string &outputFile)
{
	std::vector<Range> dels_large, dels_small;

	BamTools::BamAlignment aln;

	int refID1, refID2, start1, start2, end1, end2, insSize, support, readLength, centre, bound1, bound2, Threshold, bpRegion;
	bool isRev1, isRev2, isSC1, isSC2, isShort, isUp, at5;
	std::vector<int> clipSizes, readPositions, genomePositions;
	std::vector<int> down(5), up(5);
	Range tmpRange;

	bool nxtAlnExists = br.GetNextAlignmentCore(aln);
	while (nxtAlnExists)
	{
		//Parse
		if (aln.IsMapped() && aln.IsMateMapped())
		{
			clipSizes.clear();
			readPositions.clear();
			genomePositions.clear();
			//Large SV
			if (!aln.IsProperPair())
			{
				//First Read
				readLength = aln.Length;
				isShort = false;
				Threshold = bpRegion = mean + (3 * stdDev);
				if (aln.Position < aln.MatePosition)
				{
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
				} //Second Read
				else
				{
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
				if (!(isSC1 && isSC2))
				{
					if (refID1 == refID2 && insSize >= Threshold && !isRev1 && isRev2)
					{
						centre = std::ceil((start2 + end1) / 2);
						bound1 = std::min(centre, end1 + bpRegion);
						bound2 = std::max(centre, start2 - bpRegion);
						bound1 = end1 + bpRegion;
						bound2 = start2 - bpRegion;

						tmpRange = CreateRange(refID1, start1, bound1, isRev1, isSC1, refID2, bound2, end2, isRev2, isSC2, support, isShort);
						AddRange(dels_large, tmpRange);
					}
				}
			}
			//Small SV
			if (aln.IsProperPair())
			{
				if (aln.RefID == aln.MateRefID && aln.GetSoftClips(clipSizes, readPositions, genomePositions) && clipSizes.size() == 1 && clipSizes[0] >= minSCLen)
				{
					bpRegion = 5 * stdDev;
					isShort = true;
					support = 0;

					at5 = (aln.CigarData[0].Type == 'S');
					isUp = (aln.Position < aln.MatePosition);
					isRev1 = aln.IsReverseStrand();
					isRev2 = aln.IsMateReverseStrand();
					refID1 = refID2 = aln.RefID;

					//Going Downstream
					if (at5)
					{
						start1 = aln.Position - bpTol;
						end1 = aln.Position + bpTol;
					}
					else
					{
						start1 = aln.GetEndPosition() - bpTol;
						end1 = aln.GetEndPosition() + bpTol;
					}
					start2 = end1 + 1;
					end2 = start2 + bpRegion;
					down = {start1, end1, start2, end2};

					//Going Upstream
					if (at5)
					{
						start2 = aln.Position - bpTol;
						end2 = aln.Position + bpTol;
					}
					else
					{
						start2 = aln.GetEndPosition() - bpTol;
						end2 = aln.GetEndPosition() + bpTol;
					}
					end1 = start2 - 1;
					start1 = end1 - bpRegion;
					up = {start1, end1, start2, end2};

					//Upstream
					if (at5)
					{
						//Down bp1 ... sc/bp2
						//Up sc/bp1 ... bp2
						if ((isUp && !isRev1) || (!isUp && isRev1))
						{
							isSC1 = false;
							isSC2 = true;

							tmpRange = CreateRange(refID1, up[0], up[1], false, isSC1, refID2, up[2], up[3], true, isSC2, support, isShort);
							AddRange(dels_small, tmpRange);
						}
					}
					else
					{ //Downstream
						//Down bp1 ... bp2/sc
						//Up bp1/sc ... bp2
						if ((isUp && !isRev1) || (!isUp && isRev1))
						{
							isSC1 = true;
							isSC2 = false;

							tmpRange = CreateRange(refID1, down[0], down[1], false, isSC1, refID2, down[2], down[3], true, isSC2, support, isShort);
							AddRange(dels_small, tmpRange);
						}
					}
				}
			}
		}
		nxtAlnExists = br.GetNextAlignmentCore(aln);
	}

	std::cerr << "Candidate Large Ranges: " << dels_large.size() << '\n';
	std::cerr << "Candidate Small Ranges: " << dels_small.size() << '\n';

	std::vector<std::vector<std::string>> info_dels_large, info_dels_small;

	// DelsParse(filePath, dels_large, info_dels_large);
	// std::cerr << "Done Large deletions\n";
	// DelsParse(filePath, dels_small, info_dels_small);
	// std::cerr << "Done Small deletions\n";

	// parseOutput(outputFile, info_dels_large, info_dels_small);
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
			if (insSz < Threshold)
				continue;

			clipSizes.clear();
			readPositions.clear();
			genomePositions.clear();
			bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
			// read should not have SC
			if (isSC)
				continue;
			int upEnBound = upEn + bpRegion;
			int downStBound = downSt - bpRegion;
			DiscNode rng(upSt, upEnBound, downStBound, downEn, aln.RefID, support);
			clusterDiscRange(rng, clusters);
		}
	}
}

void processInput(const std::string discPath, const std::string softPath, const int mean, const int stdDev, const std::string folderPath)
{
	BamTools::BamReader discBr;

	if (!openInput(discPath, discBr))
	{
		return;
	}

	std::string outputFile = folderPath;
	getFileName(discPath, outputFile);

	// large deletions
	std::vector<DiscCluster> clusters, filteredClusters;
	readDiscordant(discBr, mean, stdDev, clusters);

	for (int i = 0; i < clusters.size(); ++i)
	{
		DiscNode dr = clusters[i].info;
		if (dr.support > 3 && dr.downEnd - dr.upStart <= 50000)
		{
			filteredClusters.emplace_back(clusters[i]);
		}
	}
	clusters.clear();
	std::cerr << "Candidate deletions: " << filteredClusters.size() << '\n';

	discBr.Close();

	std::vector<std::vector<std::string>> outputDelsLarge;
	DelsParse(softPath, filteredClusters, outputDelsLarge);
	std::cerr << "Done Large deletions\n";

	parseOutput(outputFile, outputDelsLarge, 0);

	// readInput(filePath, br, mean, stdDev, outputFile);
}