#include <dirent.h>

#include "Misc.hpp"

// args api
#include "args.hxx"

bool parse(int argc, char const *argv[], std::string &inpInsFile, std::string &outFolderPath)
{
	args::ArgumentParser parser("This program identifies insertions in NGS data");
	args::HelpFlag help(parser, "help", "Help menu", {'h', "help"});

	args::Group groupIO(parser, "Input Output", args::Group::Validators::All);
	args::ValueFlag<std::string> inpFileName(groupIO, "inpFile", "Input File", {'i', "inp"});
	args::ValueFlag<std::string> outFolderName(groupIO, "outFolder", "Output Folder", {'o', "out"});

	try
	{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help)
	{
		std::cout << parser;
		return false;
	}
	catch (args::ParseError e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return false;
	}
	catch (args::ValidationError e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return false;
	}

	inpInsFile = args::get(inpFileName);
	outFolderPath = args::get(outFolderName);
	if (outFolderPath.back() != '/')
	{
		outFolderPath += "/";
	}
	std::string tmpInsFolder = outFolderPath + "tmp/ins/";

	DIR *dir;
	dir = opendir(outFolderPath.c_str());
	if (dir != NULL)
	{
		closedir(dir);

		DIR *dirIns;
		dirIns = opendir(tmpInsFolder.c_str());
		if (dirIns != NULL)
		{
			closedir(dirIns);
			return true;
		}
	}
	std::cerr << "ERROR: Couldn't open output folder: " << outFolderPath << " or " << tmpInsFolder << '\n';
	return false;
}

void checkDirectIns(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<InsOutput> &out)
{
	std::string seq = aln.QueryBases;
	std::string refName = ref.at(aln.RefID).RefName;

	int pos = aln.Position; // pos on ref
	int cur = 0;			// pos on read

	for (int i = 0; i < aln.CigarData.size(); ++i)
	{
		if (aln.CigarData.at(i).Type == 'M')
		{
			pos += aln.CigarData.at(i).Length;
			cur += aln.CigarData.at(i).Length;
		}
		else if (aln.CigarData.at(i).Type == 'S')
		{
			cur += aln.CigarData.at(i).Length;
			continue;
		}
		else if (aln.CigarData.at(i).Type == 'D')
		{
			pos += aln.CigarData.at(i).Length;
			continue;
		}
		else if (aln.CigarData.at(i).Type == 'I')
		{
			if (aln.CigarData.at(i).Length >= 50)
			{
				std::string insSeq = seq.substr(cur, aln.CigarData.at(i).Length);
				out.emplace_back(InsOutput(refName, pos, insSeq, true, -1, aln.Name));
			}
			cur += aln.CigarData.at(i).Length;
		}
	}
}

void checkSC(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<InsOutput> &out)
{
	std::string seq = aln.QueryBases;
	std::string refName = ref.at(aln.RefID).RefName;

	std::vector<int> clipSizes, readPositions, genomePositions;
	bool isSC = aln.GetSoftClips(clipSizes, readPositions, genomePositions);
	// read should have SC
	if (!isSC)
		return;

	int scIdx = -1;
	bool scAtRight = false;

	if (clipSizes.size() == 1)
	{
		if (clipSizes.front() < 50)
			return;

		if (aln.CigarData.front().Type == 'S')
		{
			std::string insSeq = seq.substr(0, aln.CigarData.front().Length);
			out.emplace_back(InsOutput(refName, genomePositions.front(), insSeq, false, 0, aln.Name));
		}
		else if (aln.CigarData.back().Type == 'S')
		{
			int scLen = aln.CigarData.back().Length;
			std::string insSeq = seq.substr(seq.size() - scLen);
			out.emplace_back(InsOutput(refName, genomePositions.front(), insSeq, false, 1, aln.Name));
		}
	}
	else if (clipSizes.size() == 2)
	{
		if (clipSizes.front() > clipSizes.back())
		{
			if (clipSizes.front() < 50)
				return;
			std::string insSeq = seq.substr(0, aln.CigarData.front().Length);
			out.emplace_back(InsOutput(refName, genomePositions.front(), insSeq, false, 0, aln.Name));
		}
		else
		{
			if (clipSizes.back() < 50)
				return;
			int scLen = aln.CigarData.back().Length;
			std::string insSeq = seq.substr(seq.size() - scLen);
			out.emplace_back(InsOutput(refName, genomePositions.back(), insSeq, false, 1, aln.Name));
		}
	}
	if (scIdx == -1)
		return;
}

void checkSupplementary(BamTools::BamAlignment &aln, const BamTools::RefVector &ref, std::vector<InsOutput> &out)
{
	std::string refName = ref.at(aln.RefID).RefName;
	// 2048 is set
	if (aln.AlignmentFlag & (1 << 11))
	{
		for (int i = 0; i < out.size(); ++i)
		{
			InsOutput o = out.at(i);
			if (o.readName != aln.Name)
				continue;
			if (o.chr != refName)
				continue;
			if (aln.CigarData.front().Type == 'H' && aln.CigarData.back().Type == 'H')
				continue;

			// suffix aligned to ref
			if (aln.CigarData.front().Type == 'H' && o.loc == 1)
			{
				if (abs(o.pos - aln.Position) > 10)
					return;
				int alnLen = aln.QueryBases.size();
				std::string prvInsSeq = o.seq;
				int removeLen = prvInsSeq.size() - alnLen;
				if (removeLen < 0)
					continue;
				out.at(i).seq = prvInsSeq.substr(0, removeLen);
				out.at(i).ty = true;
				return;
			}

			// prefix aligned to ref
			if (aln.CigarData.back().Type == 'H' && o.loc == 0)
			{
				if (abs(o.pos - aln.GetEndPosition()) > 10)
					return;
				int removeLen = aln.QueryBases.size();
				std::string prvInsSeq = o.seq;
				if (removeLen < 0 || removeLen >= prvInsSeq.size())
					continue;
				out.at(i).seq = prvInsSeq.substr(removeLen);
				out.at(i).ty = true;
				return;
			}
		}
	}
}

void mergePartialPredictions(std::vector<InsOutput> &out)
{
	std::vector<InsOutput> ret;
	std::vector<int> done(out.size(), false);
	for (int i = 0; i < out.size(); ++i)
	{
		bool f = false;
		if (done.at(i))
			continue;
		for (int j = i + 1; j < out.size(); ++j)
		{
			if (done.at(j))
				continue;
			if (out.at(i).chr != out.at(j).chr)
				continue;
			if (abs(out.at(i).pos - out.at(j).pos) > 5)
				continue;
			if (out.at(i).ty || out.at(j).ty)
				continue;
			if ((out.at(i).loc == 0 && out.at(j).loc == 1) || (out.at(i).loc == 1 && out.at(j).loc == 0))
			{
				f = true;
				done.at(j) = 1;
				ret.emplace_back(out.at(i));
				ret.back().ty = true;
				if (out.at(i).loc == 0)
					ret.back().seqDown = out.at(i).seq;
				else
					ret.back().seqDown = out.at(j).seq;

				if (out.at(i).loc == 1)
					ret.back().seqUp = out.at(i).seq;
				else
					ret.back().seqUp = out.at(j).seq;

				ret.back().seq = "";
				break;
			}
		}
		if (!f)
			ret.emplace_back(out.at(i));
	}
	out.clear();
	for (int i = 0; i < ret.size(); ++i)
		out.emplace_back(ret.at(i));
}

void writeOutput(const std::vector<InsOutput> &out, const std::string &outFolderPath)
{
	std::ofstream ofs;
	std::string outFile = outFolderPath + "tmp/insertions.vcf";
	ofs.open(outFile, std::ofstream::out);

	for (InsOutput o : out)
	{
		ofs << o.chr << '\t' << o.pos << '\t' << 0 << '\t' << "N" << '\t' << "<INS>" << '\t' << "." << '\t' << "PASS" << '\t' << "SVTYPE=INS;SVLEN=" << o.seq.size() << ";";
		if (o.ty)
		{
			if (o.seq != "")
				ofs << "SEQ=" << o.seq << ";";
			else
				ofs << "SEQUP=" << o.seqUp << ";"
					<< "SEQDOWN=" << o.seqDown << ";IMPRECISE;";
		}
		else
		{
			if (o.loc == 0)
				ofs << "SEQDOWN=" << o.seq << ";";
			else if (o.loc == 1)
				ofs << "SEQUP=" << o.seq << ";";
			ofs << "IMPRECISE;";
		}
		ofs << "\tGT\t./.\n";
	}

	ofs.close();
}

void readInsInput(const std::string &inpInsFile, const std::string &outFolderPath)
{
	BamTools::BamReader br;
	if (!openInput(inpInsFile, br))
	{
		return;
	}

	BamTools::RefVector ref = br.GetReferenceData();
	BamTools::BamAlignment aln;

	std::vector<InsOutput> out;

	while (br.GetNextAlignment(aln))
	{
		// alignment quality check
		if (aln.MapQuality < MIN_MAP_QUAL)
			continue;
		// contig should be mapped
		if (!(aln.IsMapped()))
			continue;

		checkDirectIns(aln, ref, out);

		checkSC(aln, ref, out);

		checkSupplementary(aln, ref, out);
	}

	br.Close();

	mergePartialPredictions(out);

	writeOutput(out, outFolderPath);
}

int main(int argc, char const *argv[])
{
	std::string inpInsFile, outFolderPath;
	if (!parse(argc, argv, inpInsFile, outFolderPath))
		return EXIT_FAILURE;

	std::cerr << "[Step5] Insertions start\n";
	std::cerr << "        Input file: " << inpInsFile << '\n';
	std::cerr << "        Output folder: " << outFolderPath << '\n';

	readInsInput(inpInsFile, outFolderPath);

	std::cerr << "[Step5] Insertions end\n";

	return EXIT_SUCCESS;
}