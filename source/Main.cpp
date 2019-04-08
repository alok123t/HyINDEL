
#include "Parser.hpp"
#include "Detect.hpp"

int main(int argc, char const *argv[])
{
	std::chrono::steady_clock::time_point timeStart = std::chrono::steady_clock::now();

	ArgsParams ap;

	if (!parseArgs(argc, argv, ap))
	{
		return EXIT_FAILURE;
	}

	std::cerr << "[Step2] Deletions start\n";

	if (ap.verbose)
	{
		std::cerr << "Insert size: " << ap.insSz << '\n';
		std::cerr << "Standard deviation: " << ap.stdDev << '\n';
		std::cerr << "Read length: " << ap.readLen << '\n';
		std::cerr << "Coverage: " << ap.cov << '\n';
		std::cerr << "Input file: " << ap.inpFilePath << '\n';
		std::cerr << "Output folder: " << ap.outFolderPath << '\n';
		std::cerr << "Threads: " << ap.threads << '\n';
	}

	processInput(ap.inpFilePath, ap.insSz, ap.stdDev, ap.readLen, ap.cov, ap.outFolderPath, ap.verbose, ap.threads);

	std::chrono::steady_clock::time_point timeEnd = std::chrono::steady_clock::now();
	if (ap.verbose)
	{
		std::cerr << "Time taken: " << std::setprecision(1)
				  << (std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeStart).count())
				  << " secs" << '\n';
	}

	std::cerr << "[Step2] Deletions end\n";

	return EXIT_SUCCESS;
}
