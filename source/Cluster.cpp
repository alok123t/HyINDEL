
#include "ArgsParser.hpp"
#include "InputReader.hpp"

int main(int argc, char const *argv[])
{
	std::chrono::steady_clock::time_point timeStart = std::chrono::steady_clock::now();

	ArgsParams ap;

	if (!parseArgs(argc, argv, ap))
	{
		return EXIT_FAILURE;
	}

	if (ap.verbose)
	{
		std::cerr << "Insert size: " << ap.insSz << '\n';
		std::cerr << "Standard deviation: " << ap.stdDev << '\n';
		std::cerr << "Read length: " << ap.readLen << '\n';

		std::cerr << "Input files: " << ap.inpFilePath << '\n';

		std::cerr << "Output folder: " << ap.outFolderPath << '\n';

		std::cerr << "Threads: " << ap.threads << '\n';
	}

	processInput(ap.inpFilePath, ap.insSz, ap.stdDev, ap.readLen, ap.outFolderPath, ap.verbose, ap.threads);

	std::chrono::steady_clock::time_point timeEnd = std::chrono::steady_clock::now();
	if (ap.verbose)
	{
		std::cerr << "Time taken: " << std::setprecision(1)
				  << (std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeStart).count())
				  << " secs" << '\n';
	}

	return EXIT_SUCCESS;
}
