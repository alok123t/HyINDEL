
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

		std::cerr << "Input files: " << '\n';
		for (int i = 0; i < ap.discFilePaths.size(); ++i)
		{
			std::cerr << '\t' << ap.discFilePaths[i] << '\t' << ap.softFilePaths[i] << '\n';
		}

		std::cerr << "Output folder: " << ap.outFolderPath << '\n';

		std::cerr << "Threads: " << ap.threads << '\n';
	}

	transwarp::parallel executor{ap.threads};

	std::vector<std::shared_ptr<transwarp::task<void>>> tasks;

	for (int i = 0; i < ap.discFilePaths.size(); ++i)
	{
		auto discTask = transwarp::make_value_task(ap.discFilePaths[i]);
		auto softTask = transwarp::make_value_task(ap.softFilePaths[i]);
		auto meanTask = transwarp::make_value_task(ap.insSz);
		auto stdDevTask = transwarp::make_value_task(ap.stdDev);
		auto outFolderTask = transwarp::make_value_task(ap.outFolderPath);
		auto inputTask = transwarp::make_task(transwarp::consume, processInput, discTask, softTask, meanTask, stdDevTask, outFolderTask);

		tasks.emplace_back(inputTask);
	}

	for (auto task : tasks)
	{
		task->schedule(executor);
	}

	for (auto task : tasks)
	{
		task->get();
	}

	std::chrono::steady_clock::time_point timeEnd = std::chrono::steady_clock::now();
	std::cerr << "Time taken: " << std::setprecision(1)
			  << (std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeStart).count())
			  << " secs" << '\n';

	return EXIT_SUCCESS;
}
