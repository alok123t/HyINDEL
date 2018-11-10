
#include "transwarp.h"

#include "ArgsParser.hpp"
#include "InputReader.hpp"

int main(int argc, char const *argv[]) {

	std::chrono::steady_clock::time_point timeStart = std::chrono::steady_clock::now();

	ArgsParams ap;
	
	if (!parseArgs(argc, argv, &ap)) return EXIT_FAILURE;

	if (ap.verbose) {
		std::cerr << "Insert size: " << ap.insSz << '\n';
		std::cerr << "Standard deviation: " << ap.stdDev << '\n';

		std::cerr << "Input files: " << '\n';
		for (std::string fp: ap.filePaths) {
			std::cerr << '\t' << fp << '\n';
		}

		std::cerr << "Threads: " <<  ap.threads << '\n';
	}

	transwarp::parallel executor{ap.threads};

	std::vector<std::shared_ptr<transwarp::task<void>>> tasks;

	for (std::string fp: ap.filePaths) {
		auto fpTask = transwarp::make_value_task(fp);
		auto meanTask = transwarp::make_value_task(ap.insSz);
		auto stdDevTask = transwarp::make_value_task(ap.stdDev);
		auto inputTask = transwarp::make_task(transwarp::consume, processInput, fpTask, meanTask, stdDevTask);
		
		tasks.emplace_back(inputTask);
	}

	for (auto task: tasks) {
		task->schedule(executor);
	}

	for (auto task: tasks) {
		task->get();
	}

	std::chrono::steady_clock::time_point timeEnd = std::chrono::steady_clock::now();

	std::cerr << "Time taken: " << std::setprecision(1) << (std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeStart).count()) << " secs" << '\n';

	return EXIT_SUCCESS;
}
