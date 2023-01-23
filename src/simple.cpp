#include <bits/stdc++.h>

#include "cgshop2023_core/cpp_instance.hpp"
#include "cgshop2023_core/verify.hpp"
//#include "draw_solution.hpp"
#include "globals.hpp"
#include "localsearch.hpp"
#include "triangulation.hpp"
//#include <CGAL/draw_polygon_with_holes_2.h>

using namespace cgshop2023;
using namespace std;

void use_threads(vector<string> todos, size_t num_threads,
								 function<void(string)> process_file) {
	vector<thread> pool;
	mutex mtx;
	size_t i = 0;
	for (size_t t = 0; t < num_threads; ++t) {
		pool.emplace_back([&]() {
			while (true) {
				string filename;
				{
					lock_guard<mutex> l(mtx);
					if (i >= todos.size())
						break;
					filename = todos[i++];
					cout << "thread " << t << " processing file " << (i + 1) << "/"
							 << todos.size() << ": " << filename << endl;
				}
				process_file(filename);
			}
		});
	}
	for (size_t t = 0; t < num_threads; ++t) {
		pool[t].join();
	}
}

int main(int argc, char* argv[]) {
	bool orderBySize = false;
	bool init = false;
	size_t num_threads = 1;
	bool randomize = false;
	bool localsearch = false;
	size_t removal_attempts = 0;
	bool minimize = false;
	size_t replacement_choices = 0;
	bool conflict_optimize = false;
	unsigned max_attempts = 1000;
	unsigned max_iters = 10000;
	bool only_verify = false;
	for (int i = 1; i < argc; ++i) {
		string cur(argv[i]);
		auto eq = [&](const auto& a) { return cur == string(a); };
		auto next = [&]() { return string(argv[++i]); };
		if (eq("-h") || eq("--help")) {
			cerr << "Example usage (local search w/ minimize): ls instances | "
							"build/simple --order-by-size "
							"--localsearch --minimize --randomize --threads 3 "
							"--removal-attempts 100 --replacement-choices 100"
					 << endl;
			cerr << "Example usage (conflict optimizer): ls instances | build/simple "
							"--order-by-size --co --threads 3 --max-attempts 100 --max-iters "
							"100"
					 << endl;
			cerr << "Example usage (verify): ls instances | build/simple "
							"--order-by-size --verify --threads 3"
					 << endl;
		} else if (eq("--order-by-size"))
			orderBySize = true;
		else if (eq("--init"))
			init = true;
		else if (eq("--threads"))
			num_threads = stoi(next());
		else if (eq("--randomize"))
			randomize = true;
		else if (eq("--localsearch"))
			localsearch = true;
		else if (eq("--max-iters"))
			max_iters = stoi(next());
		else if (eq("--max-attempts"))
			max_attempts = stoi(next());
		else if (eq("--co"))
			conflict_optimize = true;
		else if (eq("--removal-attempts"))
			removal_attempts = stoi(next());
		else if (eq("--minimize"))
			minimize = true;
		else if (eq("--verify"))
			only_verify = true;
		else if (eq("--replacement-choices"))
			replacement_choices = stoi(next());
		else if (eq("--verbose") || eq("-v"))
			SET_VERBOSE(true);
		else {
			cerr << "Unknown command-line option: " << cur << endl;
		}
	}

	string filename;
	vector<string> files;
	while (cin >> filename) {
		files.push_back(remove_ext(filename));
	}

	if (orderBySize) {
		cerr << "Sorting inputs by total size\n";
		vector<pair<int, string>> files_sized;
		for (auto& filename : files) {
			Instance inst = Instance::read_file(filename);
			size_t size = inst.polygon().outer_boundary().size();
			for (const auto& hole : inst.polygon().holes())
				size += hole.size();
			files_sized.emplace_back(size, filename);
		}
		sort(files_sized.begin(), files_sized.end());
		files.clear();
		for (auto& [_, filename] : files_sized) {
			cout << filename << endl;
			files.push_back(filename);
		}
		cerr << "Finished sorting\n";
	}

	if (init) {
		use_threads(files, num_threads, [&](string filename) {
			Instance inst = Instance::read_file(filename);
			// Solution oldsol = Solution::read_file(filename);
			// size_t original_size = oldsol.size();
			size_t original_size = 30;
			Solution sol = basicTriangulation(inst);
			if (localsearch) {
				try_remove_all(inst, sol, randomize, removal_attempts, minimize,
											 replacement_choices);
			}
			if (!sol.write_if_better(inst, filename)) {
				cerr << "Did not see improvement to " << filename
						 << " (previous:" << original_size << ", new:" << sol.size() << ")"
						 << endl;
			}
		});
	} else if (localsearch) {
		use_threads(files, num_threads, [&](string filename) {
			Instance inst = Instance::read_file(filename);
			Solution sol = Solution::read_file(filename);
			size_t original_size = sol.size();
			try_remove_all(inst, sol, randomize, removal_attempts, minimize,
										 replacement_choices);
			if (!sol.write_if_better(inst, filename)) {
				cerr << "Did not see improvement to " << filename
						 << " (previous:" << original_size << ", new:" << sol.size() << ")"
						 << endl;
			}
		});
	} else if (conflict_optimize) {
		use_threads(files, num_threads, [&](string filename) {
			Instance inst = Instance::read_file(filename);
			Solution sol = Solution::read_file(filename);
			size_t original_size = sol.size();
			run_co(inst, sol, max_attempts, max_iters);
			if (!sol.write_if_better(inst, filename)) {
				cerr << "Did not see improvement to " << filename
						 << " (previous:" << original_size << ", new:" << sol.size() << ")"
						 << endl;
			}
		});
	} else if (only_verify) {
		Instance inst = Instance::read_file(filename);
		Solution sol = Solution::read_file(filename);
		SolutionVerifier sv(&inst, &sol);
		cerr << "Verify result: " << (sv.verify() ? "success" : "failed (invalid)")
				 << endl;
		if (sv.error_message().has_value()) {
			cerr << "Error message: " << sv.error_message().value() << endl;
		}
	}

	// old
	/*
	for (auto& filename : files) {
		auto fileloc = filename;
		ifstream ifs(fileloc);
		string out_name;
		Instance inst = Instance::read(ifs, out_name);
		cerr << "Successfully read " << out_name << endl;
		// CGAL::draw(inst.polygon());
		Solution sol = basicTriangulation(inst);
		cerr << "Successfully computed basic triangulation of " << out_name
				 << ", it has " << sol.polygons().size() << " triangles." << endl;
		// basicTriangulation(inst);
		// draw(inst, sol);
		try_remove_all(inst, sol);
		SolutionVerifier sv(&inst, &sol);
		cerr << "Verify result: " << (sv.verify() ? "success" : "failed (invalid)")
				 << endl;
		if (sv.error_message().has_value()) {
			cerr << "Error message: " << sv.error_message().value() << endl;
		} else {
			sol.write(cout, out_name);
		}
	}
	*/
}
