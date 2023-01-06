#include <bits/stdc++.h>

#include "cgshop2023_core/cpp_instance.hpp"
#include "cgshop2023_core/verify.hpp"
#include "draw_solution.hpp"
#include "localsearch.hpp"
#include "triangulation.hpp"
//#include <CGAL/draw_polygon_with_holes_2.h>

using namespace cgshop2023;
using namespace std;

int main() {
	string filename;
	vector<string> files;
	while (cin >> filename) {
		files.push_back(filename);
	}

	/*
	bool orderBySize = true;
	if (orderBySize) {
		cerr << "Sorting inputs by total size\n";
		vector<pair<int, string>> files_sized;
		for (auto& filename : files) {
			auto fileloc = filename;
			ifstream ifs(fileloc);
			string out_name;
			Instance inst = Instance::read(ifs, out_name);
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
		exit(0);
	}
	*/

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
}
