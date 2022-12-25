#include <bits/stdc++.h>

#include "cgshop2023_core/cpp_instance.hpp"
#include "cgshop2023_core/verify.hpp"
#include <CGAL/draw_polygon_with_holes_2.h>

using namespace cgshop2023;
using namespace std;

int main() {
	string filename;
	vector<string> files;
	while (cin >> filename) {
		files.push_back(filename);
	}

	for (auto& filename : files) {
		auto fileloc = filename;
		ifstream ifs(fileloc);
		string out_name;
		Instance inst = Instance::read(ifs, out_name);
		cerr << "Successfully read " << out_name << endl;
		CGAL::draw(inst.polygon());
	}
}
