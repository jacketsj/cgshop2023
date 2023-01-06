#pragma once

#include "cgshop2023_core/cpp_instance.hpp"
#include "cgshop2023_core/verify.hpp"
#include <CGAL/ch_graham_andrew.h>

using namespace cgshop2023;
using namespace std;

using Kernel = CGAL::Epeck;
using Point = CGAL::Point_2<Kernel>;
using Polygon = CGAL::Polygon_with_holes_2<Kernel>;
using SimplePolygon = CGAL::Polygon_2<Kernel>;

SimplePolygon greedy_expand(Instance& inst, SimplePolygon poly,
														const vector<Point>& desired_coverage,
														bool& allCovered) {
	allCovered = true;
	auto outerArea = area(inst.polygon());
	for (const Point& p : desired_coverage) {
		// decide wether to add p to poly
		// look at what the new polygon would be by computing convex hull (very lazy
		// method) this could be improved to O(log n) with binary search (and the
		// new edges could be singled out)
		vector<Point> points(poly.vertices_begin(), poly.vertices_end());
		points.push_back(p);
		// get the convex hull
		vector<Point> chull;
		CGAL::ch_graham_andrew(points.begin(), points.end(),
													 std::back_inserter(chull));
		SimplePolygon newPoly(chull.begin(), chull.end());

		// is it inside the polygon
		std::vector<Polygon> coverage = {};
		std::vector<Polygon> polys = {inst.polygon(), Polygon(newPoly)};
		CGAL::join(polys.begin(), polys.end(), std::back_inserter(coverage));
		if (coverage.size() == 1 && area(coverage[0]) <= outerArea) {
			// it is, so add it
			poly = newPoly;
		} else
			allCovered = false;
	}

	return poly;
}

bool try_removal(Instance& inst, Solution& sol, size_t polygon_i) {
	vector<Point> desired_coverage(sol.polygons()[polygon_i].begin(),
																 sol.polygons()[polygon_i].end());
	// try to remove polygon_i
	// do it by trying to get polygons in sol to cover it's vertices
	bool succeeded = false;
	for (size_t i = 0; i < sol.polygons().size(); ++i) {
		if (i == polygon_i)
			continue;
		bool allCovered = false;
		SimplePolygon newPoly =
				greedy_expand(inst, sol.polygons()[i], desired_coverage, allCovered);
		if (allCovered) {
			sol.polygons_m()[i] = newPoly;
			succeeded = true;
			break;
		}
	}
	return succeeded;
}

void removal_if_possible(Instance& inst, Solution& sol, size_t polygon_i) {
	if (try_removal(inst, sol, polygon_i)) {
		swap(sol.polygons_m()[polygon_i],
				 sol.polygons_m()[sol.polygons().size() - 1]);
		sol.polygons_m().pop_back();
	}
}

void try_remove_all(Instance& inst, Solution& sol) {
	cerr << "Running try_remove_all on " << sol.polygons().size()
			 << " polygons\n";
	for (int polygon_i = sol.polygons().size() - 1; polygon_i >= 0; --polygon_i) {
		removal_if_possible(inst, sol, polygon_i);
	}
	cerr << "Finished running try_remove_all, now have " << sol.polygons().size()
			 << " polygons\n";
}
