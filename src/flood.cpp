#include "flood.hpp"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using Kernel = CGAL::Epeck;
using Point = CGAL::Point_2<Kernel>;
using Line = CGAL::Line_2<Kernel>;
using Ray = CGAL::Ray_2<Kernel>;
using Polygon = CGAL::Polygon_with_holes_2<Kernel>;
using SimplePolygon = CGAL::Polygon_2<Kernel>;

void Flooder::flood_raycast(const Point& p, const Point& dir) {
	// TODO
	// compute arrangement (or maintain)
	// look at all the lines the ray intersects
	// filter for intersection points inside the polygon
	// binary search along the intersection points
	// at each one, check if the expansion to include the point (i.e. convex
	// hull w/ that point) is inside container
	vector<Point> relevant_points = container_points;
}

void Flooder::run() {
	Point interior = cur[0];
	// TODO better orders
	for (const auto& dir : container.outer_boundary()) {
		flood_raycast(interior, dir);
	}
	for (const auto& hole : container.holes())
		for (const auto& dir : hole) {
			flood_raycast(interior, dir);
		}
	auto original_area = area(cur);
	auto cur_area = original_area;
	do {
		original_area = area(cur);
		vector<pair<Point, Point>> to_run;
		for (size_t i = 0; i < cur.size(); ++i)
			to_run.emplace_back(cur[i], cur[(i + 1) % cur.size()]);
		for (auto& [p1, p2] : to_run) {
			flood_raycast(p1, p2);
			flood_raycast(p2, p1);
		}
		cur_area = area(cur);
	} while (cur_area > original_area);
}

SimplePolygon flood(const Polygon& container, SimplePolygon start) {
	Flooder flooder(container, start);
	flooder.run();
	return flooder.cur;
}

SimplePolygon flood(const Instance& inst, SimplePolygon start) {
	return flood(inst.polygon(), start);
}
