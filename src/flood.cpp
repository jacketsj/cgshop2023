#include "flood.hpp"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/intersections.h>
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

SimplePolygon extend_to_point(const vector<Point>& pointset, const Point& p) {
	vector<Point> points = pointset;
	points.push_back(p);
	vector<Point> chull;
	CGAL::ch_graham_andrew(points.begin(), points.end(),
												 std::back_inserter(chull));
	SimplePolygon newPoly(chull.begin(), chull.end());
	return newPoly;
}

void Flooder::flood_raycast(const Point& p, const Point& dir) {
	// TODO
	// compute arrangement (or maintain)
	// look at all the lines the ray intersects
	// filter for intersection points inside the polygon
	// binary search along the intersection points
	// at each one, check if the expansion to include the point (i.e. convex
	// hull w/ that point) is inside container
	vector<Point> current_points;
	vector<Line> relevant_lines = container_lines;
	for (const Point& p1 : cur) {
		current_points.push_back(p1);
		for (const Point& p2 : container_points)
			relevant_lines.emplace_back(p1, p2);
	}
	for (size_t i = 0; i < current_points.size(); ++i)
		for (size_t j = i + 1; j < current_points.size(); ++j)
			relevant_lines.emplace_back(current_points[i], current_points[j]);

	Ray r(p, dir);
	vector<pair<Kernel::FT, Point>> possible_points;
	for (auto& line : relevant_lines) {
		auto intersection = CGAL::intersection(r, line);
		if (intersection)
			if (const Point* p0 = boost::get<Point>(&intersection.get()))
				possible_points.emplace_back(CGAL::squared_distance(p, *p0), *p0);
	}

	sort(possible_points.begin(), possible_points.end());

	assert(!possible_points.empty());
	// TODO filter out points in possible points that are outside container
	// TODO filter out points in possible points that are strictly inside cur
	// (optional)
	// TODO sort points in possible points by their distance from p
	// TODO binary search on points q in possible points, checking if the convex
	// hull of current_points+q is contained within container each time (not
	// strictly), i.e. if the union of both expands the area
	auto is_outside_comp = [&](int y, const pair<Kernel::FT, Point>& qq) {
		const Point& q = qq.second;
		// compute convex hull of q merged with cur
		SimplePolygon newPoly = extend_to_point(current_points, q);
		bool inside = false;
		// CGAL oriented side and complement method
		vector<Polygon> polygon_c;
		CGAL::complement(container, std::back_inserter(polygon_c));
		// assert(polygon_c.size() == 1);
		inside = true;
		// auto generalNewPoly = GeneralPolygon(newPoly.begin(), newPoly.end());
		for (auto& outside_piece : polygon_c) {
			// does the exterior of one intersect the interior of the other
			inside = inside && CGAL::oriented_side(outside_piece, newPoly) !=
														 CGAL::ON_POSITIVE_SIDE;
		}
		return !inside;
	};
	auto first_incorrect = std::upper_bound(
			possible_points.begin(), possible_points.end(), 0, is_outside_comp);
	size_t distance = first_incorrect - possible_points.begin();
	if (distance > 0) {
		cur = extend_to_point(current_points, possible_points[distance - 1].second);
	}
}

void Flooder::run() {
	Point interior = cur[0];
	// TODO better orders
	for (const auto& dir : container.outer_boundary())
		flood_raycast(interior, dir);
	for (const auto& hole : container.holes())
		for (const auto& dir : hole)
			flood_raycast(interior, dir);
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

SimplePolygon flood(const Instance& inst, const Point& p) {
	Point points[] = {p};
	SimplePolygon start(points, points + 1);

	return flood(inst.polygon(), start);
}

Solution flood_init(const Instance& inst) {
	vector<SimplePolygon> result;
	for (const auto& p : inst.polygon().outer_boundary())
		result.push_back(flood(inst, p));
	for (const auto& hole : inst.polygon().holes())
		for (const auto& p : hole)
			result.push_back(flood(inst, p));
	return Solution(result);
}
