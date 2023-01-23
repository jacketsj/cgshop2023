#pragma once

#include "cgshop2023_core/cpp_instance.hpp"
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
using namespace cgshop2023;
using std::pair;
using std::vector;

struct Flooder {
	const Polygon& container;
	SimplePolygon cur;
	vector<Point> container_points;
	vector<Line> container_lines;
	Flooder(const Polygon& _container, const SimplePolygon& _start)
			: container(_container), cur(_start) {
		for (const auto& p : container.outer_boundary())
			container_points.push_back(p);
		for (const auto& hole : container.holes())
			for (const auto& p : hole)
				container_points.push_back(p);
		for (size_t i = 0; i < container_points.size(); ++i)
			for (size_t j = i + 1; j < container_points.size(); ++j)
				container_lines.emplace_back(container_points[i], container_points[j]);
	}
	void flood_raycast(const Point& p, const Point& dir);
	void run();
};

SimplePolygon flood(const Polygon& container, SimplePolygon start);
SimplePolygon flood(const Instance& inst, SimplePolygon start);

Solution flood_init(const Instance& inst);
