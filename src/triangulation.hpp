#pragma once

#include "cgshop2023_core/cpp_instance.hpp"
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Vector_2.h>
//#include <CGAL/draw_triangulation_2.h>
#include <iostream>
#include <utility>

using namespace cgshop2023;
using namespace std;

struct FaceInfo2 {
	FaceInfo2() {}
	int nesting_level;
	bool in_domain() { return (nesting_level + 2) % 2 == 1; }
};

typedef CGAL::Epeck K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Face_handle Face_handle;
typedef CDT::Vertex_handle Vertex_handle;
typedef CGAL::Vector_2<K> Vector;

set<pair<Vertex_handle, Vertex_handle>> constrained_edges;

pair<Vertex_handle, Vertex_handle> edge_vertices(CDT::Edge e) {
	return make_pair(e.first->vertex((e.second + 1) % 3),
									 e.first->vertex((e.second + 2) % 3));
}

bool constrained(CDT::Edge e) {
	return constrained_edges.count(edge_vertices(e)) > 0;
}

void mark_domains(CDT& ct, Face_handle start, int index,
									std::list<CDT::Edge>& border) {
	if (start->info().nesting_level != -1) {
		return;
	}
	std::list<Face_handle> queue;
	queue.push_back(start);
	while (!queue.empty()) {
		Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1) {
			fh->info().nesting_level = index;
			for (int i = 0; i < 3; i++) {
				CDT::Edge e(fh, i);
				Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1) {
					auto tri = ct.triangle(fh);
					// if (ct.is_constrained(e)) {
					if (constrained(e)) {
						// cerr << "Border edge found: e=(" << e.first->x() << endl;
						// cerr << "Border edge found (" << tri[(i + 1) % 3].x() << ','
						//		 << tri[(i + 1) % 3].y() << ")-(" << tri[(i + 2) % 3].x() <<
						//','
						//		 << tri[(i + 2) % 3].y() << ')' << endl;
						border.push_back(e);
					} else {
						// cerr << "NON Border edge found (" << tri[(i + 1) % 3].x() << ','
						//		 << tri[(i + 1) % 3].y() << ")-(" << tri[(i + 1) % 3].x() <<
						//','
						//		 << tri[(i + 1) % 3].y() << ')' << endl;
						queue.push_back(n);
					}
				}
			}
		}
	}
}

// explore set of facets connected with non constrained edges,
// and attribute to each such set a nesting level.
// We start from facets incident to the infinite vertex, with a nesting
// level of 0. Then we recursively consider the non-explored facets incident
// to constrained edges bounding the former set and increase the nesting level
// by 1. Facets in the domain are those with an odd nesting level.
void mark_domains(CDT& cdt) {
	for (CDT::Face_handle f : cdt.all_face_handles()) {
		f->info().nesting_level = -1;
	}
	std::list<CDT::Edge> border;
	mark_domains(cdt, cdt.infinite_face(), 0, border);
	while (!border.empty()) {
		CDT::Edge e = border.front();
		border.pop_front();
		Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1) {
			mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
		}
	}
}

Solution basicTriangulation(const Instance& inst) {
	auto polygon_with_holes = inst.polygon();
	CDT cdt;
	vector<vector<Vertex_handle>> boundaries;
	vector<Vertex_handle> outer_boundary;
	for (auto& vert : inst.polygon().outer_boundary()) {
		outer_boundary.emplace_back(cdt.insert(vert));
		// outer_boundary.back()->info() = 1;
	}
	boundaries.push_back(outer_boundary);
	for (auto& hole : polygon_with_holes.holes()) {
		vector<Vertex_handle> hole_boundary;
		for (auto& vert : hole) {
			hole_boundary.emplace_back(cdt.insert(vert));
			// hole_boundary.back()->info() = 2;
		}
		boundaries.push_back(hole_boundary);
	}
	for (auto& boundary : boundaries) {
		for (size_t i = 0; i < boundary.size(); ++i) {
			size_t j = (i + 1) % boundary.size();
			cdt.insert_constraint(boundary[i], boundary[j]);
			constrained_edges.emplace(boundary[i], boundary[j]);
			constrained_edges.emplace(boundary[j], boundary[i]);
		}
	}

	mark_domains(cdt);

	std::vector<SimplePolygon> polys;

	for (auto it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it) {
		auto tri = cdt.triangle(it);
		SimplePolygon poly;
		Vector v(0, 0);
		for (int i = 0; i < 3; ++i) {
			poly.push_back(tri[i]);
			v += tri[i] - Point(0, 0);
		}
		v /= 3;
		Point p = Point(0, 0) + v;
		auto correct_side = it->info().in_domain();
		if (correct_side)
			polys.emplace_back(poly);
	}
	cerr << "Total polys(triangles): " << polys.size() << endl;
	// CGAL::draw(cdt);
	return Solution(std::move(polys));

	/*
	auto insertPolygonConstraints = [&](auto& poly) {
		vector<std::pair<size_t, size_t>> pairs;
		// for (auto it = poly.begin(); it != poly.end(); ++it) {
		for (size_t i = 0; i < poly.size(); ++i) {
			pairs.emplace_back(i, (i + 1) % poly.size());
			pairs.emplace_back((i + 1) % poly.size(), i);
			cerr << "Inserting constraint: (" << poly[i].x() << "," << poly[i].y()
					 << ")-(" << poly[(i + 1) % poly.size()].x() << ","
					 << poly[(i + 1) % poly.size()].y() << ")\n";
			// auto itNext = it;
			// itNext++;
			// if (itNext == poly.end())
			//	itNext = poly.begin();
			// cdt.insert_constraint(it, itNext);
		}
		cdt.insert_constraints(poly.begin(), poly.end(), pairs.begin(),
													 pairs.end());
	};
	// insertPolygonConstraints(polygon_with_holes.outer_boundary());
	cdt.insert_constraint(polygon_with_holes.outer_boundary().begin(),
												polygon_with_holes.outer_boundary().end());
	// cdt.insert_constraint(polygon_with_holes.outer_boundary().back(),
	// polygon_with_holes.outer_boundary().front());
	for (auto hole : polygon_with_holes.holes()) {
		// hole.reverse_orientation();
		// insertPolygonConstraints(hole);
		cdt.insert_constraint(hole.begin(), hole.end());
	}

	//for (auto& e : cdt.constrained_edges()) {
	//	auto tri = cdt.triangle(e.first);
	//	auto a = tri[e.second];
	//	auto b = tri[(e.second + 1) % 3];
	//	//printf("constrained edge: (%f,%f)-(%f,%f)\n", a.x(), a.y(), b.x(),
	b.y());
	//}

	mark_domains(cdt);
	assert(cdt.is_valid());
	// CGAL::draw(cdt);

	std::vector<SimplePolygon> polys;

	for (auto it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it) {
		// cerr << "polygon found, nesting level=" << it->info().nesting_level <<
		// endl; if (it->info().in_domain()) {
		auto tri = cdt.triangle(it);
		SimplePolygon poly;
		Vector v;
		bool failed = false;
		for (int i = 0; i < 3; ++i) {
			poly.push_back(tri[i]);
			v += tri[i] - Point(0, 0); // Vector(tri[i].x(), tri[i].y());
			failed = failed || CGAL::oriented_side(tri[i], polygon_with_holes) ==
														 CGAL::NEGATIVE;
		}
		if (failed) {
			cerr << "vertex of triangle found to be outside polygon\n";
			continue;
		}
		v /= 3;
		// Point p(v.x(), v.y());
		Point p = Point(0, 0) + v;
		if (polygon_with_holes.holes().empty() ||
				CGAL::oriented_side(p, polygon_with_holes) != CGAL::NEGATIVE)
			polys.emplace_back(poly);
		//}
	}
	return Solution(std::move(polys));
	*/
}
