#pragma once

#include "cgshop2023_core/cpp_instance.hpp"
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <iostream>

using namespace cgshop2023;
using namespace std;

struct FaceInfo2 {
	FaceInfo2() {}
	int nesting_level = -2;
	bool in_domain() { return (nesting_level + 2) % 2 == 0; }
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
					if (ct.is_constrained(e))
						border.push_back(e);
					else
						queue.push_back(n);
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
	cdt.insert_constraint(polygon_with_holes.outer_boundary().begin(),
												polygon_with_holes.outer_boundary().end());
	for (auto hole : polygon_with_holes.holes()) {
    hole.reverse_orientation();
		cdt.insert_constraint(hole.begin(), hole.end());
	}

	mark_domains(cdt);
	assert(cdt.is_valid());
	// CGAL::draw(cdt);

	std::vector<SimplePolygon> polys;

	for (auto it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it) {
    cerr << "polygon found, nesting level=" << it->info().nesting_level << endl;
		if (it->info().in_domain()) {
			auto tri = cdt.triangle(it);
			SimplePolygon poly;
			for (int i = 0; i < 3; ++i)
				poly.push_back(tri[i]);
			polys.emplace_back(poly);
		}
	}
	return Solution(std::move(polys));
}
