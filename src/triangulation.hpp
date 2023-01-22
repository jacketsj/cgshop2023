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

struct Triangulator {
	set<pair<Vertex_handle, Vertex_handle>> constrained_edges;

	pair<Vertex_handle, Vertex_handle> edge_vertices(CDT::Edge e);

	bool constrained(CDT::Edge e);

	void mark_domains(CDT& ct, Face_handle start, int index,
										std::list<CDT::Edge>& border);

	// explore set of facets connected with non constrained edges,
	// and attribute to each such set a nesting level.
	// We start from facets incident to the infinite vertex, with a nesting
	// level of 0. Then we recursively consider the non-explored facets incident
	// to constrained edges bounding the former set and increase the nesting level
	// by 1. Facets in the domain are those with an odd nesting level.
	void mark_domains(CDT& cdt);
};

std::vector<SimplePolygon> getTriangles(const Polygon& polygon_with_holes);

std::vector<SimplePolygon>
getTriangles(const std::vector<Polygon>& polygons_with_holes);

Solution basicTriangulation(const Instance& inst);
