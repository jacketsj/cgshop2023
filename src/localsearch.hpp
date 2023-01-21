#pragma once

#include "cgshop2023_core/cpp_instance.hpp"
#include "cgshop2023_core/verify.hpp"
#include "globals.hpp"
#include <CGAL/ch_graham_andrew.h>

using namespace cgshop2023;
using namespace std;

using Kernel = CGAL::Epeck;
using Point = CGAL::Point_2<Kernel>;
using Polygon = CGAL::Polygon_with_holes_2<Kernel>;
using GeneralPolygon = CGAL::General_polygon_with_holes_2<Kernel>;
using SimplePolygon = CGAL::Polygon_2<Kernel>;

SimplePolygon greedy_expand(Instance& inst, SimplePolygon poly,
														const vector<Point>& desired_coverage,
														bool& allCovered);

// get the area uncovered by partial_cover in inst
vector<SimplePolygon> get_missing(const Instance& inst,
																	const vector<SimplePolygon>& partial_cover);

// get the area uncovered by full_cover[i!=polygon_i] in inst
vector<SimplePolygon>
get_missing_removal(const Instance& inst,
										const vector<SimplePolygon>& full_cover, size_t polygon_i);

// given a full_cover, and a polygon in that full cover, shrink that polygon
// greedily while maintaining the property that the cover is full
SimplePolygon minimize_to_necessary(const Instance& inst,
																		const vector<SimplePolygon>& full_cover,
																		size_t polygon_i);

// compute the total area of a set of disjoint polygons
double compute_area(const vector<SimplePolygon>& polys);

// compute a 'base score' for the removal of a polygon from a cover
double removal_score_base(const Instance& inst,
													const vector<SimplePolygon>& current_cover,
													size_t polygon_i);

// try to greedily remove a polygon i from a solution sol, by covering the area
// with a different polygon options:
// - minimize: reduce the polygon to just what's necessary first, and do the
// same with potential replacers
// - randomize: don't just go through the possible replacements one-by-one
// - replacement_choices: don't go through all of the possible replacements,
// limit to a number
bool try_removal(Instance& inst, Solution& sol, size_t polygon_i,
								 bool randomize, bool minimize, size_t replacement_choices);

// perform try_removal, and commit if it succeeds
void removal_if_possible(Instance& inst, Solution& sol, size_t polygon_i,
												 bool randomize, bool minimize,
												 size_t replacement_choices);

// greedily try to remove all polygons (with some params), up to a limit
// removal_attempts
void try_remove_all(Instance& inst, Solution& sol, bool randomize,
										size_t removal_attempts, bool minimize,
										size_t replacement_choices);
