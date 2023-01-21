#include "localsearch.hpp"

#include "cgshop2023_core/cpp_instance.hpp"
#include "cgshop2023_core/verify.hpp"
#include "globals.hpp"
#include <CGAL/ch_graham_andrew.h>
#include <random>

using namespace cgshop2023;
using namespace std;

using Kernel = CGAL::Epeck;
using Point = CGAL::Point_2<Kernel>;
using Polygon = CGAL::Polygon_with_holes_2<Kernel>;
using GeneralPolygon = CGAL::General_polygon_with_holes_2<Kernel>;
using SimplePolygon = CGAL::Polygon_2<Kernel>;

std::random_device rd;
std::mt19937 g(rd());

SimplePolygon greedy_expand(Instance& inst, SimplePolygon poly,
														const vector<Point>& desired_coverage,
														bool& allCovered) {
	allCovered = true;
	auto outerArea = area(inst.polygon());
	// decide wether to add p to poly
	// look at what the new polygon would be by computing convex hull (very lazy
	// method) this could be improved to O(log n) with binary search (and the
	// new edges could be singled out)
	vector<Point> points(poly.vertices_begin(), poly.vertices_end());
	for (const Point& p : desired_coverage)
		points.push_back(p);
	// get the convex hull
	vector<Point> chull;
	CGAL::ch_graham_andrew(points.begin(), points.end(),
												 std::back_inserter(chull));
	SimplePolygon newPoly(chull.begin(), chull.end());

	// is it inside the polygon
	bool inside = false;
	// optimal method: Triangulate the base polygon, create point location data
	// structure and compute tree-cotree decomposition. Check if new loop is
	// trivial, and if any edge intersects the outside (can check only the two
	// new edges for this, and precompute the z2-homology values along paths for
	// all) complexity becomes the complexity of checking whether a line segment
	// intersects a polygon (does a data structure exist for this? still fairly
	// efficient in practice with point location) that's a lot of work to
	// implement, so instead we use cgal's methods (probably much, much slower)

	// CGAL oriented side and complement method
	vector<Polygon> polygon_c;
	CGAL::complement(inst.polygon(), std::back_inserter(polygon_c));
	// assert(polygon_c.size() == 1);
	inside = true;
	// auto generalNewPoly = GeneralPolygon(newPoly.begin(), newPoly.end());
	for (auto& outside_piece : polygon_c) {
		// does the exterior of one intersect the interior of the other
		inside = inside && CGAL::oriented_side(outside_piece, newPoly) !=
													 CGAL::ON_POSITIVE_SIDE;
	}

	// CGAL area computation method
	/*
	std::vector<Polygon> coverage = {};
	std::vector<Polygon> polys = {inst.polygon(), Polygon(newPoly)};
	CGAL::join(polys.begin(), polys.end(), std::back_inserter(coverage));
	inside = inside || (coverage.size() == 1 && area(coverage[0]) <= outerArea);
	*/

	if (inside) {
		// it is inside, so add it
		poly = newPoly;
	} else
		allCovered = false;

	return poly;
}

vector<SimplePolygon> get_missing(const Instance& inst,
																	const vector<SimplePolygon>& partial_cover) {
	vector<Polygon> to_join;
	for (auto& poly : partial_cover) {
		to_join.push_back(Polygon(poly));
	}
	vector<Polygon> polygon_c;
	CGAL::complement(inst.polygon(), std::back_inserter(polygon_c));
	for (const auto& poly : polygon_c)
		to_join.push_back(poly);
	std::vector<Polygon> coverage = {};
	CGAL::join(to_join.begin(), to_join.end(), std::back_inserter(coverage));

	// TODO what if this is empty instead? need extra routine for that probably
	assert(coverage.size() == 1);

	vector<SimplePolygon> output;
	for (auto& poly : coverage[0].holes())
		output.push_back(poly);

	return output;
}

vector<SimplePolygon>
get_missing_removal(const Instance& inst,
										const vector<SimplePolygon>& full_cover, size_t polygon_i) {
	vector<SimplePolygon> partial_cover;
	for (size_t i = 0; i < full_cover.size(); ++i) {
		if (i == polygon_i)
			continue;
		partial_cover.push_back(full_cover[i]);
	}
	return get_missing(inst, partial_cover);
}

SimplePolygon minimize_to_necessary(const Instance& inst,
																		const vector<SimplePolygon>& full_cover,
																		size_t polygon_i) {
	auto missing = get_missing_removal(inst, full_cover, polygon_i);
	vector<Point> pointset;
	for (const auto& poly : missing) {
		for (const auto& pt : poly)
			pointset.push_back(pt);
	}
	// get the convex hull
	vector<Point> chull;
	CGAL::ch_graham_andrew(pointset.begin(), pointset.end(),
												 std::back_inserter(chull));
	SimplePolygon newPoly(chull.begin(), chull.end());
	return newPoly;
}

double compute_area(const vector<SimplePolygon>& polys) {
	double cur_area(0);
	for (auto& poly : polys) {
		cur_area += CGAL::to_double(poly.area());
	}
	return double(cur_area);
}

double removal_score_base(const Instance& inst,
													const vector<SimplePolygon>& current_cover,
													size_t polygon_i) {
	auto missing_current = get_missing(inst, current_cover);
	auto missing_next = get_missing_removal(inst, current_cover, polygon_i);
	auto total_area_current = compute_area(missing_current);
	auto total_area_next = compute_area(missing_next);
	auto count_current = missing_current.size();
	auto count_next = missing_next.size();

	// negative = worse score
	// i have no idea what I'm doing here
	auto area_delta = total_area_current - total_area_next;
	auto count_delta = count_next - count_current;
	return 8 * area_delta + 4 * count_delta;
}

bool try_removal(Instance& inst, Solution& sol, size_t polygon_i,
								 bool randomize, bool minimize, size_t replacement_choices) {
	vector<Point> desired_coverage(sol.polygons()[polygon_i].begin(),
																 sol.polygons()[polygon_i].end());
	vector<SimplePolygon> modified_cover;
	if (minimize) {
		SimplePolygon minimal_desired_coverage_full =
				minimize_to_necessary(inst, sol.polygons(), polygon_i);
		desired_coverage = vector<Point>(minimal_desired_coverage_full.begin(),
																		 minimal_desired_coverage_full.end());

		vector<Point> chull;
		CGAL::ch_graham_andrew(desired_coverage.begin(), desired_coverage.end(),
													 std::back_inserter(chull));
		SimplePolygon newPoly(chull.begin(), chull.end());
		modified_cover = sol.polygons();
		modified_cover[polygon_i] = newPoly;
	} else {
		modified_cover = sol.polygons();
	}

	// try to remove polygon_i
	// do it by trying to get polygons in sol to cover it's vertices
	vector<int> to_try;
	for (size_t i = 0; i < sol.polygons().size(); ++i)
		if (i != polygon_i)
			to_try.push_back(i);
	if (randomize)
		shuffle(to_try.begin(), to_try.end(), g);
	bool succeeded = false;
	for (size_t i = 0; i < to_try.size() &&
										 (replacement_choices == 0 || i < replacement_choices);
			 ++i) {
		if (VERBOSE()) {
			// cerr << "Doing " << i << "th replacement choice for polygon " <<
			// polygon_i
			//		 << endl;
		}
		size_t cur_i = to_try[i];
		bool allCovered = false;
		SimplePolygon to_expand = sol.polygons()[cur_i];
		if (minimize) {
			to_expand = minimize_to_necessary(inst, modified_cover, cur_i);
		}
		SimplePolygon newPoly =
				greedy_expand(inst, to_expand, desired_coverage, allCovered);
		if (allCovered) {
			sol.polygons_m()[cur_i] = newPoly;
			succeeded = true;
			break;
		}
	}
	return succeeded;
}

void removal_if_possible(Instance& inst, Solution& sol, size_t polygon_i,
												 bool randomize, bool minimize,
												 size_t replacement_choices) {
	if (try_removal(inst, sol, polygon_i, randomize, minimize,
									replacement_choices)) {
		swap(sol.polygons_m()[polygon_i],
				 sol.polygons_m()[sol.polygons().size() - 1]);
		sol.polygons_m().pop_back();
	}
}

void try_remove_all(Instance& inst, Solution& sol, bool randomize,
										size_t removal_attempts, bool minimize,
										size_t replacement_choices) {
	cerr << "Running try_remove_all on " << sol.polygons().size()
			 << " polygons\n";
	vector<int> to_remove;
	for (int polygon_i = sol.polygons().size() - 1; polygon_i >= 0; --polygon_i)
		to_remove.push_back(polygon_i);
	if (randomize)
		shuffle(to_remove.begin(), to_remove.end(), g);
	for (size_t i = 0;
			 i < to_remove.size() && (removal_attempts == 0 || i < removal_attempts);
			 ++i) {
		if (VERBOSE()) {
			cerr << "Doing " << i << "th try remove (polygon number " << i << ": "
					 << to_remove[i] << ")" << endl;
		}
		removal_if_possible(inst, sol, to_remove[i], randomize, minimize,
												replacement_choices);
	}
	cerr << "Finished running try_remove_all, now have " << sol.polygons().size()
			 << " polygons\n";
}
