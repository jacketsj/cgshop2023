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

SimplePolygon greedy_expand(const Instance& inst, SimplePolygon poly,
														const vector<Point>& desired_coverage,
														bool& allCovered);

// get the area uncovered by partial_cover in inst
vector<Polygon> get_missing(const Instance& inst,
														const vector<SimplePolygon>& partial_cover);

// get the area uncovered by partial_cover in inst, triangulated
vector<SimplePolygon>
get_missing_triangles(const Instance& inst,
											const vector<SimplePolygon>& partial_cover);

// get the area uncovered by full_cover[i!=polygon_i] in inst
vector<Polygon> get_missing_removal(const Instance& inst,
																		const vector<SimplePolygon>& full_cover,
																		size_t polygon_i);

// given a full_cover, and a polygon in that full cover, shrink that polygon
// greedily while maintaining the property that the cover is full
std::optional<SimplePolygon>
minimize_to_necessary(const Instance& inst,
											const vector<SimplePolygon>& full_cover,
											size_t polygon_i);

// compute the total area of a set of disjoint polygons
double compute_area(const vector<SimplePolygon>& polys);

// compute a 'base score' for the removal of a polygon from a cover
double removal_score_base(const Instance& inst,
													const vector<SimplePolygon>& current_cover,
													size_t polygon_i);

// TODO complete this thing!!!
SimplePolygon greedy_flood(const Instance& inst, SimplePolygon poly);

struct conflict_optimizer {
	const Instance& inst;
	vector<SimplePolygon> cover;
	vector<SimplePolygon> cover_restore;
	vector<SimplePolygon> uncovered;
	unsigned max_greedy_update = 10000;
	unsigned max_cover_size = 10000;
	conflict_optimizer(const Instance& _inst, const vector<SimplePolygon>& _cover)
			: inst(_inst), cover(_cover), cover_restore(_cover),
				uncovered(get_missing_triangles(inst, cover)) {}
	unsigned greedily_cover();
	void update_uncovered() { uncovered = get_missing_triangles(inst, cover); }
	void remove_from_cover(size_t i) {
		swap(cover[i], cover.back());
		cover.pop_back();
	}
	void remove_random(size_t count);
	void add_random(size_t count);
	void inner_iterate() {
		// TODO remove some things from the cover (using a better heuristic)
		// TODO for the number things to remove, use simulated annealing (maybe
		// combined with spacial-based conflict-optimization for choices)
		remove_random(1);
		for (unsigned expands = 0;
				 expands < max_greedy_update && !uncovered.empty();) {
			auto num_covered = greedily_cover();
			if (num_covered == 0)
				break;
			expands += num_covered;
			// update uncovered each time (up to some iteration limit)
			update_uncovered();
		}
		if (!uncovered.empty()) {
			// use heuristic to choose a triangle to add to the cover
			// TODO make better heuristic (random for now)
			// TODO use simulated annealing to choose how many, spacial-based
			// conflict-optimization for choices
			add_random(1);
			update_uncovered();
		}
	}
	void revert() {
		cover = cover_restore;
		update_uncovered();
	}
	void commit() { cover_restore = cover; }
	void run(unsigned max_attempts, unsigned max_iters) {
		for (unsigned attempt = 0; attempt < max_attempts; ++attempt) {
			for (unsigned iter = 0; iter < max_iters; ++iter) {
				inner_iterate();
				// if cover is too big or any other threshold is hit, revert
				if (cover.size() > max_cover_size) {
					revert();
					break;
				}
				// if an improvement has been found, commit it
				if (uncovered.empty() && cover.size() < cover_restore.size()) {
					commit();
					break;
				}
			}
			revert();
		}
	}
};

void run_co(Instance& inst, Solution& sol, unsigned max_attempts = 1000,
						unsigned max_iters = 10000);

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
