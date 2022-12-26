#pragma once

#include "cgshop2023_core/cpp_instance.hpp"
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <init_ogl_context.h>

using namespace cgshop2023;

using Kernel = CGAL::Epeck;
typedef CGAL::Polygon_set_2<Kernel> Polygon_set;

/*
namespace CGAL {

// Viewer class for Polygon_with_holes_2
template <class P2>
class SimplePolygonWithHoles2ViewerQt : public Basic_viewer_qt {
	typedef Basic_viewer_qt Base;
	typedef typename P2::General_polygon_2::Point_2 Point;

public:
	/// Construct the viewer without drawing anything.
	/// @param title the title of the window
	SimplePolygonWithHoles2ViewerQt(
			QWidget* parent, const char* title = "Basic Polygon_with_holes_2 Viewer")
			: Base(parent, title, true, true, true, false, false) {
		clear();
	}

	/// Construct the viewer.
	/// @param ap2 the polygon to view
	/// @param title the title of the window
	SimplePolygonWithHoles2ViewerQt(
			QWidget* parent, const P2& ap2,
			const char* title = "Basic Polygon_with_holes_2 Viewer")
			: // First draw: vertices; edges, faces; multi-color; no inverse normal
				Base(parent, title, true, true, true, false, false) {
		clear();
		compute_elements(ap2);
	}

protected:
	void compute_one_loop_elements(const typename P2::General_polygon_2& p,
																 bool hole) {
		if (hole) {
			add_point_in_face(p.vertex(p.size() - 1));
		}

		typename P2::General_polygon_2::Vertex_const_iterator prev;
		for (typename P2::General_polygon_2::Vertex_const_iterator i =
						 p.vertices_begin();
				 i != p.vertices_end(); ++i) {
			add_point(*i); // Add vertex
			if (i != p.vertices_begin()) {
				add_segment(*prev, *i);
			}											 // Add segment with previous point
			add_point_in_face(*i); // Add point in face
			prev = i;
		}

		// Add the last segment between the last point and the first one
		add_segment(*prev, *(p.vertices_begin()));
	}

	void compute_elements(const P2& p2) {
		if (p2.outer_boundary().is_empty())
			return;

		CGAL::Color c(75, 160, 255);
		face_begin(c);

		compute_one_loop_elements(p2.outer_boundary(), false);

		for (typename P2::Hole_const_iterator it = p2.holes_begin();
				 it != p2.holes_end(); ++it) {
			compute_one_loop_elements(*it, true);
			add_point_in_face(
					p2.outer_boundary().vertex(p2.outer_boundary().size() - 1));
		}

		face_end();
	}

	virtual void keyPressEvent(QKeyEvent* e) {
		// Test key pressed:
		//    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
		//    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }

		// Call: * compute_elements() if the model changed, followed by
		//       * redraw() if some viewing parameters changed that implies some
		//                  modifications of the buffers
		//                  (eg. type of normal, color/mono)
		//       * update() just to update the drawing

		// Call the base method to process others/classicals key
		Base::keyPressEvent(e);
	}
};

// Specialization of draw function.
template <class T, class C>
void draw(const CGAL::Polygon_with_holes_2<T, C>& ap2,
					const char* title = "Polygon_with_holes_2 Basic Viewer") {
	CGAL::Qt::init_ogl_context(4, 3);
	int argc = 1;
	const char* argv[2] = {"t2_viewer", nullptr};
	QApplication app(argc, const_cast<char**>(argv));
	SimplePolygonWithHoles2ViewerQt<CGAL::Polygon_with_holes_2<T, C>> mainwindow(
			app.activeWindow(), ap2, title);
	mainwindow.show();
	app.exec();
}

template <class PS2>
class SimplePolygonSet2ViewerQt : public SimplePolygonWithHoles2ViewerQt<
																			typename PS2::Polygon_with_holes_2> {
	typedef SimplePolygonWithHoles2ViewerQt<typename PS2::Polygon_with_holes_2>
			Base;

public:
	SimplePolygonSet2ViewerQt(QWidget* parent, const PS2& aps2,
														const char* title = "Basic Polygon_set_2 Viewer")
			: Base(parent, title) {
		std::vector<typename PS2::Polygon_with_holes_2> polygons;
		aps2.polygons_with_holes(std::back_inserter(polygons));

		for (typename PS2::Polygon_with_holes_2& P : polygons) {
			Base::compute_elements(P);
		}
	}
};

// Specialization of draw function.
template <class T, class C>
void draw(const CGAL::Polygon_set_2<T, C>& aps2,
					const char* title = "Polygon_set_2 Basic Viewer") {
	CGAL::Qt::init_ogl_context(4, 3);
	int argc = 1;
	const char* argv[2] = {"t2_viewer", nullptr};
	QApplication app(argc, const_cast<char**>(argv));
	SimplePolygonSet2ViewerQt<CGAL::Polygon_set_2<T, C>> mainwindow(
			app.activeWindow(), aps2, title);
	mainwindow.show();
	app.exec();
}

} // End namespace CGAL
*/

using Polygon = CGAL::Polygon_with_holes_2<Kernel>;
using SimplePolygon = CGAL::Polygon_2<Kernel>;

void draw(Instance& inst, Solution& sol) {
	// TODO use polygon set draw, need all the includes for stuff

	Polygon combined(inst.polygon().outer_boundary(), sol.polygons().begin(),
									 sol.polygons().end());
	CGAL::draw(combined);

	// Polygon_set S;
	// Polygon for (auto& poly : s.polygons()) { S.insert(poly); }

	// CGAL::draw(S);
}
