#include "cgshop2023_core/verify.hpp"
#include "cpp_instance.hpp"
#include <exception>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <type_traits>

std::string remove_ext(std::string s) {
	size_t last_slash = s.find_last_of("/");
	if (last_slash != std::string::npos)
		s = s.substr(last_slash + 1);

	size_t last_dot = s.find_first_of(".");
	if (last_dot == std::string::npos)
		return s;
	return s.substr(0, last_dot);
}

const std::string OUTPATH = "solutions/";
const std::string INPATH = "instances/";
const std::string INEXT = ".instance.json";
const std::string OUTEXT = ".solution.json";

std::string in_file_full(std::string name) {
	return INPATH + remove_ext(name) + INEXT;
}

std::string out_file_full(std::string name) {
	return OUTPATH + remove_ext(name) + OUTEXT;
}

namespace cgshop2023 {

using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;

template <typename K, typename V>
static void write_kv(std::ostream& output, const K& key, const V& value) {
	output << '\"' << key << "\": ";
	if constexpr (std::is_integral_v<V> || std::is_floating_point_v<V>) {
		output << value;
	} else {
		output << '\"' << value << '\"';
	}
}

template <typename P>
static void write_point(std::ostream& output, const P& p) {
	output << "{\"x\": " << int(std::round(CGAL::to_double(p.x())))
				 << ", \"y\": " << int(std::round(CGAL::to_double(p.y()))) << "}";
}

template <typename V>
static void write_num_exact(std::ostream& output, const V& v) {
	auto gmpq = CGAL::Gmpq(v.exact());
	output << "{\"num\": ";
	output << gmpq.numerator();
	output << ",\"den\": ";
	output << gmpq.denominator();
	output << "}";
}

template <typename P>
static void write_point_exact(std::ostream& output, const P& p) {
	output << "{\"x\": ";
	write_num_exact(output, p.x());
	output << ", \"y\": ";
	write_num_exact(output, p.y());
	output << "}";
}

template <typename PC>
static void write_point_container(std::ostream& out, const PC& c) {
	out << '[';
	auto beg = c.begin();
	auto last = c.end();
	--last;
	for (; beg != last; ++beg) {
		write_point(out, *beg);
		out << ", ";
	}
	write_point(out, *last);
	out << ']';
}

template <typename PC>
static void write_point_container_exact(std::ostream& out, const PC& c) {
	out << '[';
	auto beg = c.begin();
	auto last = c.end();
	--last;
	for (; beg != last; ++beg) {
		write_point_exact(out, *beg);
		out << ", ";
	}
	write_point_exact(out, *last);
	out << ']';
}

template <typename PCC>
static void write_point_container_container_exact(std::ostream& out,
																									const PCC& c) {
	out << "\"polygons\": [\n";
	auto beg = c.begin();
	auto last = c.end();
	--last;
	for (; beg != last; ++beg) {
		write_point_container_exact(out, *beg);
		out << ",\n";
	}
	write_point_container_exact(out, *last);
	out << "\n]";
}

void Instance::write(std::ostream& output, const std::string& name) {
	output << '{';
	write_kv(output, "type", "CGSHOP2023_Instance");
	output << ", ";
	write_kv(output, "name", name);
	output << ", ";
	write_kv(output, "n", num_vertices());
	output << ", \"outer_boundary\": ";
	write_point_container(output, m_polygon.outer_boundary().container());
	output << ", ";
	const auto& hc = m_polygon.holes();
	if (!hc.empty()) {
		output << "\"holes\": [";
		auto hlast = hc.end();
		--hlast;
		for (auto hi = hc.begin(); hi != hlast; ++hi) {
			write_point_container(output, hi->container());
			output << ", ";
		}
		write_point_container(output, hlast->container());
		output << ']';
	} else {
		output << "\"holes\": []";
	}
	output << "}\n";
}

void Solution::write(std::ostream& output, const std::string& name) {
	output << '{';
	write_kv(output, "type", "CGSHOP2023_Solution");
	output << ",\n";
	write_kv(output, "instance", name);
	output << ",\n";
	write_point_container_container_exact(output, polygons());
	output << "}";
}

static Kernel::FT int64_to_cgal_exact(std::int64_t v) {
	double lo32 = v & 0xffffffff;
	double hi32 = double(v >> 32) * 4294967296.0;
	return Kernel::FT(lo32) + Kernel::FT(hi32);
}

static SimplePolygon json_to_points(const nlohmann::json& plist) {
	std::vector<Point> points;
	for (const auto& p : plist) {
		std::int64_t x = p.at("x").get<std::int64_t>();
		std::int64_t y = p.at("y").get<std::int64_t>();
		points.emplace_back(int64_to_cgal_exact(x), int64_to_cgal_exact(y));
	}
	return SimplePolygon(points.begin(), points.end());
}

static Kernel::FT rational_to_cgal_exact(const nlohmann::json& val) {
	Kernel::FT num = int64_to_cgal_exact(val.at("num").get<std::int64_t>());
	Kernel::FT dem = int64_to_cgal_exact(val.at("dem").get<std::int64_t>());
	return num / dem;
}

static SimplePolygon json_to_rational_points(const nlohmann::json& plist) {
	std::vector<Point> points;
	for (const auto& p : plist) {
		const auto& x = p.at("x");
		const auto& y = p.at("y");
		points.emplace_back(rational_to_cgal_exact(x), rational_to_cgal_exact(y));
	}
	return SimplePolygon(points.begin(), points.end());
}

Instance Instance::read(std::istream& input, std::string& out_name) {
	nlohmann::json jsdata;
	input >> jsdata;
	if (jsdata.at("type") != "CGSHOP2023_Instance") {
		throw std::runtime_error("Not a CGSHOP 2023 instance file!");
	}
	out_name = jsdata.at("name").get<std::string>();
	const auto& ob = jsdata.at("outer_boundary");
	const auto& holes = jsdata.at("holes");
	std::vector<SimplePolygon> out_holes;
	SimplePolygon boundary = json_to_points(ob);
	for (const auto& h : holes) {
		out_holes.emplace_back(json_to_points(h));
	}
	return Instance(
			Polygon(std::move(boundary), out_holes.begin(), out_holes.end()));
}

Solution Solution::read(std::istream& input) {
	nlohmann::json jsdata;
	input >> jsdata;
	if (jsdata.at("type") != "CGSHOP2023_Solution") {
		throw std::runtime_error("Not a CGSHOP 2023 solution file!");
	}
	const auto& polygons = jsdata.at("polygons");
	std::vector<SimplePolygon> out_polygons;
	for (const auto& p : polygons) {
		out_polygons.emplace_back(json_to_rational_points(p));
	}
	return Solution(std::move(out_polygons));
}

Instance Instance::read_file(const std::string& name) {
	auto fileloc = in_file_full(name);
	ifstream ifs(fileloc);
	string out_name;
	Instance inst = Instance::read(ifs, out_name);
	return inst;
}

Solution Solution::read_file(const std::string& name) {
	auto fileloc = out_file_full(name);
	ifstream ifs(fileloc);
	string out_name;
	Solution sol = Solution::read(ifs);
	return sol;
}

bool solution_exists(const std::string& name) {
	string filename = out_file_full(name);
	std::ifstream ifs(filename);
	return ifs.is_open();
}

bool Solution::write_if_better(const Instance& inst, const std::string& name) {
	SolutionVerifier svN(&inst, this);
	if (!solution_exists(name) && svN.verify()) {
		cerr << "No existing saved solution for " << name
				 << ". This one is valid, writing it (size=" << size() << ")." << endl;
		string out_name = out_file_full(name);
		ofstream ofs(out_name);
		write(ofs, remove_ext(name));
		return true;
	}
	Solution sol_old = Solution::read_file(name);
	SolutionVerifier svO(&inst, &sol_old);
	if (!svN.verify()) {
		cerr << "Warning: Tried to save invalid solution for " << name << endl;
		return false;
	}
	if (!svO.verify() || size() < sol_old.size()) {
		cerr << "Found solution improvement for " << name << ": " << sol_old.size()
				 << "->" << size() << endl;
		string out_name = out_file_full(name);
		ofstream ofs(out_name);
		write(ofs, remove_ext(name));
		return true;
	}
	return false;
}

} // namespace cgshop2023
