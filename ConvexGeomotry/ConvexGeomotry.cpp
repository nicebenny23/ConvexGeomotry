// ConvexGeomotry.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <vector>
#include <iostream>
#pragma once
#include <format>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <optional>
#include <algorithm> 
#include <flat_set>
#include <ranges>
#include "random.h"
#include <format>
#include <boost/multiprecision/cpp_dec_float.hpp>


using precise = boost::multiprecision::cpp_dec_float_50;


template <>
struct std::formatter<precise, char>
	: std::formatter<double, char>
{
	auto format(const precise& value, std::format_context& ctx) const
	{
		return std::formatter<double, char>::format(
			value.convert_to<double>(), ctx
		);
	}
};

#include <unordered_set>
template<typename ...Args>
void print(std::format_string<Args...> fmt,Args&&... args) {
	std::cout << std::format(fmt, std::forward<Args>(args)...)<<'\n';
}
struct Vector {

	std::vector<precise> coords;

	Vector() = default;

	explicit Vector(const std::vector<precise>& v)
		: coords(v) {
	}
	explicit Vector(std::vector<precise>&& v)
		: coords(std::move(v)) {
	}

	explicit Vector(size_t n)
		: coords(n, 0.0) {
	}

	size_t dim() const {
		return coords.size();
	}

	static void check_dimensions(const Vector& a, const Vector& b) {
		if (a.dim() != b.dim())
			throw std::invalid_argument("Point dimensions do not match.");
	}

	static precise dot_product(const Vector& a, const Vector& b) {
		check_dimensions(a, b);

		precise sum = 0.0;
		for (size_t i = 0; i < a.dim(); i++)
			sum += a[i] * b[i];

		return sum;
	}
	//one when perfectly aligned
	static precise aligment0to1(const Vector& a, const Vector& b) {
		return (dot_product(a, b) + 1) / 2.0;

	}static precise distance(const Vector& a, const Vector& b) {
		return (a - b).magnitude();

	}
	precise magnitude() const {
		return boost::multiprecision::sqrt(dot_product(*this, *this));
	}
	Vector normal() const {
		if (dim()==0)
		{
			return *this;
		}
		return *this / magnitude();
	}
	Vector& normalize() {
		*this /= magnitude();
		return *this;
	}
	precise& operator[](size_t i) {
		return coords[i];
	}

	const precise& operator[](size_t i) const {
		return coords[i];
	}

	Vector operator+(const Vector& other) const {
		check_dimensions(*this, other);

		Vector result(dim());
		for (size_t i = 0; i < dim(); i++)
			result[i] = coords[i] + other[i];

		return result;
	}

	void abs_each() {
		for (size_t i = 0; i < dim(); i++)
		{
			coords[i] = boost::multiprecision::abs(coords[i]);
		}
	}
	Vector& operator+=(const Vector& other) {
		for (size_t i = 0; i < dim(); i++)
			coords[i] += other[i];

		return *this;
	}

	Vector operator-() const {
		return *this * -1.0;
	}

	Vector operator-(const Vector& other) const {
		Vector self = *this;
		self -= other;
		return self;
	}

	Vector& operator-=(const Vector& other) {
		for (size_t i = 0; i < dim(); i++) {
			coords[i] -= other.coords[i];
		}
		return *this;
	}


	Vector operator*(precise scalar) const {
		Vector result(dim());

		for (size_t i = 0; i < dim(); i++)
			result[i] = coords[i] * scalar;

		return result;
	}

	Vector& operator*=(precise scalar) {
		for (size_t i = 0; i < dim(); i++)
			coords[i] *= scalar;

		return *this;
	}

	Vector operator/(precise scalar) const {
		if (scalar == 0.0)
			throw std::invalid_argument("Division by zero.");

		return *this * (1.0 / scalar);
	}

	Vector& operator/=(precise scalar) {
		return *this *= 1 / scalar;
	}

	bool operator==(const Vector& other) const {
		return coords == other.coords;
	}

	bool operator!=(const Vector& other) const {
		return !(*this == other);
	}
	static Vector random(size_t dim) {
		Vector v;
		v.coords.reserve(dim);
		for (size_t i = 0; i < dim; i++)
		{
			v.coords.push_back(gaussian());
		}
		return v.normal();

	}
	Vector project_onto_this_normalized_already(const Vector& projected_onto) const {
		return *this * Vector::dot_product(projected_onto, *this);
	}
	Vector project_onto_this(const Vector& projected_onto) const {
		return normal().project_onto_this_normalized_already(projected_onto);
	}

};

//<=n vectors and returns orthagonal ones with the same span
inline std::optional<std::pair<std::vector<Vector>, precise>> gramm_shmitt_process(std::vector<Vector> vectors) {
	precise scale = 1.0;

	for (size_t i = 0; i < vectors.size(); i++) {
		Vector& v = vectors[i];
		

		
		for (int pass = 0; pass < 2; ++pass) {
			for (size_t j = 0; j < i; ++j) {
				v -= vectors[j] * Vector::dot_product(v, vectors[j]);
			}
		}
		precise mag = v.magnitude();
		scale *= mag;

		if (mag >= 1e-8) {
			v /= mag;
		}
		else
		{
			return std::nullopt;
		}
	}
	
	return std::optional(std::pair<std::vector<Vector>, precise>({ vectors, scale }));
}

struct Span {
	precise determinant() const {
		return gramm_shmitt_process(points).transform([](auto&& value) {return value.second;}).value_or(0);
	}
	std::vector<Vector> points;
	size_t dim() const {
		return points.size();
	}

	Vector to_basis(const Vector& point) const {
		Vector pnt(dim());
		for (size_t i = 0; i < dim(); i++)
		{
			precise mag = points[i].magnitude();
			pnt[i] = Vector::dot_product(points[i], (point))/(mag*mag);
		}
		return pnt;
	}
};
struct Ray {
	Vector start;
	Vector direction;
	size_t dim() const {
		return start.dim();
	}
	Ray(const Vector& pnt, const Vector& direction) :start(pnt), direction(direction.normal()) {

	}
	Vector normal() const {
		return direction.normal();
	}
	Vector along(precise dist) const {
		return start + direction*dist;
	}
};
struct HalfSpace {
	Vector normal;
	precise bound;
	HalfSpace(const Vector& norm, precise bnd) :normal(norm), bound(bnd) {
		bound /= norm.magnitude();
		normal.normalize();
	}
	HalfSpace translated(const Vector& p) const {
		return HalfSpace(normal, bound + Vector::dot_product(p, normal));
	}
	void translate(const Vector& p) {
		*this = translated(p);
	}
	HalfSpace operator-() {
		return HalfSpace(-normal, -bound);
	}
	bool operator==(const HalfSpace& other) const = default;

	bool operator!=(const HalfSpace& other) const = default;
	bool contains(const Vector& other) const {
		return contains_exactly(other)||suffieciently_close(other);
	}
	enum class ContainmentState {
		boundry,
		strongly_in,
		strongly_out,
	};
	ContainmentState containment(const Vector& other) const {

		precise state = bound - Vector::dot_product(normal, other);
		precise abs_dist = abs(state);
		if (abs_dist > 1e-8 && abs_dist < 1e-5) {
	//		print("[DEBUG Containment] Point near boundary gray-zone!");
			//print("  Signed distance: {:.10e} | Target threshold: 1e-8", abs_dist);
		}
		
		if (abs(state)<1e-8)
		{
			return ContainmentState::boundry;
		}
		if (state>=0)
		{
			return ContainmentState::strongly_in;
		}
		return ContainmentState::strongly_out;
	}
	bool contains_exactly(const Vector& other) const {
		return Vector::dot_product(normal, other) <= bound ;
	}
	bool suffieciently_close(const Vector& other) const {
		return boost::multiprecision::abs(Vector::dot_product(normal, other) - bound) <= 1e-8;
	}
	Vector project_onto_plane(const Vector& other) const {
		return other + normal * (bound - Vector::dot_product(normal, other));
	}
	Vector project_onto_halfspace(const Vector& other) const {
		return other + normal * min(precise(0.0),(bound - Vector::dot_product(normal, other)));
	}
	size_t dim() const {
		return normal.dim();
	}
	precise dist(const Vector& point) const {
		return (project_onto_plane(point) - point).magnitude();
	}
	//negitive if in
	precise signed_dist(const Vector& pnt) const {
		return (Vector::dot_product(normal, pnt)- bound);
	}
	std::optional<precise> hit_time(const Ray& ray) const {
		precise ht = dist(ray.start) / Vector::dot_product(normal, ray.direction);
		if (ht > 0)
		{
			return ht;
		}
		return std::nullopt;
	}
	Span span() const {
		if (dim() == 0)
		{
			return Span();
		}
		std::vector<Vector> points;
		points.push_back(normal);

		Vector center = project_onto_plane(Vector(dim()));
		for (size_t i = 1; i < dim(); i++)
		{
			Vector unit = center;
			unit[i] += 1;
			points.push_back(project_onto_plane(unit) - center);
		}
		points = gramm_shmitt_process(points).transform([](auto&& val) {return val.first;}).value_or(std::vector<Vector>());
		std::swap(points[0], points.back());
		points.pop_back();
		return Span{ points };
	}
};
namespace std {

	template<>
	struct std::formatter<Vector> {
		constexpr auto parse(std::format_parse_context& ctx) {
			return ctx.begin();
		}

		auto format(const Vector& v, std::format_context& ctx) const {
			auto out = ctx.out();
			*out++ = '(';

			for (size_t i = 0; i < v.dim(); i++)
			{
				if (i != 0)
				{
					*out++ = ',';
					*out++ = ' ';
				}

				out = std::format_to(out, "{:.10}", v[i]);
			}

			*out++ = ')';
			
			return out;
		}
	};
	template<>
	struct std::formatter<HalfSpace> {
		constexpr auto parse(std::format_parse_context& ctx) {
			return ctx.begin();
		}

		auto format(const HalfSpace& v, std::format_context& ctx) const {
			return std::format_to(ctx.out(), "(normal:{})({:.10})", v.normal, v.bound);
		}
	};
	template<typename T>
	struct std::formatter<std::vector<T>>
	{
		constexpr auto parse(std::format_parse_context& ctx)
		{
			return ctx.begin();
		}

		auto format(const std::vector<T>& vec, std::format_context& ctx) const
		{
			auto out = ctx.out();

			out = std::format_to(out, "len({})", vec.size());
			*out++ = '[';

			for (size_t i = 0; i < vec.size(); i++)
			{
				if (i != 0)
				{
					*out++ = ',';
					*out++ = ' ';
				}

				out = std::format_to(out, "{}", vec[i]);
			}

			*out++ = ']';

			return out;
		}
	};
	template<>
	struct hash<Vector> {
		size_t operator()(const Vector& v) const noexcept {
			size_t seed = 0;

			for (precise x : v.coords) {
				size_t h = std::hash<precise>{}(x);

				// hash combine (boost-style)
				seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}

			return seed;
		}
	};
	template<>
	struct hash<HalfSpace> {
		size_t operator()(const HalfSpace& v) const noexcept {
			size_t seed = 0;
			seed=std::hash<Vector>()(v.normal);
				seed ^= std::hash<precise>()(v.bound) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			return seed;
		}
	};
	template<typename T>
	struct hash<std::vector<T>> {
		size_t operator()(const std::vector<T>& vec) const noexcept {
			size_t seed = 0;

			for (const T& v : vec) {
				seed ^= std::hash<T>()(v)
					+ 0x9e3779b9
					+ (seed << 6)
					+ (seed >> 2);
			}

			return seed;
		}
	};
}
	struct GaussianMatrix {
		GaussianMatrix(std::vector<HalfSpace> halfspaces) {
			for (size_t i = 0; i < halfspaces.size(); i++) { 
				Vector next_normal = halfspaces[i].normal;

				next_normal.coords.push_back(halfspaces[i].bound);
				rows.push_back(std::move(next_normal)); 
			}
		}
		std::vector<Vector> rows;
		std::optional<Vector> solve() {
			std::vector<Vector> original_rows = rows;
			for (size_t i = 0; i < rows.size(); i++)
			{
				size_t pivot = i;
				precise best = boost::multiprecision::abs(rows[i][i]);
				//partial pivot
				for (size_t r = i; r < rows.size(); ++r)
				{
					precise val = boost::multiprecision::abs(rows[r][i]);
					if (val >=best)
					{
						best = val;
						pivot = r;
					}
				}
				std::swap(rows[pivot], rows[i]);
				if (abs(rows[i][i])<=1e-8)
				{
					return std::nullopt;
				}
				rows[i] /= rows[i][i];
				Vector row = rows[i];

				for (size_t j = i + 1; j < rows.size(); j++)
				{
					rows[j] -= row * rows[j][i];
				}
			}
			size_t dim = rows.size();
			Vector res(dim);
			for (int i = static_cast<int>(dim) - 1; i >= 0; i--) {
				precise sum = rows[i][dim];
				for (size_t j = i + 1; j < dim; j++) {
					sum -= rows[i][j] * res[j];
				}
				res[i] = sum;
			}
			return res;
		}

	};
	//n vectors
	inline std::optional<HalfSpace> space_from_vectors(const std::vector<Vector>& vectors) {
		if (vectors.empty() || vectors.size() != vectors[0].dim()) {
			throw std::logic_error("cannot construct a half space with this amount of vectors");
		}

		std::vector<Vector> edges;
		edges.reserve(vectors.size() - 1);
		for (size_t i = 1; i < vectors.size(); i++) {
			edges.push_back(vectors[i] - vectors[0]);
		}

		auto basis = gramm_shmitt_process(edges);
		if (!basis) {
			return std::nullopt;
		}

		auto& orthogonal_basis = basis.value().first;

		for (size_t i = 0; i < vectors[0].dim(); i++) {
			Vector normal(vectors[0].dim());
			normal[i] = 1.0;

			for (int pass = 0; pass < 2; pass++) {
				for (const auto& b : orthogonal_basis) {
					normal -= b.project_onto_this_normalized_already(normal);
				}
			}

			if (normal.magnitude() >= 1e-2) {
				normal /= normal.magnitude(); // Ensure normalized normal
				HalfSpace hs(normal, Vector::dot_product(vectors[0], normal));

				// --- AFFINE COMBINATION TEST ---
				// Construct a random point inside the plane span: P = sum(alpha_i * v_i) with sum(alpha_i) = 1
				Vector affine_point(vectors[0].dim());
				precise idk = 0;
				for (size_t k = 0; k < vectors.size(); ++k) {
					precise w = (k!=vectors.size()-1)? precise(k*-9.3): (1-idk);
					idk +=w;
					affine_point += vectors[k] * w;
				}

				precise dist = hs.dist(affine_point);
				if (dist > 1e-10) {
					print( "[DEBUG PLANE] Affine combination test failed! Point dist to plane{}",dist );
					
				}
				return hs;
			}
		}
		throw std::logic_error("not possible");
	}

	struct Vrep {
		void add_unchecked(const Vector& pnt)
		{
			points.push_back(pnt);
		}
		precise support(const Vector& direction) const {
			
			precise spt = -std::numeric_limits<precise>().infinity();
			for (const Vector& v:*this)
			{
				spt = std::max(spt,Vector::dot_product(direction, v));
			}
			return spt;
		}
		bool full_span(size_t trials = 1000) const {
			size_t d = dim();
			if (points.size() < d + 1)
			{
				return false;
			}
			for (size_t t = 0; t < trials; t++)
			{
				std::vector<Vector> edges;
				edges.reserve(d);
				std::vector<Vector> chosen;
				while (chosen.size() < d + 1)
				{
					Vector base = points[random() * points.size()];
					if (!std::ranges::contains(chosen, base))
					{
						chosen.push_back(base);

					}
				}
				for (size_t i = 0; i < d; i++)
				{
					edges.push_back(chosen[i] - chosen.back());
				}
				if (gramm_shmitt_process(edges).has_value())
				{
					return true;
				}
			}

			return false;
		}
		void translate(const Vector& p) {
			for (Vector& v : points)
			{
				v += p;
			}
		}
		using iterator = std::vector<Vector>::iterator;
		iterator begin() {
			return points.begin();
		}
		iterator end() {
			return points.end();
		}
		using const_iterator = std::vector<Vector>::const_iterator;
		const_iterator begin() const {
			return points.begin();
		}
		const_iterator end() const {
			return points.end();
		}
		size_t dim() const {
			return points.empty() ? 0 : points[0].dim();
		}
		bool redudant(const Vector& point) {
			for (const Vector& pnt : points)
			{
				if (Vector::distance(point, pnt) <=1e-8)
				{
					return true;
				}
			}
			return false;
		}
		void add(const Vector& point) {
			for (const Vector& pnt : points)
			{
				if (Vector::distance(point, pnt) < 1e-8)
				{
					return;
				}

				if (Vector::distance(point, pnt) <= 1e-4) {
					precise u = Vector::distance(point, pnt);
					int l = 4;
				}
			}
			points.push_back(point);
		}

		std::vector<Vector> points;
		std::vector<Vector> points_on(const HalfSpace& space) const {
			std::vector<Vector> point_list;
			for (const Vector& v : points)
			{
				HalfSpace::ContainmentState state = space.containment(v);

				if (space.suffieciently_close(v))
				{

					point_list.push_back(v);
				}
			}
			return point_list;
		}
		bool contains_all(HalfSpace space) const {
			for (const Vector& v : points)
			{
				if (!space.contains(v))
				{
					return false;
				}
			}
			return true;
		}
		bool is_supporting_hyperplane(const HalfSpace& space) const {
			size_t passes = 0;
			size_t fails = 0;
			for (const Vector& v : points)
			{
				HalfSpace::ContainmentState state = space.containment(v);

				if (state == HalfSpace::ContainmentState::strongly_in)
				{
					passes++;
				}
				else if (state == HalfSpace::ContainmentState::strongly_out)
				{
					fails++;
				}


				if (fails != 0 && passes != 0)
				{
					return false;
				}
			}
			return fails+passes!=size();
		}
		bool is_neccesary_supporting_hyperplane(const HalfSpace& space) const {
			size_t passes = 0;
			size_t fails = 0;
			for (const Vector& v : points)
			{
				HalfSpace::ContainmentState state = space.containment(v);

				if (state==HalfSpace::ContainmentState::strongly_in)
				{
					passes++;
				}
				else if (state == HalfSpace::ContainmentState::strongly_out)
				{
					fails++;
				}
				

				if (fails != 0 && passes != 0)
				{
					return false;
				}
			}
			return dim()<=size()-passes-fails;
		}

		bool contained_in(const HalfSpace& space) const {
			for (const Vector& v : points)
			{
				if (!space.contains(v))
				{
					return false;
				}
			}
			return true;
		}
		size_t size() const {
			return points.size();
		}
		Vector point_in() {
			Vector total(dim());
			for (size_t i = 0; i < size(); i++)
			{
				total += points[i];
			}
			return total / size();
		}
		std::vector<Vector> select(const std::vector<size_t>& indices) const {
			std::vector<Vector> vertices;
			vertices.reserve(indices.size());
			for (size_t i : indices) {
				vertices.push_back(points[i]);
			}
			return vertices;
		}
	};
	struct Frep {
		void translate(const Vector& p) {
			for (HalfSpace& h : *this)
			{
				h.translate(p);
			}
		}
		Frep() {

		}
		Vector project_onto(Vector pnt,size_t break_count) const {
			size_t times = 0;

			while (true) {
				const HalfSpace* worst = nullptr;
				precise worst_violation = 0;

				for (const HalfSpace& h : planes) {
					precise violation = h.signed_dist(pnt);

					if (violation > worst_violation) {
						worst_violation = violation;
						worst = &h;
					}
				}

				if (worst_violation < 1e-8)
				{
					break;
				}
				pnt = worst->project_onto_halfspace(pnt);
				if (times>= break_count) {
					break;
				}
				times++;

			}

			
			return pnt;
		}
		
		Vector chebeshev_center(size_t iterations=100) const {
			Vector pnt = project_onto(Vector(dim()),1000);
			
			precise slack = 0;
			for(int i=0;i< iterations;i++){

				std::optional<HalfSpace> space;
				for (const HalfSpace& h : planes)
				{
					//greatest error plane 
					if (!space || h.signed_dist(pnt) >= space.value().signed_dist(pnt))
					{
						slack = h.signed_dist(pnt);
						space = h;
					}
				}

				if (space)
				{
					HalfSpace max= space.value();
					
					Vector dir = Vector::random(dim());
					while(Vector::dot_product(max.normal,dir)<=0)
					{
						dir = Vector::random(dim());
					}
					precise move_dist = -std::numeric_limits<precise>().infinity();
					for (const HalfSpace& h : planes)
					{
						precise h_m = Vector::dot_product(dir, max.normal);
						precise h_o = Vector::dot_product(dir, h.normal);
						//our_dp*x+our_amt=other_dp*x+other_amt
						precise normal_ratio= h_m-h_o;

						if (abs(normal_ratio)>=1e-5)
						{
							precise dp_max= max.signed_dist(pnt);
							precise dp_other= h.signed_dist(pnt);
							precise t = -(dp_max-dp_other) / (normal_ratio);
							
							//get this working later
							if (Vector::dot_product(dir, h.normal)<0)
							{

								move_dist = std::max(move_dist, t);
							}
								
						}
							
						
					}
					if (move_dist==-std::numeric_limits<precise>().infinity())
					{
						return pnt;
					}
					pnt += dir*move_dist;
				}


			}
			return pnt;
		}

		
			
		std::optional<precise> hit_time(const Ray& r) const {
			std::optional<precise> min_time;
			for (const HalfSpace& space : planes)
			{
				std::optional<precise> d = space.hit_time(r);
				if (d.has_value() && (!min_time.has_value() || d.value() <= min_time.value()))
				{
					min_time = d;
				}
			}
			return min_time;
		}
		std::optional<Vector> hit(const Ray& r) const {
			std::optional<precise> ht = hit_time(r);
			if (!ht)
			{
				return std::nullopt;
			}
			return r.along(ht.value());
		}
		bool corner(const Vector& point) const
		{
			size_t bnd_count = 0;
			for (const HalfSpace& h : planes)
			{
				if (h.suffieciently_close(point))
				{
					bnd_count++;
				}
				else
				{
					if (!h.contains_exactly(point))
					{
						return false;
					}
				}

			}
			return (bnd_count >= dim());
		}
		bool boundry(const Vector& point) const {
			bool bnd = false;
			for (const HalfSpace& h : planes)
			{
				if (h.suffieciently_close(point))
				{
					bnd = true;
				}
				else
				{
					if (!h.contains_exactly(point))
					{
						return false;
					}
				}

			}
			return bnd;
		}
		using iterator = std::vector<HalfSpace>::iterator;
		iterator begin() {
			return planes.begin();
		}
		iterator end() {
			return planes.end();
		}
		using const_iterator = std::vector<HalfSpace>::const_iterator;
		const_iterator begin() const {
			return planes.begin();
		}
		const_iterator end() const {
			return planes.end();
		}
		std::vector<HalfSpace> planes;
		
		size_t size() const {
			return planes.size();
		}
		size_t dim() const {
			return planes.size() == 0 ? 0:planes[0].dim();
		}
		std::vector<HalfSpace> faces_on(const Vector& point) const {
			std::vector<HalfSpace> plane_list;
			for (const HalfSpace& space: planes)
			{
				if (space.suffieciently_close(point))
				{
					plane_list.push_back(space);
				}

			}
			return plane_list;
		}
		bool contains(const Vector& point) const {
			for (const HalfSpace& h : planes)
			{
				if (!h.contains(point))
				{
					return false;
				}
			}
			return true;
		}

		void add_unchecked(const HalfSpace& space)
		{
			planes.push_back(space);
		}
		void add(const HalfSpace& space)
		{
			for (const auto& old : planes)
			{
				precise angle = Vector::dot_product(old.normal, space.normal);
				if (boost::multiprecision::abs(angle - 1) <= 1e-8 &&
					boost::multiprecision::abs(old.bound - space.bound) <= 1e-8)
				{
					return;
				}
			}

			planes.push_back(space);
		}
	};
	inline void combinations_helper(
		int n,
		int k,
		int start,
		std::vector<size_t>& current,
		std::vector<std::vector<size_t>>& result)
	{
		if (current.size() == k) {
			result.push_back(current);
			return;
		}

		for (int i = start; i <= n - (k - current.size()); ++i) {
			current.push_back(i);
			combinations_helper(n, k, i + 1, current, result);
			current.pop_back();
		}
	}

	inline std::vector<std::vector<size_t>> combinations(int n, int k)
	{
		if (n < k)
		{
			return { };
		}
		std::vector<std::vector<size_t>> result;
		std::vector<size_t> current;
		combinations_helper(n, k, 0, current, result);
		return result;
	}
	std::optional<Frep> to_frep(const Vrep& vrep) {
		Frep rep;
		if (!vrep.full_span(100))
		{
			return std::nullopt;
		}
		std::vector<std::vector<size_t>> combines{ combinations(vrep.size(), vrep.dim()) };
		for (std::vector<size_t> inds : combines) {
			std::vector<Vector> points;
			for (size_t ind : inds) {
				points.push_back(vrep.points[ind]);
			}
			auto space = space_from_vectors(points);
			if (space)
			{


				if (!vrep.contains_all(space.value()))
				{
					space.value()= -space.value();
				}
				if (vrep.is_neccesary_supporting_hyperplane(space.value())) {
					rep.add(space.value());
				}
			}
		}

		return rep;
	}
	bool bounded(Frep frep){
		if (frep.size()<=frep.dim())
		{
			return false;
		}
		Vrep dual;
		frep.translate(-frep.chebeshev_center());
		for (auto& h : frep) {
			dual.add(h.normal / h.bound);
		}
		std::optional<Frep> dual_rep= to_frep(dual);
		Vector orgin = Vector(frep.dim());
		return dual_rep&& dual_rep.value().contains(orgin)&& !dual_rep.value().boundry(orgin);
	}

	std::optional<Vrep> to_vrep(const Frep& frep) {

		if (!bounded(frep))
		{
			return std::nullopt;
		}
		Vrep rep;

		//std::cout << std::format("{}", frep.planes);
		std::vector<std::vector<size_t>> combines{ combinations(frep.size(), frep.dim()) };
		for (std::vector<size_t> inds : combines) {
			std::vector<HalfSpace> planes;
			for (size_t ind : inds) {
				planes.push_back(frep.planes[ind]);
			}
			std::optional<Vector> pnt = GaussianMatrix(planes).solve();
			if (pnt&&frep.boundry(*pnt)) {
				rep.add(*pnt);
			}
		}

		return rep;
	}	
	
	precise sphere_sa(size_t n) {
		return (pow(3.1415926, n / 2.0)) / tgamma(1 + (n / 2.0));
	}
	precise radial_volume_apx(const Frep& frep, size_t amt) {
		Vector pnt = frep.chebeshev_center();
		precise value = 0;
		for (size_t i = 0; i < amt; i++)
		{
			Ray look = Ray(pnt, Vector::random(pnt.dim()));
			precise val = frep.hit_time(look).value();
			value += pow(val, frep.dim()) / amt;
		}
		value *= sphere_sa(frep.dim());
		return value;
	}
	template<typename T>
	std::vector<T> unordered_intersection(
		const std::vector<T>& a,
		const std::vector<T>& b)
	{
		std::vector<T> result;

		std::unordered_set<T> set_a(
			a.begin(),
			a.end()
		);

		for (const T& item : b)
		{
			if (set_a.contains(item))
			{
				result.push_back(item);
			}
		}

		return result;
	}
	struct Incremental {
		size_t plane_count() const {

			return faces.size();
		}
		Vrep points;
		Frep faces;
		size_t dim() {
			return std::max(points.dim(),faces.dim());
		}
		
		bool body() {
			return points.size() >dim()&& faces.size() > faces.dim();

		}
		void well_defined_invariant() {
			for (const HalfSpace& space : faces)
			{
				if (!points.is_neccesary_supporting_hyperplane(space)) {
					
					print("points{},\nface{}", points.points, space);

					print("{}", points.points_on(space));
					points.is_neccesary_supporting_hyperplane(space);
					
					throw std::logic_error("all planes must be supporting");
				}
			}
			for (const Vector& point:points)
			{
				if (!faces.corner(point))
				{
					throw std::logic_error("all points must be on the boundry");
				}
			}
		}
		bool contains(const Vector& point) const {
			return faces.contains(point);
		}
		bool contained_in(const HalfSpace& space) const {
			return points.contained_in(space);
		}
		void add(const HalfSpace& space) {

			if (!body())
			{
				faces.add(space);

				points= to_vrep(faces).value_or(Vrep());

				if (body())
				{
					faces= to_frep(points).value();
				}
				return;
			}
			if (contained_in(space)) {
				return;
			}
			Vrep in;
			Vrep out;
			for (const Vector& v:points)
			{
				if (space.contains(v))
				{
					in.add_unchecked(v);
				}
				else
				{
					out.add_unchecked(v);
				}
			}
			points = in;
			for (const Vector& lost : out)
			{
				auto p1 = faces.faces_on(lost);
				for (const Vector& kept : in)
				{
					
					std::vector<HalfSpace> ihl;
					for (const HalfSpace& face:p1)
					{
						if (face.suffieciently_close(kept))
						{
							ihl.push_back(face);
						}
					}
					if (ihl.size() >= dim() - 1)
					{
						if (ihl.size() >= dim())
						{
							ihl.erase(ihl.begin() + dim() - 1, ihl.end());

						}
						ihl.push_back(space);
						std::optional<Vector> pnt=GaussianMatrix(ihl).solve();
						if (!pnt)
						{
							continue;
						}
						if (!contains(pnt.value()))
						{
							throw std::logic_error("point must be contained");
						}
							points.add(pnt.value());
							for (const HalfSpace& s : ihl)
							{
								if (!s.suffieciently_close(pnt.value()))
								{
									throw std::logic_error("error");
								}
							}
							for (const HalfSpace& h : ihl) {
								precise dist = boost::multiprecision::abs(Vector::dot_product(h.normal, *pnt) - h.bound);
								
							}
					}
				}
			}
			faces.add(space);
			Frep kept_faces;
			for (const HalfSpace& space : faces)
			{
				//only check those at risk;
				if (out.is_supporting_hyperplane(space))
				{
					if (!points.is_neccesary_supporting_hyperplane(space))
					{
						continue;
					}
				}
					kept_faces.add_unchecked(space);
				
			}
			faces = kept_faces;
			well_defined_invariant();
		}

		void add(const Vector& point) {
			if (!body())
			{

				points.add(point);
				faces = to_frep(points).value_or(Frep());
				if (body())
				{
					points = to_vrep(faces).value_or(points);
				}
				return;
			}

			if (contains(point)) {
				return;
			}
			Frep contained;
			Frep lost_rep;
			for (const HalfSpace& space : faces.planes)
			{
				if (space.contains(point))
				{
					contained.add(space);
				}
				else
				{
					lost_rep.add(space);
				}
			}

			faces = contained;
			std::vector<Vector> new_points;
			std::vector<Vector> rip;

			for (const HalfSpace& lost : lost_rep.planes)
			{
				auto p1 = points.points_on(lost);
				for (const HalfSpace& kept : contained.planes)
				{
				

					std::vector<Vector> f_p;
					for (const Vector& point_in: p1)
					{
						if (kept.suffieciently_close(point_in))
						{
							f_p.push_back(point_in);
						}
					}
					if (f_p.size()>=dim()-1)
					{
						if (f_p.size() >= dim())
						{
							f_p.erase(f_p.begin() + dim() - 1, f_p.end());

						}
						
						f_p.push_back(point);
						auto new_space = space_from_vectors(f_p);
						if (new_space)
						{
							for (const Vector& v : f_p)
							{
								if (!new_space.value().suffieciently_close(v))
								{
									new_space = space_from_vectors(f_p);
									print("{}", f_p);
									throw std::logic_error("a");
								}

							};

							add_face_trivial(new_space.value());
						}
						else {
							int l = 4;
						}
					}
				}
			}
			points.add(point);
			Vrep kept_points;
			Vrep lost_points;
			for (const Vector& pnt : points)
			{
				if (lost_rep.boundry(pnt))
				{
					if (!faces.corner(pnt))
					{
						lost_points.add_unchecked(pnt);
						continue;
					}
				}
				kept_points.add_unchecked(pnt);

			}
			points = kept_points;

			well_defined_invariant();
		}
		void add_face_trivial(HalfSpace space) {
			if (!points.contains_all(space))
			{
				space = (-space);
			}
			if (!points.contains_all(space))
			{
				print("points{}\n faces{}", points.points, faces.planes);
				print("generators{}", points.points_on(space));
				print("{}", space);
				for (Vector v:points)
				{

					if (!space.contains(v))
					{
					print("{}", v);
					print("dist{}", space.dist(v));
					}
				}

				throw std::logic_error("how");
			}
			faces.add(space);

		}
		precise support(const Vector& direction) const {
			return points.support(direction);
		}
		static std::optional<Incremental> try_build(Vrep point_list) {
			Incremental built;
			if (point_list.dim() == 0)
			{
				return Incremental{};

			}
			if (point_list.size() <= point_list.dim())
			{
				throw std::logic_error("no sub simplexes");
			}
			for (Vector pnt : point_list)
			{
				built.add(pnt);
			}
			return built;
		}static std::optional<Incremental> try_build(const Frep& plane_list) {
			Incremental built;
			if (plane_list.dim() == 0)
			{
				return Incremental{};

			}
			if (plane_list.size() <= plane_list.dim())
			{
				throw std::logic_error("no sub simplexes");
			}
			for (HalfSpace plane: plane_list)
			{
				built.add(plane);
			}
			return built;
		}

		precise volume() {
			
			well_defined_invariant();
			if (dim() == 0)
			{
				return 1;
			}
			Vector center = points.point_in();
			
			precise volume = 0;
			for (const HalfSpace& face : faces.planes)
			{


				std::vector<Vector> points_on_plane = points.points_on(face);
				Span spn = face.span();
				Vrep unleashed;
				for (Vector pnt : points_on_plane)
				{
					unleashed.add(spn.to_basis(pnt));
					
				}
				precise dist = face.dist(center);
			
				
				std::optional<Incremental> built = Incremental::try_build(unleashed);
				if (built)
				{
					volume += built->volume() * dist / dim();
				}
			}
			
			well_defined_invariant();
			return volume;
		}
	};
	struct SupportFunction {
		virtual size_t dim() const = 0;
		virtual precise operator()(Vector amt) const = 0;
		HalfSpace half_space_from(const Vector& amt) const {
			return HalfSpace(amt, this->operator()(amt));
		}
	
		Incremental f_rep(size_t size) const {
			Incremental p{};


			while (p.plane_count() < size)
			{
				std::vector<Vector> unique_candidates;
				size_t unique_direction_count = 10;
				// Layer 1: create 20 diverse candidates
				while (unique_candidates.size() < unique_direction_count)
				{
					unique_candidates.push_back(Vector::random(dim()));
				}


				Vector best;
				precise max_discrepancy = -std::numeric_limits<precise>::infinity();

				for (const Vector& candidate : unique_candidates)
				{
					precise discrepancy =abs(this->operator()(candidate) - p.support(candidate));
					if (discrepancy > max_discrepancy)
					{
						max_discrepancy = discrepancy;
						best = candidate;
					}
				}
				p.add(half_space_from(best));
			}
			return p;
		}
		bool contains(Vector pnt, size_t evals) const {
			for (size_t i = 0; i < evals; i++)
			{
				
				Vector next = Vector::random(dim());
				if (!half_space_from(next).contains(pnt))
				{
					return false;
				}
			}
			return true;
		}

		precise volume(size_t cont) const {
			
			Incremental inc = f_rep(cont);
		/*size_t sed = seed();
			while (inc.plane_count()>7) {
				for (size_t i = 0; i < inc.plane_count(); i++)
				{
					Incremental unleashed;
			
					for (size_t u = 0;u < inc.plane_count();u++)
					{

						if (i != u)
						{
							unleashed.add(inc.faces.planes[u]);
						}
					}
					try {
						if (unleashed.faces.size()>40||bounded(unleashed.faces))
						{
							print("vol{}", unleashed.volume());
						}
					}
					catch (std::logic_error& fail) {
						if (unleashed.body())
						{
							inc = unleashed;
						}
					}
				}
				
			}

			print("points{}\nfaces{}", inc.points.points, inc.faces.planes);
			*/
			return inc.volume();
		}

	};
	struct Sphere :SupportFunction {
		size_t dimention =4;
		size_t dim() const {
			return dimention;
		}

		precise operator()(Vector amt) const {
			return 1;
		}
	};
	Frep make_cube_frep()
	{
		Frep cube;

		for (size_t i = 0; i < 3; i++)
		{
			Vector n(3);
			n[i] = 1;
			cube.add(HalfSpace(n, 1));

			n[i] = -1;
			cube.add(HalfSpace(n, 1));
		}

		return cube;
	}
	void test_3d_micro_cluster_failure() {
		std::vector<Vector> points = {
			Vector({-4.0715805488172,     3.12477044164411,   -0.685762622210356}),
			Vector({ 0.47260663030706,    2.09108820497028,   -1.42089973493547}),
			Vector({ 2.6483251770572,     0.371189105956377,   0.727897650238536}),
			Vector({-0.00544572097927758,-0.261263158616372, -0.272989611968719}), // Micro-pair A
			Vector({-0.00544561832631156,-0.261263148637472, -0.272989643073816})  // Micro-pair B
		};

		Incremental hull;
		for (size_t i = 0; i < points.size(); ++i) {
			std::cout << "Inserting point [" << i << "]...\n";
			
			if (i==4)
			{
				int l = 4;
			}
			hull.add(points[i]);
			
		}
	}
	int main()
	{
		print("vol{}",Sphere().volume(200));
	}

	// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
	// Debug program: F5 or Debug > Start Debugging menu

	// Tips for Getting Started: 
	//   1. Use the Solution Explorer window to add/manage files
	//   2. Use the Team Explorer window to connect to source control
	//   3. Use the Output window to see build output and other messages
	//   4. Use the Error List window to view errors
	//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
	//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
