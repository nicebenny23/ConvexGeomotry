#include "vector.h"
#pragma once

template<typename Number=Real>
struct Ray {
	Vector<Number> start;
	Vector<Number> direction;
	size_t dim() const {
		return start.dim();
	}
	Ray(const Vector<Number>& pnt, const Vector<Number>& direction) :start(pnt), direction(direction.normal()) {

	}
	Vector<Number> normal() const {
		return direction.normal();
	}
	Vector<Number> along(Number dist) const {
		return start + direction * dist;
	}
};
using std::min;

enum class ContainmentState {
	boundry,
	strongly_in,
	strongly_out,
};
template<typename Number = Real>
struct HalfSpace {
	template<typename OtherNumber>
		requires std::constructible_from<Number, OtherNumber>
	explicit HalfSpace(const HalfSpace<OtherNumber>& other):normal(other.normal),bound(other.bound) {

	}
	using vec_type = Vector<Number>;
	vec_type normal;
	Number bound;
	HalfSpace(const vec_type& norm, Number bnd) :normal(norm), bound(bnd) {
		bound /= norm.magnitude();
		normal.normalize();
	}
	HalfSpace translated(const vec_type& p) const {
		return HalfSpace(normal, bound + vec_type::dot_product(p, normal));
	}
	void translate(const vec_type& p) {
		*this = translated(p);
	}
	HalfSpace operator-() const{
		return HalfSpace(-normal, -bound);
	}
	bool operator==(const HalfSpace& other) const = default;
	bool operator!=(const HalfSpace& other) const = default;
	bool contains(const vec_type& other,double eps=1e-8) const {
		return containment(other,eps)!=ContainmentState::strongly_out;
	}
	ContainmentState containment(const vec_type& other,double eps=1e-8) const {

		Number state = bound - vec_type::dot_product(normal, other);
		Number abs_dist = abs(state);
		if (abs(state) < eps)
		{
			return ContainmentState::boundry;
		}
		if (state >= 0)
		{
			return ContainmentState::strongly_in;
		}
		return ContainmentState::strongly_out;
	}
	bool contains_exactly(const vec_type& other) const {
		return vec_type::dot_product(normal, other) <= bound;
	}
	bool suffieciently_close(const vec_type& other,double eps=1e-8) const {
		return containment(other, eps) == ContainmentState::boundry;
	}
	vec_type project_onto_plane(const vec_type& other) const {
		return other + normal * (bound - vec_type::dot_product(normal, other));
	}
	vec_type project_onto_halfspace(const vec_type& other) const {
		return other + normal * min(Number(0.0), (bound - vec_type::dot_product(normal, other)));
	}
	size_t dim() const {
		return normal.dim();
	}
	Number dist(const vec_type& point) const {
		return (project_onto_plane(point) - point).magnitude();
	}
	//negitive if in
	Number signed_dist(const vec_type& pnt) const {
		return (vec_type::dot_product(normal, pnt) - bound);
	}
	std::optional<Number> hit_time(const Ray<Number>& ray) const {
		Number ht = dist(ray.start) / vec_type::dot_product(normal, ray.direction);
		if (ht > 0)
		{
			return ht;
		}
		return std::nullopt;
	}
	
};
namespace std {
	template<typename Number>
	struct std::formatter<HalfSpace<Number>> {
		constexpr auto parse(std::format_parse_context& ctx) {
			return ctx.begin();
		}

		auto format(const HalfSpace< Number>& v, std::format_context& ctx) const {
			return std::format_to(ctx.out(), "(normal:{})({:.10})", v.normal, v.bound);
		}
	};

	template<typename Number>
	struct hash<HalfSpace<Number>> {
		size_t operator()(const HalfSpace< Number>& v) const noexcept {
			size_t seed = 0;
			seed = std::hash<Vector<Number>>()(v.normal);
			seed ^= std::hash<Number>()(v.bound) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			return seed;
		}
	};
}