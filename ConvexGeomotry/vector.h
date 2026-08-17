#include <vector>
#include <iostream>
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
#include <boost/multiprecision/cpp_bin_float.hpp>
#pragma once
using QuadBackend =
boost::multiprecision::backends::cpp_bin_float<
	113,
	boost::multiprecision::backends::digit_base_2,
	void,
	std::int16_t,
	-16382,
	16383
>;

namespace boost::multiprecision::detail {

	template<>
	struct is_lossy_conversion<double, QuadBackend>
		: std::false_type
	{
	};

	template<>
	struct is_explicitly_convertible<double, QuadBackend>
		: std::true_type
	{
	};

} // namespace boost::multiprecision::detail

using Real = boost::multiprecision::number<
	QuadBackend,
	boost::multiprecision::et_on
>;
namespace std {
	template <>
	struct std::formatter<Real, char>
		: std::formatter<double, char>
	{
		auto format(const Real& value, std::format_context& ctx) const
		{
			return std::formatter<double, char>::format(
				value.convert_to<double>(), ctx
			);
		}
	};
}
#include <unordered_set>
template<typename ...Args>
void print(std::format_string<Args...> fmt, Args&&... args) {
	std::cout << std::format(fmt, std::forward<Args>(args)...) << '\n';
}
template<typename Number=Real>
struct Vector {
	template<typename OtherNumber>
		requires std::constructible_from<Number, OtherNumber>
	explicit Vector(const Vector<OtherNumber>& other)
		: coords(other.coords.begin(), other.coords.end())
	{
	}	std::vector<Number> coords;

	Vector() = default;

	explicit Vector(const std::vector<Number>& v)
		: coords(v) {
	}
	explicit Vector(std::vector<Number>&& v)
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

	static Number dot_product(const Vector& a, const Vector& b) {
		check_dimensions(a, b);

		Number sum = 0.0;
		for (size_t i = 0; i < a.dim(); i++)
			sum += a[i] * b[i];

		return sum;
	}
	//one when perfectly aligned
	static Number aligment0to1(const Vector& a, const Vector& b) {
		return (dot_product(a, b) + 1) / 2.0;

	}
	static Number distance(const Vector& a, const Vector& b) {
		return (a - b).magnitude();
	}
	static Number distance_squared(const Vector& a, const Vector& b) {
		return (a - b).magnitude_squared();

	}
	Number magnitude() const {

		using boost::multiprecision::sqrt;
		using std::sqrt;
		return sqrt(magnitude_squared());
	}
	Number magnitude_squared() const {
		return dot_product(*this, *this);
	}
	Vector normal() const {
		if (dim() == 0)
		{
			return *this;
		}
		return *this / magnitude();
	}
	Vector& normalize() {
		*this /= magnitude();
		return *this;
	}
	Number& operator[](size_t i) {
		return coords[i];
	}

	const Number& operator[](size_t i) const {
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


	Vector operator*(Number scalar) const {
		Vector result(dim());

		for (size_t i = 0; i < dim(); i++)
			result[i] = coords[i] * scalar;

		return result;
	}

	Vector& operator*=(Number scalar) {
		for (size_t i = 0; i < dim(); i++)
			coords[i] *= scalar;

		return *this;
	}

	Vector operator/(Number scalar) const {
		if (scalar == 0.0)
			throw std::invalid_argument("Division by zero.");

		return *this * (1.0 / scalar);
	}

	Vector& operator/=(Number scalar) {
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


namespace std {

	template<typename Number>
	struct std::formatter<Vector<Number>> {
		constexpr auto parse(std::format_parse_context& ctx) {
			return ctx.begin();
		}

		auto format(const Vector<Number> & v, std::format_context& ctx) const {
			auto out = ctx.out();
			*out++ = '(';

			for (size_t i = 0; i < v.dim(); i++)
			{
				if (i != 0)
				{
					*out++ = ',';
					*out++ = ' ';
				}

				out = std::format_to(out, "{:.9}", v[i]);
			}

			*out++ = ')';

			return out;
		}
	};

	template<typename Number >
	struct hash<Vector<Number>> {
		size_t operator()(const Vector< Number>& v) const noexcept {
			size_t seed = 0;

			for (Number x : v.coords) {
				size_t h = std::hash<Number>{}(x);

				// hash combine (boost-style)
				seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}

			return seed;
		}
	};
}