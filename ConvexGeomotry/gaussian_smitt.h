#include "half_space.h"
#pragma once

template<typename Number>
struct GaussianMatrix {
	GaussianMatrix(std::vector<HalfSpace<Number>> halfspaces) {
		for (size_t i = 0; i < halfspaces.size(); i++) {
			Vector next_normal = halfspaces[i].normal;

			next_normal.coords.push_back(halfspaces[i].bound);
			rows.push_back(std::move(next_normal));
		}
	}
	using vec_type = Vector<Number>;

	using half_space_type= HalfSpace<Number>;
	std::vector<vec_type > rows;
	std::optional<vec_type > solve() {
		std::vector<vec_type> original_rows = rows;
		for (size_t i = 0; i < rows.size(); i++)
		{
			size_t pivot = i;
			Number best = abs(rows[i][i]);
			//partial pivot
			for (size_t r = i; r < rows.size(); ++r)
			{
				Number val = abs(rows[r][i]);
				if (val >= best)
				{
					best = val;
					pivot = r;
				}
			}
			std::swap(rows[pivot], rows[i]);
			if (abs(rows[i][i]) <= 1e-8)
			{
				return std::nullopt;
			}
			rows[i] /= rows[i][i];
			vec_type row = rows[i];

			for (size_t j = i + 1; j < rows.size(); j++)
			{
				rows[j] -= row * rows[j][i];
			}
		}
		size_t dim = rows.size();
		vec_type res(dim);
		for (int i = static_cast<int>(dim) - 1; i >= 0; i--) {
			Number sum = rows[i][dim];
			for (size_t j = i + 1; j < dim; j++) {
				sum -= rows[i][j] * res[j];
			}
			res[i] = sum;
		}
		return res;
	}

};
//<=n vectors and returns orthagonal ones with the same span
template<typename Number = Real>
inline std::optional<std::pair<std::vector<Vector<Number>>, Number>> gramm_shmitt_process(const std::vector<Vector<Number>>& vectors) {
	Real scale = 1.0;
	std::vector<Vector<Real>> real_converted{ vectors.begin(),vectors.end() };
	for (size_t i = 0; i < real_converted.size(); i++) {
		Vector<Real>& v = real_converted[i];



		for (size_t j = 0; j < i; ++j) {
			v -= real_converted[j] * Vector<Real>::dot_product(v, real_converted[j]);
		}
		Real mag = v.magnitude();
		scale *= mag;

		if (mag >= 1e-8) {
			v /= mag;
		}
		else
		{
			return std::nullopt;
		}
	}

	return std::optional(std::pair<std::vector<Vector<Number>>, Number>({ std::vector<Vector<Number>>(real_converted.begin(),real_converted.end()), Number(scale) }));
}
template<typename Number = Real>
struct Span {
	Number determinant() const {
		return gramm_shmitt_process(points).transform([](auto&& value) {return value.second;}).value_or(0);
	}
	std::vector<Vector<Number>> points;
	size_t dim() const {
		return points.size();
	}

	Vector<Number> to_basis(const Vector<Number>& point) const {
		Vector<Number> pnt(dim());
		for (size_t i = 0; i < dim(); i++)
		{
			Number mag = points[i].magnitude();
			pnt[i] = Vector<Number>::dot_product(points[i], (point)) / (mag * mag);
		}
		return pnt;
	}
};
template<typename Number>
inline Span<Number> span(const HalfSpace<Number>& space) {
	if (space.dim() == 0)
	{
		return Span<Number>();
	}
	using vec_type = Vector<Real>;
	std::vector<vec_type> points;
	points.push_back(Vector<Real>(space.normal));

	vec_type center = vec_type(HalfSpace<Real>(space).project_onto_plane(vec_type(space.dim())));
	for (size_t i = 1; i < space.dim(); i++)
	{
		vec_type unit = center;
		unit[i] += 1;
		points.push_back(HalfSpace<Real>(space).project_onto_plane(unit) - center);
	}
	points = gramm_shmitt_process(points).transform([](auto&& val) {return val.first;}).value_or(std::vector<vec_type>());
	std::swap(points[0], points.back());
	points.pop_back();
	return Span<Number>{ std::vector<Vector<Number>>(points.begin(), points.end())};
}
//n vectors
template<typename Number>
inline std::optional<HalfSpace<Number>> space_from_vectors(const std::vector<Vector<Number>>& vectors) {
	if (vectors.empty() || vectors.size() != vectors[0].dim()) {
		throw std::logic_error("cannot construct a half space with this amount of vectors");
	}

	using half_space_type = HalfSpace<Number>;
	using vec_type = Vector<Number>;
	std::vector<Vector<Real>> edges;
	edges.reserve(vectors.size() - 1);
	for (size_t i = 1; i < vectors.size(); i++) {
		edges.push_back(Vector<Real>(vectors[i]) - Vector<Real>(vectors[0]));
	}

	auto basis = gramm_shmitt_process<Real>(edges);
	if (!basis) {
		return std::nullopt;
	}

	auto& orthogonal_basis = basis.value().first;

	for (size_t i = 0; i < vectors[0].dim(); i++) {
		Vector<Real> normal(vectors[0].dim());
		normal[i] = 1.0;

		for (int pass = 0; pass < 1; pass++) {
			for (const auto& b : orthogonal_basis) {
				normal -= b.project_onto_this_normalized_already(normal);
			}
		}

		if (normal.magnitude() >= 1e-1) {
			normal /= normal.magnitude(); 
			HalfSpace<Real> hs(normal, Vector<Real>::dot_product(Vector<Real>(vectors[0]), normal));
			return half_space_type(hs);
		}
	}
	throw std::logic_error("not possible");
}
