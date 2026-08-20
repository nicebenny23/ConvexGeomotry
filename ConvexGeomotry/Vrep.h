#include "half_space.h"
#pragma once
template<typename Number>
struct Vrep {

	using half_space_type= HalfSpace< Number>;
	using vec_type = Vector < Number>;
	void add_unchecked(const vec_type& pnt)
	{
		points.push_back(pnt);
	}
	Number support(const vec_type& direction) const {

		Number spt = -std::numeric_limits<Number>().infinity();
		for (const vec_type& v : *this)
		{
			spt = std::max(spt, vec_type::dot_product(direction, v));
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
			std::vector<vec_type> edges;
			edges.reserve(d);
			std::vector<vec_type> chosen;
			while (chosen.size() < d + 1)
			{
				vec_type base = points[random() * points.size()];
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
	void translate(const vec_type& p) {
		for (vec_type& v : points)
		{
			v += p;
		}
	}
	using iterator = std::vector<vec_type>::iterator;
	iterator begin() {
		return points.begin();
	}
	iterator end() {
		return points.end();
	}
	using const_iterator = std::vector<vec_type>::const_iterator;
	const_iterator begin() const {
		return points.begin();
	}
	const_iterator end() const {
		return points.end();
	}
	size_t dim() const {
		return points.empty() ? 0 : points[0].dim();
	}
	bool redudant(const vec_type& point) {
		for (const vec_type& pnt : points)
		{
			if (vec_type::distance(point, pnt) <= 1e-8)
			{
				return true;
			}
		}
		return false;
	}
	void add(const vec_type& point) {
		for (const vec_type& pnt : points)
		{
			double eps = 1e-8;
			double eps_squared = eps * eps;
			if (vec_type::distance_squared_le_than(point, pnt,eps_squared))
			{
				return;
			}

		}
		points.push_back(point);
	}

	std::vector<vec_type> points;
	std::vector<vec_type> points_on(const half_space_type& space) const {
		std::vector<vec_type> point_list;
		for (const vec_type& v : points)
		{
			ContainmentState state = space.containment(v);

			if (space.suffieciently_close(v))
			{

				point_list.push_back(v);
			}
		}
		return point_list;
	}
	bool contains_all(half_space_type space) const {
		for (const vec_type& v : points)
		{
			if (!space.contains(v))
			{
				return false;
			}
		}
		return true;
	}
	bool is_supporting_hyperplane(const half_space_type& space) const {
		size_t passes = 0;
		size_t fails = 0;
		for (const vec_type& v : points)
		{
			ContainmentState state = space.containment(v);

			if (state == ContainmentState::strongly_in)
			{
				passes++;
			}
			else if (state == ContainmentState::strongly_out)
			{
				fails++;
			}


			if (fails != 0 && passes != 0)
			{
				return false;
			}
		}
		return fails + passes != size();
	}
	bool is_neccesary_supporting_hyperplane(const half_space_type& space) const {
		size_t passes = 0;
		size_t fails = 0;
		for (const vec_type& v : points)
		{
			ContainmentState state = space.containment(v);

			if (state == ContainmentState::strongly_in)
			{
				passes++;
			}
			else if (state == ContainmentState::strongly_out)
			{
				fails++;
			}


			if (fails != 0 && passes != 0)
			{
				return false;
			}
		}
		return dim() <= size() - passes - fails;
	}

	bool contained_in(const half_space_type& space,double eps=1e-8) const {
		for (const vec_type& v : points)
		{
			if (!space.contains(v,eps))
			{
				return false;
			}
		}
		return true;
	}
	size_t size() const {
		return points.size();
	}
	vec_type point_in() {
		vec_type total(dim());
		for (size_t i = 0; i < size(); i++)
		{
			total += points[i];
		}
		return total / size();
	}
	std::vector<vec_type> select(const std::vector<size_t>& indices) const {
		std::vector<vec_type> vertices;
		vertices.reserve(indices.size());
		for (size_t i : indices) {
			vertices.push_back(points[i]);
		}
		return vertices;
	}
};