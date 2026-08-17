#include "half_space.h"
#pragma once
template<typename Number>
struct Frep {

	using half_space_type = HalfSpace< Number>;
	using vec_type = Vector< Number>;
	void translate(const vec_type& p) {
		for (half_space_type& h : *this)
		{
			h.translate(p);
		}
	}
	Frep() {

	}
	vec_type project_onto(vec_type pnt, size_t break_count) const {
		size_t times = 0;

		while (true) {
			const half_space_type* worst = nullptr;
			Number worst_violation = 0;

			for (const half_space_type& h : planes) {
				Number violation = h.signed_dist(pnt);

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
			if (times >= break_count) {
				break;
			}
			times++;

		}


		return pnt;
	}

	vec_type chebeshev_center(size_t iterations = 100) const {
		vec_type pnt = project_onto(vec_type(dim()), 1000);

		Number slack = 0;
		for (int i = 0;i < iterations;i++) {

			std::optional<half_space_type> space;
			for (const half_space_type& h : planes)
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
				half_space_type max = space.value();

				vec_type dir = vec_type::random(dim());
				while (vec_type::dot_product(max.normal, dir) <= 0)
				{
					dir = vec_type::random(dim());
				}
				Number move_dist = -std::numeric_limits<Number>().infinity();
				for (const half_space_type& h : planes)
				{
					Number h_m = vec_type::dot_product(dir, max.normal);
					Number h_o = vec_type::dot_product(dir, h.normal);
					//our_dp*x+our_amt=other_dp*x+other_amt
					Number normal_ratio = h_m - h_o;

					if (abs(normal_ratio) >= 1e-5)
					{
						Number dp_max = max.signed_dist(pnt);
						Number dp_other = h.signed_dist(pnt);
						Number t = -(dp_max - dp_other) / (normal_ratio);

						//get this working later
						if (vec_type::dot_product(dir, h.normal) < 0)
						{

							move_dist = std::max(move_dist, t);
						}

					}


				}
				if (move_dist == -std::numeric_limits<Number>().infinity())
				{
					return pnt;
				}
				pnt += dir * move_dist;
			}


		}
		return pnt;
	}



	std::optional<Number> hit_time(const Ray<Number>& r) const {
		std::optional<Number> min_time;
		for (const half_space_type& space : planes)
		{
			std::optional<Number> d = space.hit_time(r);
			if (d.has_value() && (!min_time.has_value() || d.value() <= min_time.value()))
			{
				min_time = d;
			}
		}
		return min_time;
	}
	std::optional<vec_type> hit(const Ray<Number>& r) const {
		std::optional<Number> ht = hit_time(r);
		if (!ht)
		{
			return std::nullopt;
		}
		return r.along(ht.value());
	}

	bool corner(const vec_type& point) const
	{
		size_t bnd_count = 0;
		for (const half_space_type& h : planes)
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
	bool boundry(const vec_type& point) const {
		bool bnd = false;
		for (const half_space_type& h : planes)
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
	using iterator = std::vector<half_space_type>::iterator;
	iterator begin() {
		return planes.begin();
	}
	iterator end() {
		return planes.end();
	}
	using const_iterator = std::vector<half_space_type>::const_iterator;
	const_iterator begin() const {
		return planes.begin();
	}
	const_iterator end() const {
		return planes.end();
	}
	std::vector<half_space_type> planes;

	size_t size() const {
		return planes.size();
	}
	size_t dim() const {
		return planes.size() == 0 ? 0 : planes[0].dim();
	}
	std::vector<half_space_type> faces_on(const vec_type& point) const {
		std::vector<half_space_type> plane_list;
		for (const half_space_type& space : planes)
		{
			if (space.suffieciently_close(point))
			{
				plane_list.push_back(space);
			}

		}
		return plane_list;
	}
	bool contains(const vec_type& point,double eps=1e-8) const {
		for (const half_space_type& h : planes)
		{
			if (!h.contains(point,eps))
			{
				return false;
			}
		}
		return true;
	}

	void add_unchecked(const half_space_type& space)
	{
		planes.push_back(space);
	}
	void add(const half_space_type& space)
	{
		for (const auto& old : planes)
		{
			Number angle = vec_type::dot_product(old.normal, space.normal);
			if (abs(angle - 1) <= 1e-8 &&
				abs(old.bound - space.bound) <= 1e-8)
			{
				return;
			}
		}

		planes.push_back(space);
	}
};