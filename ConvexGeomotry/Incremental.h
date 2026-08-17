#include "trivial_construction.h"
#pragma once
template<typename Number>
Number sphere_sa(size_t n) {
	return (pow(3.1415926, n / 2.0)) / tgamma(1 + (n / 2.0));
}
template<typename Number>
Real radial_volume_apx(const Frep< Number>& frep, size_t amt) {
	Vector< Number> pnt = frep.chebeshev_center();
	Number value = 0;
	for (size_t i = 0; i < amt; i++)
	{
		Ray< Number> look = Ray< Number>(pnt, Vector< Number>::random(pnt.dim()));
		Real val = frep.hit_time(look).value();
		value += pow(val, frep.dim()) / amt;
	}
	value *= sphere_sa< Number>(frep.dim());
	return value;
}

enum class polytope_degeneracy{
	Unbounded,
	Hyperplanar,
	Unconstructed,
};
template<typename Number>
struct Incremental {
	size_t plane_count() const {

		return faces.size();
	}

	using vec_type = Vector<Number>;
	using half_space_type = HalfSpace<Number>;

	using frep_type = Frep<Number>;
	using vrep_type = Vrep<Number>;
	std::optional<polytope_degeneracy> degeneracy=polytope_degeneracy::Unconstructed;
	vrep_type points;
	frep_type faces;
	size_t dim() {
		return std::max(points.dim(), faces.dim());
	}
	bool body() {
		return degeneracy==std::nullopt;

	}
	void well_defined_invariant() {
		for (const half_space_type& space : faces)
		{
			if (!points.is_neccesary_supporting_hyperplane(space)) {

				print("points{},\nface{}", points.points, space);

				print("{}", points.points_on(space));
				points.is_neccesary_supporting_hyperplane(space);

				throw std::logic_error("all planes must be supporting");
			}
		}
		for (const vec_type& point : points)
		{
			if (!faces.corner(point))
			{
				throw std::logic_error("all points must be on the boundry");
			}
		}
	}
	bool contains(const vec_type& point,double eps=1e-8) const {
		return faces.contains(point,eps);
	}
	bool contained_in(const half_space_type& space, double eps = 1e-8) const {
		return points.contained_in(space,eps);
	}
	void add(const half_space_type& space) {

		if (!body())
		{
			
			faces.add(space);

			points = to_vrep(faces).value_or(points);

			if (points.size()!=0)
			{
				degeneracy = std::nullopt;
				faces = to_frep(points).value();
			}
			else {
				degeneracy = polytope_degeneracy::Unbounded;
			}
			return;
		}
		if (contained_in(space,1e-6)) {
			return;
		}
		vrep_type in;
		vrep_type out;
		for (const vec_type& v : points)
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
		for (const vec_type& lost : out)
		{
			auto p1 = faces.faces_on(lost);
			for (const vec_type& kept : in)
			{

				std::vector<half_space_type> ihl;
				for (const half_space_type& face : p1)
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
					std::optional<vec_type> pnt = GaussianMatrix(ihl).solve();
					if (!pnt)
					{
						continue;
					}
					if (!contains(pnt.value()))
					{
						throw std::logic_error("point must be contained");
					}
					points.add(pnt.value());
					for (const half_space_type& s : ihl)
					{
						if (!s.suffieciently_close(pnt.value()))
						{
							throw std::logic_error("error");
						}
					}
					for (const half_space_type& h : ihl) {
						Number dist = abs(vec_type::dot_product(h.normal, *pnt) - h.bound);

					}
				}
			}
		}
		faces.add(space);
		if (points.size()==0) {
			return;
		}
		if (points.contains_all(space) && points.contains_all(-space)) {
			degeneracy = polytope_degeneracy::Hyperplanar;
			faces = Frep<Number>();
			return;
		}

		
		frep_type kept_faces;
		for (const half_space_type& space : faces)
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
	}

	void add(const vec_type& point) {
		if (!body())
		{
			if (degeneracy == polytope_degeneracy::Unbounded)
			{
				throw std::logic_error("cannot add plane while unbounded");
			}
			points.add(point);
			faces = to_frep(points).value_or(frep_type());
			if (faces.size()!=0)
			{
				degeneracy = std::nullopt;
				points = to_vrep(faces).value_or(points);

			}
			else {
				degeneracy = polytope_degeneracy::Hyperplanar;
			}
		}

		if (contains(point,1e-6)) {
			return;
		}

		frep_type contained;
		frep_type lost_rep;
		for (const half_space_type& space : faces.planes)
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
		std::vector<vec_type> new_points;
		std::vector<vec_type> rip;

		for (const half_space_type& lost : lost_rep.planes)
		{

			auto p1 = points.points_on(lost);

			for (const half_space_type& kept : contained.planes)
			{


				std::vector<vec_type> f_p;
				for (const vec_type& point_in : p1)
				{
					if (kept.suffieciently_close(point_in))
					{
						f_p.push_back(point_in);
					}
				}
				if (f_p.size() >= dim() - 1)
				{
					if (f_p.size() >= dim())
					{
						f_p.erase(f_p.begin() + dim() - 1, f_p.end());

					}

					f_p.push_back(point);

					auto new_space = space_from_vectors(f_p);

					bool ss = false;
					if (new_space)
					{
						auto space = new_space.value();
						
						
						add_face_trivial(space);
					}
					else {
						int l = 4;
					}
				}
			}
		}
		points.add(point);
		vrep_type kept_points;
		vrep_type lost_points;
		for (const vec_type& pnt : points)
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

	}
	void add_face_trivial(half_space_type space)
	{
		if (!points.contains_all(space))
		{
			space = -space;
		}

		if (!points.contains_all(space))
		{
			print("points{}\n faces{}", points.points, faces.planes);
			print("generators{}", points.points_on(space));
			print("{}", space);

			// First: test using the existing Number-space.
			for (const vec_type& v : points)
			{
				Number fake_residual = (-space).signed_dist(v);

				if (fake_residual > Number(1e-8))
				{
					print("point = {}", v);
					auto on = points.points_on(space);
					if (on.size()>=dim())
					{

						// Reconstruct the plane from the points in Real.
						std::vector<Vector<Real>> real_generators{ on.begin(),on.end() };
						auto fake_again = HalfSpace<Real>(space);
						auto real_space = space_from_vectors(real_generators).value();
						Vector<Real> real_point(v);
						Real real_residual = fake_again.signed_dist(real_point);
						print("REAL signed residual = {:.15e}", real_residual);

					}
				}
			}

			throw std::logic_error("how");
		}

		faces.add(space);
	}
	Number support(const vec_type& direction) const {
		return points.support(direction);
	}
	static std::optional<Incremental> try_build(vrep_type point_list) {
		Incremental built;
		if (point_list.dim() == 0)
		{
			return Incremental{};

		}
		if (point_list.size() <= point_list.dim())
		{
			throw std::logic_error("no sub simplexes");
		}
		for (vec_type pnt : point_list)
		{
			built.add(pnt);
		}
		return built;
	}static std::optional<Incremental> try_build(const frep_type& plane_list) {
		Incremental built;
		if (plane_list.dim() == 0)
		{
			return Incremental{};

		}
		if (plane_list.size() <= plane_list.dim())
		{
			throw std::logic_error("no sub simplexes");
		}
		for (half_space_type plane : plane_list)
		{
			built.add(plane);
		}
		return built;
	}

	Number volume() {

		well_defined_invariant();
		if (dim()==0)
		{
			return 1;
		}
		if (degeneracy==polytope_degeneracy::Hyperplanar)
		{
			return 0;
		}
		if (degeneracy == polytope_degeneracy::Unconstructed)
		{
			return 0;
		}
		if (degeneracy==polytope_degeneracy::Unbounded)
		{
			return std::numeric_limits<Number>().infinity();
		}
		vec_type center = points.point_in();

		Number volume = 0;
		for (const half_space_type& face : faces.planes)
		{


			std::vector<vec_type> points_on_plane = points.points_on(face);
			Span spn = span(face);
			vrep_type unleashed;
			for (vec_type pnt : points_on_plane)
			{
				unleashed.add(spn.to_basis(pnt));

			}
			Number dist = face.dist(center);


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