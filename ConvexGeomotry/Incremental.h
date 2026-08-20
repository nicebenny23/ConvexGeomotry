#include "trivial_construction.h"
#include "id.h"
#pragma once
template<typename Number>
Number sphere_sa(size_t n) {
	return (pow(3.1415926, n / 2.0)) / tgamma(1 + (n / 2.0));
}
template<typename Number>
Number radial_volume_apx(const Frep< Number>& frep, size_t amt) {
	Vector< Number> pnt = frep.chebeshev_center();
	Number value = 0;
	for (size_t i = 0; i < amt; i++)
	{
		Ray< Number> look = Ray< Number>(pnt, Vector< Number>::random(pnt.dim()));
		Number val = frep.hit_time(look).value();
		value += pow(val, frep.dim()) / amt;
	}
	value *= sphere_sa< Number>(frep.dim());
	return value;
}

enum class polytope_degeneracy {
	Unbounded,
	Hyperplanar,
	Unconstructed,
};
struct face_marker {
};
using FacetID = stn::typed_id<face_marker>;
struct point_marker {

};
using VertexId = stn::typed_id<point_marker>;

template<typename Number>
struct Facet {
	bool operator==(const Facet& other) const {
		return other.id == id;
	}

	Facet(FacetID id, HalfSpace<Number> pln) :id(id), plane(pln) {

	}
	FacetID id;
	HalfSpace<Number> plane;
	std::unordered_set<VertexId> vertices;
	void add_vertex(VertexId vertex_id) {
		vertices.emplace(vertex_id);
	}
	bool prunable() {
		return vertices.size() == 0;
	}

	void erase_vertex(VertexId id) {
		vertices.erase(id);
	}
};

template<typename Number>
struct Vertex {
	bool operator==(const Vertex& other) const {
		return other.id == id;
	}
	Vertex(VertexId id, Vector<Number> pnt) :id(id), point(pnt) {

	}
	VertexId id;
	Vector<Number> point;
	void add_face(FacetID facet_id) {
		faces.emplace(facet_id);
	}
	std::unordered_set<FacetID> faces;
	bool prunable() {
		return faces.size() <= point.dim();
	}
	void erase_face(FacetID id) {
		faces.erase(id);
	}
};

template<typename T>
struct std::hash<Vertex<T>> {
	size_t operator()(Vertex<T> t) const noexcept {
		return std::hash<VertexId>{}(t.id);
	}
};
template<typename T>
struct std::hash<Facet<T>> {
	size_t operator()(Facet<T> t) const noexcept {
		return std::hash<FacetID>{}(t.id);
	}
};

template<typename Number>
struct Incremental {
	size_t plane_count() const {

		return facets.size();
	}

	using vec_type = Vector<Number>;
	using half_space_type = HalfSpace<Number>;
	using vertex_type = Vertex<Number>;

	using facet_type = Facet<Number>;
	using frep_type = Frep<Number>;
	using vrep_type = Vrep<Number>;
	std::unordered_map<VertexId, Vertex<Number>> vertices;

	std::unordered_map<FacetID, Facet<Number>> facets;
	std::optional<polytope_degeneracy> degeneracy = polytope_degeneracy::Unconstructed;
	size_t face_counter;
	size_t vertex_counter;
	std::optional<VertexId> insert_vertex_id(const Vector<Number>& v) {
		for (auto& vertex:vertices)
		{
			if (vec_type::distance_squared_le_than(vertex.second.point, v, 1e-12))
			{
				return std::nullopt;
			}

		}
		VertexId id = VertexId(vertex_counter++);
		vertices.emplace(id, Vertex<Number>(id, v));
		return id;
	}
	std::optional<FacetID> insert_facet_id(const HalfSpace<Number>& h) {
		for (const auto& face_t: facets)
		{
			const facet_type& face = face_t.second;
			Number angle = vec_type::dot_product(face.plane.normal, h.normal);
			if (abs(angle - 1) <= 1e-8 &&abs(face.plane.bound - h.bound) <= 1e-8)
			{
				return std::nullopt;
			}
		}
		FacetID id = FacetID(face_counter++);
		facets.emplace(id, Facet<Number>(id, h));
		return id;
	}
	size_t dim() {
		return std::max<size_t>(0, facets.begin()->second.plane.dim());
	}
	bool body() {
		return degeneracy == std::nullopt;

	}
	bool contains(const vec_type& point, double eps = 1e-8) const {
		for (auto& u : facets)
		{
			if (!u.second.plane.contains(point, eps))
			{
				return false;
			}
		}
		return true;
	}
	bool contained_in(const half_space_type& space, double eps = 1e-8) const {
		for (auto& u : vertices)
		{
			if (!space.contains(u.second.point, eps))
			{
				return false;
			}
		}
		return true;
	}
	void remove_vertex_connections(VertexId id) {
		vertex_type& v = vertices.at(id);
		for (FacetID face_id : v.faces)
		{
			facets.at(face_id).vertices.erase(id);
		}
	}
	void remove_facet_connections(FacetID id) {
		facet_type& f = facets.at(id);
		for (VertexId vertice : f.vertices)
		{
			vertices.at(vertice).faces.erase(id);
		}
	}

	void prune_facets() {
		std::unordered_map<FacetID, Facet<Number>> facets_kept;
		for (auto& f : facets)
		{
			if (!f.second.prunable())
			{
				facets_kept.emplace(f.first, f.second);
			}
			else {
				remove_facet_connections(f.second.id);

			}
		}

		facets = facets_kept;
	}

	void prune_vertices() {
		std::unordered_map<VertexId, Vertex<Number>> vertices_kept;
		for (auto& v: vertices)
		{
			if (!v.second.prunable())
			{
				vertices_kept.emplace(v.first, v.second);
			}
			else {
				remove_vertex_connections(v.second.id);
			}
		}

		vertices= vertices_kept;

	}
	Frep<Number> into_frep() {

		Frep<Number> faces;
		for (auto& u : facets)
		{
			faces.add(u.second.plane);
		}
		return faces;
	}
	

	Vrep<Number> into_vrep() {

		Vrep<Number> points;
		for (auto& u : vertices)
		{
			points.add(u.second.point);
		}
		return points;
	}void add(const half_space_type& new_space) {

		auto new_id_opt = insert_facet_id(new_space);
		if (!new_id_opt)
		{
			return;
		}
		FacetID new_id = new_id_opt.value();
		if (degeneracy != std::nullopt)
		{
			if (bounded(into_frep()))
			{

				auto combination_list = combinations(facets.size(), dim());
				for (std::vector<size_t> facet_id_list : combination_list) {
					std::vector<HalfSpace<Number>> half_spaces;
					for (size_t u : facet_id_list)
					{
						half_spaces.push_back(facets.at(FacetID(u)).plane);
					}
					auto pnt = GaussianMatrix(half_spaces).solve();
					if (pnt && contains(pnt.value()))
					{
						auto vertex_id_opt = insert_vertex_id(pnt.value());
						if (!vertex_id_opt)
						{
							continue;
						}
						auto vertex_id = vertex_id_opt.value();
						for (size_t facet_indices : facet_id_list)
						{
							vertices.at(vertex_id).add_face(FacetID(facet_indices));
							facets.at(FacetID(facet_indices)).add_vertex(vertex_id);
						}

					}
				}
				degeneracy = std::nullopt;
			}
			return;
		}
		std::unordered_set<Vertex<Number>> out;
		for (const auto& v : vertices)
		{
			if (!new_space.contains(v.second.point))
			{
				remove_vertex_connections(v.second.id);
				out.emplace(v.second);
			}
		}

		for (const Vertex<Number>& lost : out)
		{
			std::unordered_set<VertexId> potential_incidences;
			for (FacetID u : lost.faces)
			{
				//only kept are still in
				potential_incidences.insert_range(facets.at(u).vertices);
			}
			for (VertexId vertex_id : potential_incidences)
			{

				std::vector<FacetID> faces_in_both;

				vertex_type& vertex = vertices.at(vertex_id);
				faces_in_both.push_back(new_id);
				for (FacetID face_id : lost.faces)
				{
					if (vertex.faces.contains(face_id))
					{
						faces_in_both.push_back(face_id);
					}
				}
				if (faces_in_both.size() >= dim())
				{
					std::vector<HalfSpace<Number>> face_vector_in_both;
					for (FacetID id : faces_in_both)
					{
						if (face_vector_in_both.size() < dim())
						{
							face_vector_in_both.push_back(facets.at(id).plane);
						}
					}
					std::optional<vec_type> pnt = GaussianMatrix(face_vector_in_both).solve();
					if (pnt)
					{
						auto vertex_id_opt = insert_vertex_id(pnt.value());
						if (!vertex_id_opt)
						{
							continue;
						}
						auto new_vertex_id = vertex_id_opt.value();
						vertices.at(new_vertex_id).faces = std::unordered_set<FacetID>{ faces_in_both.begin(),faces_in_both.end() };
						for (FacetID touching : faces_in_both)
						{
							facets.at(touching).vertices.emplace(new_vertex_id);
						}
					}
				}

			}
		}
		for (Vertex v : out)
		{
			vertices.erase(v.id);
		}
		prune_facets();
		if (facets.size() <= dim())
		{

			degeneracy = polytope_degeneracy::Unbounded;
		}
	}

	void add(const vec_type& new_vertex) {

		auto new_id_opt = insert_vertex_id(new_vertex);
		if (!new_id_opt)
		{
			return;
		}
		VertexId new_id = new_id_opt.value();
		if (degeneracy != std::nullopt)
		{
			if (into_vrep().full_span())
			{

				auto combination_list = combinations(vertices.size(), dim());
				for (std::vector<size_t> vertex_id_list : combination_list) {
					std::vector<vec_type> points;
					for (size_t ind : vertex_id_list)
					{
						points.push_back(vertices.at(VertexId(ind)).point);
					}
					auto face= space_from_vectors(points);
					if (face && contained_in(face.value()))
					{
						auto face_id_opt=insert_facet_id(face.value());
						if (!face_id_opt)
						{
							continue;
						}
						auto face_id = face_id_opt.value();
						for (size_t vertex_indices : vertex_id_list)
						{
							facets.at(face_id).add_vertex(VertexId(vertex_indices));
							vertices.at(VertexId(vertex_indices)).add_face(face_id);
						}

					}
				}
				degeneracy = std::nullopt;
			}
			return;
		}
		std::unordered_set<facet_type> out;
		for (const auto& f: facets)
		{
			if (!f.second.plane.contains(new_vertex))
			{
				remove_facet_connections(f.second.id);
				out.emplace(f.second);
			}
		}

		for (const facet_type& lost : out)
		{
			std::unordered_set<FacetID> potential_incidences;
			for (VertexId u : lost.vertices)
			{
				//only kept are still in
				potential_incidences.insert_range(vertices.at(u).faces);
			}
			for (FacetID facet_id : potential_incidences)
			{
				std::vector<VertexId> points_in_both;

				facet_type& face= facets.at(facet_id);
				points_in_both.push_back(new_id);
				for (VertexId vertex_id : lost.vertices)
				{
					if (face.vertices.contains(vertex_id))
					{
						points_in_both.push_back(vertex_id);
					}
				}
				if (points_in_both.size() >= dim())
				{
					std::vector<Vector<Number>> point_vector_in_both;
					for (VertexId id : points_in_both)
					{
						if (point_vector_in_both.size() < dim())
						{
							point_vector_in_both.push_back(vertices.at(id).point);
						}
					}
					std::optional<half_space_type> space= space_from_vectors(point_vector_in_both);
					if (space)
					{
						auto facet_id_opt = insert_facet_id(space.value());
						if (!facet_id_opt)
						{
							continue;
						}
						auto new_facet_id = facet_id_opt.value();
						facets.at(new_facet_id).vertices= std::unordered_set<VertexId>{ points_in_both.begin(),points_in_both.end() };
						for (VertexId touching : points_in_both)
						{
							vertices.at(touching).add_face(new_facet_id);
						}
					}
				}

			}
		}
		for (Facet f : out)
		{
			facets.erase(f.id);
		}
		prune_facets();
		
	}
	

	Number support(const vec_type& direction) const {

		Number spt = -std::numeric_limits<Number>().infinity();
		for (const auto& v : vertices)
		{
			spt = std::max(spt, vec_type::dot_product(direction, v.second.point));
		}
		return spt;
	}
	Number volume() {

		if (dim() == 0)
		{
			return 1;
		}
		if (degeneracy == polytope_degeneracy::Hyperplanar)
		{
			return 0;
		}
		if (degeneracy == polytope_degeneracy::Unconstructed)
		{
			return 0;
		}
		if (degeneracy == polytope_degeneracy::Unbounded)
		{
			return std::numeric_limits<Number>().infinity();
		}
		/*
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

		return volume;
	*/
	}
};