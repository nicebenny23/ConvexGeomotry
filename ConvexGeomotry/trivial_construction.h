#include "Vrep.h"
#include "Frep.h"
#include "gaussian_smitt.h"
namespace std{
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
template<typename Number>
std::optional<Frep< Number>> to_frep(const Vrep< Number>& vrep) {
	Frep< Number> rep;
	if (!vrep.full_span(100))
	{
		return std::nullopt;
	}
	std::vector<std::vector<size_t>> combines{ combinations(vrep.size(), vrep.dim()) };
	for (std::vector<size_t> inds : combines) {
		std::vector<Vector< Number>> points;
		for (size_t ind : inds) {
			points.push_back(vrep.points[ind]);
		}
		auto space = space_from_vectors(points);
		if (space)
		{


			if (!vrep.contains_all(space.value()))
			{
				space.value() = -space.value();
			}
			if (vrep.is_neccesary_supporting_hyperplane(space.value())) {
				rep.add(space.value());
			}
		}
	}

	return rep;
}
template<typename Number>
bool bounded(Frep<Number> frep) {
	if (frep.size() <= frep.dim())
	{
		return false;
	}
	Vrep<Number> dual;
	frep.translate(-frep.chebeshev_center());
	for (auto& h : frep) {
		if (h.bound==0)
		{
			return false;
		}
		dual.add(h.normal / h.bound);
	}
	std::optional<Frep< Number>> dual_rep = to_frep(dual);
	Vector< Number> orgin = Vector< Number>(frep.dim());
	return dual_rep && dual_rep.value().contains(orgin) && !dual_rep.value().boundry(orgin);
}
template<typename Number>
std::optional<Vrep<Number>> to_vrep(const Frep< Number>& frep) {

//	print("{}", frep.size());
	if (!bounded(frep))
	{
		return std::nullopt;
	}
	Vrep<Number> rep;

	//std::cout << std::format("{}", frep.planes);
	std::vector<std::vector<size_t>> combines{ combinations(frep.size(), frep.dim()) };
	for (std::vector<size_t> inds : combines) {
		std::vector<HalfSpace<Number>> planes;
		for (size_t ind : inds) {
			planes.push_back(frep.planes[ind]);
		}
		std::optional<Vector<Number>> pnt = GaussianMatrix(planes).solve();
		if (pnt && frep.boundry(*pnt)) {
			rep.add(*pnt);
		}
	}

	return rep;
}
