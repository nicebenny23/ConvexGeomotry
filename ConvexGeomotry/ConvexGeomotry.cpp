#include "Incremental.h"
using CommonType = double;
	struct SupportFunction {
		virtual size_t dim() const = 0;
		virtual CommonType operator()(Vector<CommonType> amt) const = 0;
		HalfSpace<CommonType> half_space_from(const Vector<CommonType>& amt) const {
			return HalfSpace<CommonType>(amt, this->operator()(amt));
		}
	
		Incremental<CommonType> f_rep(size_t size) const {
			Incremental<CommonType> p{};


			while (p.plane_count() < size)
			{
				std::vector<Vector<CommonType>> unique_candidates;
				size_t unique_direction_count = 10;
				// Layer 1: create 20 diverse candidates
				while (unique_candidates.size() < unique_direction_count)
				{
					unique_candidates.push_back(Vector<CommonType>::random(dim()));
				}


				Vector<CommonType> best;
				CommonType max_discrepancy = -std::numeric_limits<CommonType>::infinity();

				for (const Vector<CommonType>& candidate : unique_candidates)
				{
					CommonType discrepancy =abs(this->operator()(candidate) - p.support(candidate));
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
		bool contains(Vector<CommonType> pnt, size_t evals) const {
			for (size_t i = 0; i < evals; i++)
			{
				
				Vector<CommonType> next = Vector<CommonType>::random(dim());
				if (!half_space_from(next).contains(pnt))
				{
					return false;
				}
			}
			return true;
		}

		CommonType volume(size_t cont) const {
			
			Incremental<CommonType> inc = f_rep(cont);
		/*	size_t sed = seed();
			while (inc.plane_count()>7) {
				for (size_t i = 0; i < inc.plane_count(); i++)
				{
					Incremental<CommonType> unleashed;
			
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

		CommonType operator()(Vector<CommonType> amt) const {
			return 1;
		}
	};
	Frep<CommonType> make_cube_frep()
	{
		Frep<CommonType> cube;

		for (size_t i = 0; i < 3; i++)
		{
			Vector<CommonType> n(3);
			n[i] = 1;
			cube.add(HalfSpace< CommonType>(n, 1));

			n[i] = -1;
			cube.add(HalfSpace< CommonType>(n, 1));
		}

		return cube;
	}
	void test_3d_micro_cluster_failure() {
		
		std::vector< Vector< CommonType>> points = {
		Vector< CommonType>({1,0,0}),
		Vector< CommonType>({0,0,1}),
		Vector< CommonType>({ .5,.99999995,  .5}),
		Vector< CommonType>({ .499999994,1.0000000,  .5000000000}),
		Vector< CommonType>({ .500000000,1.00000005,  .499999999})

		};
		Incremental<CommonType> hull;
		for (size_t i = 0; i < points.size(); ++i) {
			std::cout << "Inserting point [" << i << "]...\n";
			if (i+1==points.size())
			{
				int l = 4;
			}
			hull.add(points[i]);
		print("vertices{}", hull.points.points);
		print("faces{}", hull.faces.planes);
			
		}
	}
	int main()
	{
		test_3d_micro_cluster_failure();
		while (true) {
			print("vol{}", Sphere().volume(10));
		}
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
