#include "random.h"

uint64_t seed64 = 17840520000464551223;
	inline void randomize_64(uint64_t& seed_64) {

		seed64 ^= seed64 << 13;
		seed64 ^= seed64 >> 7;
		seed64 ^= seed64 << 17;

	}
	size_t seed() {
		return seed64;
	}
double random()
{
	randomize_64(seed64);
	return static_cast<double>(seed64) / UINT64_MAX;
}

