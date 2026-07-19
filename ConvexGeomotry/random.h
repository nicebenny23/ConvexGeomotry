#pragma once

#include <stdint.h>
#include <math.h>
#include <random>
#include <numbers>
#pragma once
	double random();inline double gaussian() {
        double u1 = random();
        double u2 = random();

        return std::sqrt(-2.0 * std::log(u1)) *
            std::cos(2.0 * std::numbers::pi * u2);
    }