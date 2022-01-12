#pragma once

#include <random>
#include <chrono>
#include <vector>
#include <global.h>

ns_easyquantum

class RandomEngine {
public:

	/* Generate a random number */
	virtual double operator()() = 0;

	/* Generate a group of numbers */
	std::vector<double> operator()(size_t n);
};

/* Default random engine. Uses <random> in STL. */
class DefaultRandomEngine : public RandomEngine {
private:
	std::mt19937_64 engine;
	std::uniform_real_distribution<double> dist;
public:
	DefaultRandomEngine() :
		engine(std::chrono::system_clock::now().time_since_epoch().count()),
		dist(0, 1) { }

	DefaultRandomEngine(uint64_t seed)
		:engine(seed) {	}

	/* Generate a random number */
	double operator()();
};

/* 16807 Random number generator. Written by XC. */
class XC_RandomEngine16807 : public RandomEngine {
	int irandseed = 0;
	int ia = 16807;
	int im = 2147483647;
	int iq = 127773;
	int ir = 2836;
	int irandnewseed = 0;
public:
	XC_RandomEngine16807() {
		irandseed = (int)std::chrono::system_clock::now().time_since_epoch().count();
	}
	XC_RandomEngine16807(int64_t _seed)
		: irandseed((int)_seed)
	{ }

	double operator()();
};

/* Use default random number generator, outputs a random number*/
double _default_random_generator();

using default_random_generator = XC_RandomEngine16807;

/* Make several engines */
RandomEngine** make_engines(size_t shots);

/* Make several engines with various seeds */
RandomEngine** make_engines(size_t shots, std::vector<size_t> seeds);

/* Free engines created by make_engines */
void free_engines(RandomEngine** rngs, size_t shots);

ns_end