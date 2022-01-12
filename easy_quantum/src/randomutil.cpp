#include <randomutil.h>

ns_easyquantum

std::vector<double> RandomEngine::operator()(size_t n) {
	std::vector<double> ret;
	ret.reserve(n);
	for (size_t i = 0u; i < n; ++i) {
		ret.push_back((*this)());
	}
	return ret;
}

double DefaultRandomEngine::operator()() {
	return dist(engine);
}

double XC_RandomEngine16807::operator()() {
	if (ia * (irandseed % iq) - ir * (irandseed / iq) >= 0)
		irandnewseed = ia * (irandseed % iq) - ir * (irandseed / iq);
	else
		irandnewseed = ia * (irandseed % iq) - ir * (irandseed / iq) + im;
	irandseed = irandnewseed;
	return (double)irandnewseed / im;
}

double _default_random_generator() {
	static XC_RandomEngine16807 engine;
	return engine();
}

RandomEngine** make_engines(size_t shots) {
	RandomEngine** rngs = (RandomEngine**)malloc(sizeof(RandomEngine*)*shots);
	for (size_t i = 0; i < shots; ++i) {
		rngs[i] = new default_random_generator(
			std::chrono::system_clock::now().time_since_epoch().count() + 665544 * i + 44
		);
	}
	return rngs;
}

RandomEngine** make_engines(size_t shots, std::vector<size_t> seeds) {
	RandomEngine** rngs = (RandomEngine**)malloc(sizeof(RandomEngine*)*shots);
	for (size_t i = 0; i < shots; ++i) {
		rngs[i] = new default_random_generator(seeds[i]);
	}
	return rngs;
}

void free_engines(RandomEngine** rngs, size_t shots) {
	for (size_t i = 0; i < shots; ++i) {
		delete rngs[i];
	}
}

ns_end;