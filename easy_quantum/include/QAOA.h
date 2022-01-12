#pragma once

#include "quantum.h"
#include "Differential.h"

ns_easyquantum

struct pair_hash {
	template<class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& p) const {
		auto h1 = std::hash<T1>{}(p.first);
		auto h2 = std::hash<T2>{}(p.second);
		return h1 ^ h2;
	}
};

template<size_t n_vertices_, typename qtraits_t = default_qtraits>
class QAOA_Graph {
public:
	using qtraits = qtraits_t;
	static constexpr size_t n_vertices = n_vertices_;

	QAOA_Graph() {
		for (size_t i = 0; i < n_vertices; ++i)
			for (size_t j = 0; j < n_vertices; ++j) {
				connection[i][j] = -1;
				weight[i][j] = 0;
			}
	}
	int connection[n_vertices][n_vertices]; // -1 not connected, 1 connected
	double weight[n_vertices][n_vertices];

	QAOA_Graph(const QAOA_Graph<n_vertices, qtraits> &g2) {
		for (size_t i = 0; i < n_vertices; ++i)
			for (size_t j = 0; j < n_vertices; ++j) {
				connection[i][j] = g2.connection[i][j];
				weight[i][j] = g2.weight[i][j];
			}
	}

	void add_edge(size_t i, size_t j, double weight_ = 0) {
		if (i >= n_vertices || j >= n_vertices) throw std::runtime_error("Bad index!");

		connection[i][j] = 1;
		connection[j][i] = 1;
		weight[i][j] = weight_;
		weight[j][i] = weight_;
	}

	std::string to_string() const {
		std::stringstream ss;
		for (size_t i = 0; i < n_vertices; ++i)
			for (size_t j = i + 1; j < n_vertices; ++j) {
				if (connection[i][j] == 1)
					ss << "(" << i << "," << j << "):" << weight[i][j] << std::endl;
			}
		return ss.str();
	}

	std::string to_string_connection() const {
		std::stringstream ss;
		for (size_t i = 0; i < n_vertices; ++i) {
			ss << i << " : ";
			for (size_t j = 0; j < n_vertices; ++j) {
				if (connection[i][j] == 1) {
					ss << j << " ";
				}
			}
			ss << std::endl;
		}
		return ss.str();
	}

	void refresh_weight_normal(double mean = 0.0, double errors = 1.0) {
		std::default_random_engine reng;
		std::normal_distribution<double> nd(mean, errors);
		for (size_t i = 0; i < n_vertices; ++i) {
			for (size_t j = i; j < n_vertices; ++j) {
				if (connection[i][j] == 1) {
					auto weight_ = nd(reng);
					weight[i][j] = weight_;
					weight[j][i] = weight_;
				}
			}
		}
	}

	void refresh_weight_uniform(double min = -1.0, double max = 1.0) {
		std::default_random_engine reng;
		std::uniform_real_distribution<double> nd(min, max);
		for (size_t i = 0; i < n_vertices; ++i) {
			for (size_t j = i; j < n_vertices; ++j) {
				if (connection[i][j] == 1) {
					auto weight_ = nd(reng);
					weight[i][j] = weight_;
					weight[j][i] = weight_;
				}
			}
		}
	}

	size_t get_n_edges() const {
		size_t n_edges = 0;
		for (size_t i = 0; i < n_vertices; ++i) {
			for (size_t j = i + 1; j < n_vertices; ++j) {
				if (connection[i][j] == 1) {
					n_edges++;
				}
			}
		}
		return n_edges;
	}

	std::vector<double> obtain_diag_hamiltonian() const {
		std::vector<double> s((1ull << n_vertices), 0);
		for (size_t p = 0; p < (1ull << n_vertices); ++p) {
			for (size_t i = 0; i < n_vertices; ++i) {
				for (size_t j = i; j < n_vertices; ++j) {
					if (connection[i][j] == 1) {
						auto weight_ = weight[i][j];
						if (((p >> i) % 2 + (p >> j) % 2) % 2) {
							s[p] -= weight_;
						}
						else
							s[p] += weight_;
					}
				}
			}
		}
		return s;
	}

	size_t get_max_cut() const {
		double max_cut = 0;
		size_t cut_cnf = 0;
		for (size_t i = 0; i < (pow2(n_vertices)); ++i) {
			double cut = get_cut_value(i);
			if (cut >= max_cut) {
				cut_cnf = i;
				max_cut = cut;
			}
		}
		return cut_cnf;
	}

	double get_cut_value(int cut_config) const {
		double cut_value = 0;
		for (size_t i = 0; i < n_vertices; ++i) {
			for (size_t j = i; j < n_vertices; ++j) {
				if (connection[i][j] == 1) {
					auto weight_ = weight[i][j];
					if ((cut_config >> i) % 2
						!=
						(cut_config >> j) % 2
						) {
						cut_value += weight_;
					}
				}
			}
		}
		return cut_value;
	}
};

template<typename QAOA_Graph, size_t n_step>
struct QAOA_Problem {
	using qtraits = typename QAOA_Graph::qtraits;
	using fp_t = typename qtraits::value_t;
	using qid = typename qtraits::qidx_t;
	static constexpr size_t n_vertices = QAOA_Graph::n_vertices;

	std::unique_ptr<QAOA_Graph> graph = nullptr;
	std::unique_ptr<ObjectiveFunction<qtraits>> obj = nullptr;

	QAOA_Problem(const QAOA_Graph &g) :graph(g) {	}

	
};

ns_end