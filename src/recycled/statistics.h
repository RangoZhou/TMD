#pragma once

#include <algorithm> // sort
#include <cmath> // sqrt

namespace tmd {

// simple but inefficient implementations

inline double mean(const std::vector<double>& v) {
	double acc = 0;
	for(std::size_t i = 0; i < v.size(); ++i)
		acc += v[i];
	return v.empty() ? 0 : (acc/v.size());
}

inline double deviation(const std::vector<double>& v) {
	double m = mean(v);
	double acc = 0;
	for(std::size_t i = 0; i < v.size(); ++i)
		acc += (v[i] - m)*(v[i] - m);
	return v.empty() ? 0 : std::sqrt(acc/v.size());
}

inline double rmsd(const std::vector<double>& a, const std::vector<double>& b) {
	assert(a.size() == b.size());
	double acc = 0;
	for(std::size_t i = 0; i < a.size(); ++i)
		acc += (a[i] - b[i])*(a[i] - b[i]);
	return a.empty() ? 0 : std::sqrt(acc/a.size());
}

inline double average_difference(const std::vector<double>& b, const std::vector<double>& a) { // b - a
	assert(a.size() == b.size());
	double acc = 0;
	for(std::size_t i = 0; i < a.size(); ++i)
		acc += b[i] - a[i];
	return a.empty() ? 0 : (acc/a.size());
}

inline double pearson(const std::vector<double>& x, const std::vector<double>& y) {
	std::size_t n = x.size();
	assert(n == y.size());
	if(n == 0) return 0; // ?
	double sum_x = 0;
	double sum_y = 0;
	double sum_x_sq = 0;
	double sum_y_sq = 0;
	double sum_prod = 0;
	for(std::size_t i = 0; i < n; ++i) {
		sum_x += x[i];
		sum_y += y[i];
		sum_x_sq += x[i]*x[i];
		sum_y_sq += y[i]*y[i];
		sum_prod += x[i] * y[i];
	}
	double sd_x = std::sqrt(sum_x_sq/n - (sum_x/n)*(sum_x/n)); // FIXME the argument is supposed to be >= 0, but ...
	double sd_y = std::sqrt(sum_y_sq/n - (sum_y/n)*(sum_y/n));
	double cov = sum_prod/n - (sum_x/n) * (sum_y/n);
	double tmp = sd_x * sd_y;
	if(std::abs(tmp) < epsilon_fl) return 0; // ?
	return cov / tmp;
}

struct spearman_aux {
	double x;
	std::size_t i;
	spearman_aux(double x, std::size_t i) : x(x), i(i) {}
};

inline bool operator<(const spearman_aux& a, const spearman_aux& b) {
	return a.x < b.x;
}

inline std::vector<double> get_rankings(const std::vector<double>& x) {
	std::vector<spearman_aux> to_sort;
	for(std::size_t i = 0; i < x.size(); ++i)
		to_sort.push_back(spearman_aux(x[i], i));
	std::sort(to_sort.begin(), to_sort.end());

	std::vector<double> tmp(x.size(), 0);
	for(std::size_t rank = 0; rank < to_sort.size(); ++rank)
		tmp[to_sort[rank].i] = rank;
	return tmp;
}

inline double spearman(const std::vector<double>& x, const std::vector<double>& y) {
	return pearson(get_rankings(x), get_rankings(y));
}

}