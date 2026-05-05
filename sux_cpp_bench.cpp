#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <sux/bits/EliasFano.hpp>

using timer = std::chrono::high_resolution_clock;

static std::vector<uint32_t> read_data_binary(const std::string &path) {
	std::ifstream in(path, std::ios::binary);
	if (!in) {
		std::cerr << "Error opening " << path << std::endl;
		return {};
	}
	uint64_t size;
	in.read(reinterpret_cast<char *>(&size), sizeof(size));
	std::vector<uint32_t> data(size);
	in.read(reinterpret_cast<char *>(data.data()), size * sizeof(uint32_t));
	if (!in) {
		std::cerr << "Error reading " << path << std::endl;
		return {};
	}
	return data;
}

struct BenchResult {
	double build_ns, select_ns, bpk, rank_ns;
};

template <int LOG2_ONES, int LOG2_WORDS>
BenchResult bench_ef(const std::vector<uint64_t> &ones, uint64_t u, size_t n,
					 const std::vector<size_t> &rands1, const std::vector<uint32_t> &rands2) {
	size_t times_test = rands1.size();

	auto build_start = timer::now();
	sux::bits::EliasFano<sux::util::AllocType::MALLOC, LOG2_ONES, LOG2_WORDS> ef(ones, u);
	double build_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - build_start).count() / static_cast<double>(n);

	size_t cnt = 0;
	auto start = timer::now();
	for (size_t i = 0; i < times_test; ++i)
		cnt += ef.select(rands1[i]);
	const volatile size_t _sink1 = cnt;
	double select_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - start).count() / static_cast<double>(times_test);

	double bpk = static_cast<double>(ef.bitCount()) / static_cast<double>(n);

	cnt = 0;
	start = timer::now();
	for (size_t i = 0; i < times_test; ++i)
		cnt += ef.rank(rands2[i]);
	const volatile size_t _sink2 = cnt;
	double rank_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - start).count() / static_cast<double>(times_test);

	return {build_ns, select_ns, bpk, rank_ns};
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: sux_cpp_bench <file1> [file2 ...]" << std::endl;
		return 1;
	}

	std::cout << std::scientific << std::setprecision(3);

	std::cout << "filename,n,u,ratio";
	for (int log2ones : {9, 10, 11}) {
		for (int log2words : {2, 3}) {
			std::string tag = "sux_cpp_ef_" + std::to_string(log2ones) + "_" + std::to_string(log2words);
			std::cout << "," << tag << "_time_build"
					  << "," << tag << "_time_select"
					  << "," << tag << "_bpk"
					  << "," << tag << "_time_rank";
		}
	}
	std::cout << std::endl;

	for (int h = 1; h < argc; ++h) {
		std::string filepath(argv[h]);
		auto data = read_data_binary(filepath);
		if (data.size() < 5) {
			std::cerr << "Skipping " << filepath << ": too few elements (" << data.size() << ")" << std::endl;
			continue;
		}

		size_t n = data.size();
		uint64_t u = static_cast<uint64_t>(data.back()) + 1;
		double ratio = static_cast<double>(n) / static_cast<double>(u);

		std::size_t found = filepath.find_last_of('/');
		std::string filename = filepath.substr(found + 1);

		std::cout << "'" << filename << "'," << n << "," << u << "," << ratio;

		size_t TIMES_TEST = n / 5;

		std::mt19937 mt1(2323);
		std::uniform_int_distribution<size_t> dist1(0, n - 2);
		std::vector<size_t> rands1(TIMES_TEST);
		for (size_t i = 0; i < TIMES_TEST; ++i)
			rands1[i] = dist1(mt1);

		std::mt19937 mt2(4242);
		std::uniform_int_distribution<uint32_t> dist2(data.front(), data.back() - 1);
		std::vector<uint32_t> rands2(TIMES_TEST);
		for (size_t i = 0; i < TIMES_TEST; ++i)
			rands2[i] = dist2(mt2);

		std::vector<uint64_t> ones(data.begin(), data.end());

#define BENCH_AND_PRINT(L2O, L2W) \
		{ \
			auto r = bench_ef<L2O, L2W>(ones, u, n, rands1, rands2); \
			std::cout << "," << r.build_ns << "," << r.select_ns << "," << r.bpk << "," << r.rank_ns; \
		}

		BENCH_AND_PRINT(9, 2)
		BENCH_AND_PRINT(9, 3)
		BENCH_AND_PRINT(10, 2)
		BENCH_AND_PRINT(10, 3)
		BENCH_AND_PRINT(11, 2)
		BENCH_AND_PRINT(11, 3)

#undef BENCH_AND_PRINT

		std::cout << std::endl;
	}

	return 0;
}
