// Lab 3 - significance of the mean closeness centrality.
//
// Usage:
//   ./lab3 info   <network.txt>
//   ./lab3 metric <network.txt> [M/N ratio]           (no ratio = exact)
//   ./lab3 test   <network.txt> <er|switching> <T> [Q] [seed]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <string>

#include "closeness.hpp"
#include "graph.hpp"
#include "null_models.hpp"

using Clock = std::chrono::steady_clock;
static double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

static int usage() {
    std::fprintf(stderr,
                 "usage:\n"
                 "  lab3 info   <network.txt>\n"
                 "  lab3 metric <network.txt> [M/N ratio]\n"
                 "  lab3 test   <network.txt> <er|switching> <T> [Q] [seed]\n");
    return 1;
}

int main(int argc, char** argv) {
    if (argc < 3) return usage();
    std::string cmd = argv[1];
    Graph g = load_graph(argv[2]);

    if (cmd == "info") {
        int isolated = 0;
        for (int i = 0; i < g.N; ++i) isolated += g.degree(i) == 0;
        double k = 2.0 * g.E / g.N;
        double delta = 2.0 * g.E / ((double)g.N * (g.N - 1));
        std::printf("N=%d E=%lld <k>=%.4f delta=%.3e isolated=%d\n", g.N, g.E, k, delta, isolated);

    } else if (cmd == "metric") {
        auto t0 = Clock::now();
        if (argc >= 4) {
            double ratio = std::atof(argv[3]);
            int M = std::max(1, (int)(ratio * g.N));
            std::mt19937_64 rng(1);
            double c = mean_closeness_estimate(g, M, rng);
            std::printf("closeness ~ %.6f  (M=%d, M/N=%.3f, %.2fs)\n", c, M, ratio, seconds_since(t0));
        } else {
            double c = mean_closeness(g);
            std::printf("closeness = %.6f  (exact, %.2fs)\n", c, seconds_since(t0));
        }

    } else if (cmd == "test") {
        if (argc < 5) return usage();
        std::string model = argv[3];
        int T = std::atoi(argv[4]);
        int Q = argc >= 6 ? std::atoi(argv[5]) : 10;
        unsigned long long seed = argc >= 7 ? std::strtoull(argv[6], nullptr, 10) : 1;
        std::mt19937_64 rng(seed);

        auto t0 = Clock::now();
        double x = mean_closeness(g);
        std::printf("real network: closeness = %.6f (%.2fs)\n", x, seconds_since(t0));

        // Order in which vertices are visited when bounding C_NH.
        // TODO: try a cleverer ordering (see section 3.1 of the guide).
        std::vector<int> order(g.N);
        std::iota(order.begin(), order.end(), 0);

        int f = 0;  // number of null graphs with C_NH >= x
        for (int t = 0; t < T; ++t) {
            auto t1 = Clock::now();
            Graph nh = model == "er" ? erdos_renyi(g.N, g.E, rng) : switching_model(g, Q, rng);
            int M;
            bool ge = closeness_at_least(nh, x, order, &M);
            f += ge;
            std::printf("  graph %d/%d: C_NH %s x  (M/N=%.3f, %.2fs)\n", t + 1, T,
                        ge ? ">=" : "< ", (double)M / g.N, seconds_since(t1));
            std::fflush(stdout);
        }
        std::printf("p-value (%s, T=%d%s) ~ %d/%d = %.4f\n", model.c_str(), T,
                    model == "er" ? "" : (", Q=" + std::to_string(Q)).c_str(), f, T, (double)f / T);

    } else {
        return usage();
    }
    return 0;
}
