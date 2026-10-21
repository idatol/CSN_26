#include "null_models.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <unordered_set>

// Key for an undirected edge, independent of the order of u and v.
static inline uint64_t edge_key(int u, int v) {
    if (u > v) std::swap(u, v);
    return (uint64_t)u << 32 | (uint32_t)v;
}

Graph erdos_renyi(int N, long long E, std::mt19937_64& rng) {
    std::uniform_int_distribution<int> pick(0, N - 1);
    std::unordered_set<uint64_t> present;
    present.reserve(E * 2);
    Graph g(N);
    while (g.E < E) {
        int u = pick(rng), v = pick(rng);
        if (u == v) continue;                             // loop
        if (!present.insert(edge_key(u, v)).second) continue;  // multiedge
        g.add_edge(u, v);
    }
    return g;
}

Graph switching_model(const Graph& g, int Q, std::mt19937_64& rng) {
    // TODO: implement. Suggested outline:
    //   edges   = g.edge_list();                  // vector of (u, v)
    //   present = set of edge_key(u, v)           // O(1) multiedge check
    //   repeat Q * E times:
    //     pick two random edges u~v and s~t from `edges`
    //     propose the switch u~t, s~v (also consider the other orientation)
    //     reject it if it creates a loop or a multiedge (see report questions)
    //     otherwise update `edges` and `present`
    //   build and return a Graph from `edges`
    (void)g; (void)Q; (void)rng;
    std::fprintf(stderr, "switching_model is not implemented yet\n");
    std::exit(1);
}
