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
    std::vector<std::pair<int, int>> edges = g.edge_list();
    std::unordered_set<uint64_t> present;
    present.reserve(edges.size() * 2);
    for (auto [u, v] : edges) present.insert(edge_key(u, v));

    std::uniform_int_distribution<size_t> pick(0, edges.size() - 1);
    std::bernoulli_distribution coin(0.5);
    long long tries = (long long)Q * g.E;
    for (long long k = 0; k < tries; ++k) {
        size_t a = pick(rng), b = pick(rng);
        auto [u, v] = edges[a];
        auto [s, t] = edges[b];
        if (coin(rng)) std::swap(s, t);

        // u~v, s~t  ->  u~t, s~v
        if (u == t || s == v) continue;
        if (present.count(edge_key(u, t)) || present.count(edge_key(s, v)))
            continue;

        present.erase(edge_key(u, v));
        present.erase(edge_key(s, t));
        present.insert(edge_key(u, t));
        present.insert(edge_key(s, v));
        edges[a] = {u, t};
        edges[b] = {s, v};
    }

    Graph h(g.N);
    for (auto [u, v] : edges) h.add_edge(u, v);
    return h;
}
