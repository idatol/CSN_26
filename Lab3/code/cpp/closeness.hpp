#pragma once
#include <random>
#include <vector>

#include "graph.hpp"

// Reusable BFS buffers, so we do not allocate N ints for every vertex.
struct BFS {
    std::vector<int> dist, queue;
    explicit BFS(int n) : dist(n, -1), queue(n) {}
};

// C_i = 1/(N-1) * sum_{j != i} 1/d_ij  (unreachable j contribute 0).
double closeness_vertex(const Graph& g, int i, BFS& bfs);

// Exact mean closeness, Eq. 3. O(N * (N + E)).
double mean_closeness(const Graph& g);

// Monte Carlo estimate, Eq. 8: mean of C_i over M uniformly random vertices.
double mean_closeness_estimate(const Graph& g, int M, std::mt19937_64& rng);

// Decides whether mean closeness of g is >= x, visiting vertices in `order`
// and stopping as soon as the bounds of Eqs. 5-6 settle the answer.
// `visited` (optional) returns how many vertices were needed (the M of Eq. 4).
bool closeness_at_least(const Graph& g, double x, const std::vector<int>& order,
                        int* visited = nullptr);
