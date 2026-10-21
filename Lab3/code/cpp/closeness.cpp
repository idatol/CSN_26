#include "closeness.hpp"

#include <algorithm>
#include <numeric>

double closeness_vertex(const Graph& g, int i, BFS& bfs) {
    if (g.degree(i) == 0 || g.N < 2) return 0.0;

    std::fill(bfs.dist.begin(), bfs.dist.end(), -1);
    int head = 0, tail = 0;
    bfs.queue[tail++] = i;
    bfs.dist[i] = 0;
    double sum = 0.0;
    while (head < tail) {
        int u = bfs.queue[head++];
        int du = bfs.dist[u] + 1;
        for (int v : g.adj[u]) {
            if (bfs.dist[v] != -1) continue;
            bfs.dist[v] = du;
            sum += 1.0 / du;
            bfs.queue[tail++] = v;
        }
    }
    return sum / (g.N - 1);
}

double mean_closeness(const Graph& g) {
    BFS bfs(g.N);
    double sum = 0.0;
    for (int i = 0; i < g.N; ++i) sum += closeness_vertex(g, i, bfs);
    return sum / g.N;
}

double mean_closeness_estimate(const Graph& g, int M, std::mt19937_64& rng) {
    std::vector<int> perm(g.N);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    M = std::min(M, g.N);
    BFS bfs(g.N);
    double sum = 0.0;
    for (int k = 0; k < M; ++k) sum += closeness_vertex(g, perm[k], bfs);
    return sum / M;
}

bool closeness_at_least(const Graph& g, double x, const std::vector<int>& order,
                        int* visited) {
    BFS bfs(g.N);
    double sum = 0.0;  // sum of C_i over the first M vertices of `order`
    const double N = g.N;
    for (int M = 1; M <= g.N; ++M) {
        sum += closeness_vertex(g, order[M - 1], bfs);
        double cmin = sum / N;
        double cmax = sum / N + 1.0 - M / N;
        if (cmin >= x || cmax < x) {
            if (visited) *visited = M;
            return cmin >= x;
        }
    }
    if (visited) *visited = g.N;
    return sum / N >= x;  // not reached in practice: at M = N, cmin == cmax
}
