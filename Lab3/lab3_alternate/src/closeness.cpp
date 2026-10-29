// Lab 3 (CSN): is the mean closeness centrality significantly large?
// Null models: Erdos-Renyi G(N,E) and degree-preserving edge switching.
//
// Usage: ./closeness <network.txt> [T=100] [Q=0] [mfrac=1.0] [seed=1]
//   T      number of random graphs per null model
//   Q      switching model tries Q*E switches (0 = max(10, ceil(ln E)))
//   mfrac  fraction of vertices used to estimate C of the random graphs
//          (1.0 = exact, with early stop). The real network's C is always exact.
// Output: one TSV line: language N E <k> delta C p_binomial p_switching
// Replicates and BFS sources run in parallel (OpenMP); results do not depend on
// the number of threads (each replicate has its own RNG stream).
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

typedef pair<int, int> Edge;

static const double EPS = 1e-12;

// ---------- graph loading / building ----------

static uint64_t key(int u, int v, int N) {
    if (u > v) swap(u, v);
    return (uint64_t)u * N + v;
}

// Edge list file: header "N E", then one edge per line: word1<TAB>word2.
// A tab separates the words because some word forms contain spaces.
// The header is checked against the file: E must equal the number of edge lines, and
// N must be at least the number of distinct words (extra vertices are isolated).
// Loops and repeated edges are dropped after the check.
static void fail(const string& path, const string& msg) {
    cerr << path << ": " << msg << "\n";
    exit(1);
}

static int load(const string& path, vector<Edge>& edges) {
    ifstream in(path);
    if (!in) fail(path, "cannot open file");
    string line;
    long hN = -1, hE = -1;
    if (!getline(in, line) || sscanf(line.c_str(), "%ld %ld", &hN, &hE) != 2)
        fail(path, "first line must be \"N E\"");
    unordered_map<string, int> id;
    vector<Edge> raw;
    for (long ln = 2; getline(in, line); ln++) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        size_t tab = line.find('\t');
        if (tab == string::npos || line.find('\t', tab + 1) != string::npos)
            fail(path, "line " + to_string(ln) + ": expected exactly two words separated by a tab");
        int u = id.emplace(line.substr(0, tab), id.size()).first->second;
        int v = id.emplace(line.substr(tab + 1), id.size()).first->second;
        raw.push_back({u, v});
    }
    if ((long)raw.size() != hE)
        fail(path, "header says E=" + to_string(hE) + " but the file has " + to_string(raw.size()) + " edges");
    if ((long)id.size() > hN)
        fail(path, "header says N=" + to_string(hN) + " but the file has " + to_string(id.size()) + " distinct words");
    int N = hN;  // N - id.size() isolated vertices (words that only formed loops)
    unordered_set<uint64_t> seen;
    for (auto& e : raw) {
        if (e.first == e.second) continue;
        if (seen.insert(key(e.first, e.second, N)).second) edges.push_back(e);
    }
    return N;
}

// Compressed sparse row adjacency: neighbours of v are nbr[off[v] .. off[v+1])
struct Graph {
    int N;
    vector<int> off, nbr;
    int deg(int v) const { return off[v + 1] - off[v]; }
};

static Graph buildGraph(int N, const vector<Edge>& edges) {
    Graph g;
    g.N = N;
    g.off.assign(N + 1, 0);
    for (auto& e : edges) { g.off[e.first + 1]++; g.off[e.second + 1]++; }
    for (int i = 0; i < N; i++) g.off[i + 1] += g.off[i];
    g.nbr.resize(2 * edges.size());
    vector<int> pos(g.off.begin(), g.off.end() - 1);
    for (auto& e : edges) {
        g.nbr[pos[e.first]++] = e.second;
        g.nbr[pos[e.second]++] = e.first;
    }
    return g;
}

// ---------- closeness ----------

struct Workspace {
    vector<int> dist, queue;
    Workspace(int N) : dist(N, -1), queue(N) {}
};

// C_i = 1/(N-1) * sum_{j != i} 1/d_ij (BFS; unreachable vertices add 0)
static double closeness(const Graph& g, int s, Workspace& w) {
    if (g.deg(s) == 0) return 0.0;
    int head = 0, tail = 0;
    w.queue[tail++] = s;
    w.dist[s] = 0;
    double sum = 0;
    while (head < tail) {
        int u = w.queue[head++];
        int du = w.dist[u] + 1;
        for (int e = g.off[u]; e < g.off[u + 1]; e++) {
            int v = g.nbr[e];
            if (w.dist[v] < 0) {
                w.dist[v] = du;
                sum += 1.0 / du;
                w.queue[tail++] = v;
            }
        }
    }
    for (int i = 0; i < tail; i++) w.dist[w.queue[i]] = -1;
    return sum / (g.N - 1);
}

// Exact mean closeness, parallel over BFS sources: used for the real network (one graph).
static double meanClosenessExact(const Graph& g) {
    int N = g.N;
    double sum = 0;
#pragma omp parallel
    {
        Workspace w(N);
#pragma omp for schedule(dynamic, 16) reduction(+ : sum)
        for (int i = 0; i < N; i++) sum += closeness(g, i, w);
    }
    return sum / N;
}

// Mean closeness over a random fraction of vertices (frac=1 -> exact value).
// Sequential: used inside the parallel loop over replicates.
static double meanCloseness(const Graph& g, double frac, Workspace& w, mt19937_64& rng) {
    int N = g.N;
    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    int M = N;
    if (frac < 1.0) {
        shuffle(order.begin(), order.end(), rng);
        M = max(1, (int)(frac * N));
    }
    double sum = 0;
    for (int i = 0; i < M; i++) sum += closeness(g, order[i], w);
    return sum / M;
}

// Decides C(g) >= thr, stopping as soon as the bounds settle it (exact).
//   lower bound: uncomputed vertices contribute 0
//   upper bound: vertex of degree k has at most k vertices at distance 1 and
//                the rest at distance >= 2, so C_i <= (k + (N-1-k)/2)/(N-1)
// Isolated vertices have C_i = 0 and are never visited.
static bool atLeast(const Graph& g, double thr, Workspace& w, mt19937_64& rng) {
    int N = g.N;
    vector<int> order;
    vector<double> ub(N, 0.0);
    double remaining = 0;
    for (int i = 0; i < N; i++) {
        int k = g.deg(i);
        if (k == 0) continue;
        ub[i] = (k + (N - 1 - k) / 2.0) / (N - 1);
        remaining += ub[i];
        order.push_back(i);
    }
    shuffle(order.begin(), order.end(), rng);
    double sum = 0;
    if (remaining / N < thr - EPS) return false;
    for (int v : order) {
        remaining -= ub[v];
        sum += closeness(g, v, w);
        if (sum / N >= thr - EPS) return true;
        if ((sum + remaining) / N < thr - EPS) return false;
    }
    return sum / N >= thr - EPS;
}

// ---------- null models ----------

// G(N,E): E distinct non-loop edges chosen uniformly.
static vector<Edge> binomialGraph(int N, size_t E, mt19937_64& rng) {
    uniform_int_distribution<int> pick(0, N - 1);
    unordered_set<uint64_t> used;
    used.reserve(2 * E);
    vector<Edge> edges;
    edges.reserve(E);
    while (edges.size() < E) {
        int u = pick(rng), v = pick(rng);
        if (u == v || !used.insert(key(u, v, N)).second) continue;
        edges.push_back({u, v});
    }
    return edges;
}

// Edge switching: take u~v, s~t and rewire to u~t, s~v (the second edge is
// randomly flipped so both pairings are reachable). Rejected when it would
// create a loop or a multiedge; rejected tries still count towards Q*E.
static vector<Edge> switching(vector<Edge> edges, int N, double Q, mt19937_64& rng) {
    size_t E = edges.size();
    unordered_set<uint64_t> used;
    used.reserve(2 * E);
    for (auto& e : edges) used.insert(key(e.first, e.second, N));
    uniform_int_distribution<size_t> pick(0, E - 1);
    size_t tries = (size_t)(Q * E);
    for (size_t t = 0; t < tries; t++) {
        size_t a = pick(rng), b = pick(rng);
        if (a == b) continue;
        int u = edges[a].first, v = edges[a].second;
        int s = edges[b].first, x = edges[b].second;
        if (rng() & 1) swap(s, x);
        // new edges: u~x and s~v
        if (u == x || s == v) continue;                         // loop
        if (used.count(key(u, x, N)) || used.count(key(s, v, N))) continue;  // multiedge
        used.erase(key(u, v, N));
        used.erase(key(s, x, N));
        used.insert(key(u, x, N));
        used.insert(key(s, v, N));
        edges[a] = {u, x};
        edges[b] = {s, v};
    }
    return edges;
}

// ---------- main ----------

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "usage: " << argv[0] << " network.txt [T] [Q] [mfrac] [seed]\n";
        return 1;
    }
    string path = argv[1];
    int T = argc > 2 ? atoi(argv[2]) : 100;
    double Q = argc > 3 ? atof(argv[3]) : 0;
    double frac = argc > 4 ? atof(argv[4]) : 1.0;
    unsigned long seed = argc > 5 ? atol(argv[5]) : 1;

    vector<Edge> edges;
    int N = load(path, edges);
    size_t E = edges.size();
    if (Q <= 0) Q = max(10.0, ceil(log((double)E)));
    Graph g = buildGraph(N, edges);

    string lang = path.substr(path.find_last_of('/') + 1);
    size_t cut = lang.find("_dependency_network");
    if (cut != string::npos) lang = lang.substr(0, cut);

    double C = meanClosenessExact(g);  // always exact, even when frac < 1
    double k = 2.0 * E / N, delta = 2.0 * E / ((double)N * (N - 1));
    cerr << lang << ": C=" << C << " Q=" << Q << " threads=" << omp_get_max_threads() << "\n";

    // With frac<1 the null values are estimates (no early stop); the real C is exact.
    // Replicate r of model m uses its own RNG stream -> independent of thread count.
    int fBin = 0, fSw = 0;
#pragma omp parallel
    {
        Workspace w(N);
#pragma omp for schedule(dynamic, 1) reduction(+ : fBin, fSw)
        for (int r = 0; r < 2 * T; r++) {
            bool sw = r >= T;
            mt19937_64 rng(seed * 1000003UL + 2 * (r % T) + sw + 1);
            vector<Edge> es = sw ? switching(edges, N, Q, rng) : binomialGraph(N, E, rng);
            Graph h = buildGraph(N, es);
            bool ge = frac < 1.0 ? meanCloseness(h, frac, w, rng) >= C - EPS
                                 : atLeast(h, C, w, rng);
            (sw ? fSw : fBin) += ge;
        }
    }

    printf("language\tN\tE\t<k>\tdelta\tC\tp_binomial\tp_switching\n");
    printf("%s\t%d\t%zu\t%.4f\t%.6f\t%.6f\t%.4f\t%.4f\n", lang.c_str(), N, E, k, delta, C,
           (double)fBin / T, (double)fSw / T);
    return 0;
}
