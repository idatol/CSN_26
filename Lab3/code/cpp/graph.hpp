#pragma once
#include <string>
#include <utility>
#include <vector>

// Undirected simple graph stored as a list.
struct Graph {
    int N = 0;
    long long E = 0;
    std::vector<std::vector<int>> adj;
    std::vector<std::string> names;  // names[i] = word of vertex i (may be empty)

    explicit Graph(int n = 0) : N(n), adj(n), names(n) {}

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
        ++E;
    }
    int degree(int u) const { return (int)adj[u].size(); }

    // Every edge exactly once as (u, v) with u < v.
    std::vector<std::pair<int, int>> edge_list() const;
};

// Reads a file written by data_preparation.py
// Words that do not appear in any edge become isolated vertices.
Graph load_graph(const std::string& path);
