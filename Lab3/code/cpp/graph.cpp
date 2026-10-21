#include "graph.hpp"

#include <fstream>
#include <stdexcept>
#include <unordered_map>

std::vector<std::pair<int, int>> Graph::edge_list() const {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(E);
    for (int u = 0; u < N; ++u)
        for (int v : adj[u])
            if (u < v) edges.emplace_back(u, v);
    return edges;
}

Graph load_graph(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);

    int N;
    long long E;
    in >> N >> E;
    std::string line;
    std::getline(in, line);

    Graph g(N);
    std::unordered_map<std::string, int> id;
    id.reserve(N * 2);
    auto get_id = [&](const std::string& w) {
        auto it = id.find(w);
        if (it != id.end()) return it->second;
        int i = (int)id.size();
        if (i >= N) throw std::runtime_error(path + ": more distinct words than N");
        id.emplace(w, i);
        g.names[i] = w;
        return i;
    };

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        size_t tab = line.find('\t');
        if (tab == std::string::npos) throw std::runtime_error(path + ": bad line: " + line);
        g.add_edge(get_id(line.substr(0, tab)), get_id(line.substr(tab + 1)));
    }
    if (g.E != E) throw std::runtime_error(path + ": header says E=" + std::to_string(E) +
                                           " but read " + std::to_string(g.E));
    return g;
}
