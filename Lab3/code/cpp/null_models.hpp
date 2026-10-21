#pragma once
#include <random>

#include "graph.hpp"

// Erdos-Renyi G(N, E): N vertices and exactly E edges chosen uniformly at
// random, with no loops and no multiedges.
Graph erdos_renyi(int N, long long E, std::mt19937_64& rng);

// Switching model: starts from g and tries Q*E random switchings (failed
// attempts count too). Keeps the degree sequence, no loops, no multiedges.
Graph switching_model(const Graph& g, int Q, std::mt19937_64& rng);
