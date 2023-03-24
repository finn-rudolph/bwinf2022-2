#include <bits/stdc++.h>
using namespace std;

vector<unsigned> min_operations_bfs(vector<unsigned> const &p);

vector<unsigned> min_operations_astar(vector<unsigned> const &p);

vector<unsigned> min_operations_bnb(vector<unsigned> const &p);

vector<unsigned> reconstruct_operations(
    vector<unordered_map<uint64_t, uint64_t>> const &pre, unsigned length,
    uint64_t index);

void print_operations(vector<unsigned> p, vector<unsigned> const &op);
