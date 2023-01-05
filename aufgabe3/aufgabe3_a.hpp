#include <bits/stdc++.h>
using namespace std;

vector<unsigned> min_operations_bfs(vector<unsigned> const &p);

vector<unsigned> min_operations_astar(vector<unsigned> const &p);

vector<unsigned> min_operations_bnb(vector<unsigned> const &p);

vector<unsigned> reconstruct_operations(
    vector<unordered_map<size_t, size_t>> const &pre, size_t m, size_t i);

void print_operations(vector<unsigned> p, vector<unsigned> const &op);
