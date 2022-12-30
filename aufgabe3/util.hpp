#ifndef UTIL_HPP
#define UTIL_HPP 1

#include <bits/stdc++.h>
using namespace std;

void precalc_factorial(size_t n);

size_t ind(vector<unsigned> const &p);

size_t ind_gamma(vector<unsigned> const &p, unsigned i);

vector<unsigned> gamma_inv(vector<unsigned> const &p, unsigned i, unsigned r);

size_t ind_gamma_inv(vector<unsigned> const &p, unsigned i, unsigned j);

void ith_permutation(size_t n, size_t i, vector<unsigned> &p);

vector<unsigned> ith_permutation(size_t n, size_t i);

#endif