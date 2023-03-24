#ifndef UTIL_H
#define UTIL_H 1

#include <bits/stdc++.h>
using namespace std;

uint64_t factorial(uint64_t n);

uint64_t ind(vector<unsigned> const &p);

vector<unsigned> gamma(vector<unsigned> const &p, unsigned i);

uint64_t ind_gamma(vector<unsigned> const &p, unsigned i);

vector<unsigned> ith_permutation(unsigned n, uint64_t i);

vector<unsigned> reverse_and_eat(vector<unsigned> const &p, unsigned i);

#endif