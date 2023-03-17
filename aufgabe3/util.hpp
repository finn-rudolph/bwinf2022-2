#ifndef UTIL_H
#define UTIL_H 1

#include <bits/stdc++.h>
using namespace std;

size_t factorial(size_t n);

size_t ind(vector<unsigned> const &p);

vector<unsigned> gamma(vector<unsigned> const &p, size_t i);

size_t ind_gamma(vector<unsigned> const &p, size_t i);

vector<unsigned> ith_permutation(size_t n, size_t i);

vector<unsigned> reverse_and_eat(vector<unsigned> const &p, size_t i);

#endif