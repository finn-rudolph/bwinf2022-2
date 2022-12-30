#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

vector<__uint128_t> factorial;

void precalc_factorial(size_t n)
{
    factorial = vector<__uint128_t>(n + 1);
    factorial[0] = 1;
    for (size_t i = 1; i <= n; i++)
        factorial[i] = factorial[i - 1] * i;
}

size_t ind(vector<unsigned> const &p)
{
    uint64_t mask = 0; // Bit i ist 1, wenn n - i - 1 bereits aufgetreten ist.
    size_t k = 0;

    for (size_t j = 0; j < p.size(); j++)
    {
        k += (p[j] - __builtin_popcount(mask >> (p.size() - p[j] - 1))) *
             factorial[p.size() - j - 1];
        mask ^= 1 << (p.size() - p[j] - 1);
    }

    assert(mask == (1ULL << p.size()) - 1);
    return k;
}

size_t ind_gamma(vector<unsigned> const &p, unsigned i)
{
    assert(factorial.size() >= p.size());
    assert(i < p.size());

    uint64_t mask = 0;
    size_t k = 0;

    // Um die Permutation in richtiger Reihenfolge von links nach rechts
    // abzuarbeiten, wird der umgekehrte Teil umgekehrt durchgegangen.
    for (size_t j = i - 1; j < p.size(); j--)
    {
        unsigned const x = p[j] - (p[j] > p[i]);
        k += (x - __builtin_popcount(mask >> (p.size() - x - 2))) *
             factorial[p.size() - (i - j) - 1];
        mask ^= 1 << (p.size() - x - 2);
    }

    for (size_t j = i + 1; j < p.size(); j++)
    {
        unsigned const x = p[j] - (p[j] > p[i]);
        k += (x - __builtin_popcount(mask >> (p.size() - x - 2))) *
             factorial[p.size() - j - 1];
        mask ^= 1 << (p.size() - x - 2);
    }

    return k;
}

vector<unsigned> gamma_inv(vector<unsigned> const &p, unsigned i, unsigned r)
{
    assert(i < p.size() + 1);
    assert(r <= p.size() + 1);

    vector<unsigned> q(p.size() + 1);

    for (size_t j = 0; j < i; j++)
        q[j] = p[i - j - 1] + (p[i - j - 1] >= r);
    q[i] = r;
    for (size_t j = i + 1; j < q.size(); j++)
        q[j] = p[j - 1] + (p[j - 1] >= r);

    return q;
}

size_t ind_gamma_inv(vector<unsigned> const &p, unsigned i, unsigned r)
{
    assert(i < p.size() + 1);
    assert(r <= p.size() + 1);

    size_t k = 0;
    uint64_t mask = 0;

    for (size_t j = 0; j < i; j++)
    {
        unsigned const x = p[i - j - 1] + (p[i - j - 1] >= r);
        k += (x - __builtin_popcount(mask >> (p.size() - x))) *
             factorial[p.size() - j];
        mask ^= 1 << (p.size() - x);
    }

    k += (r - __builtin_popcount(mask >> (p.size() - r))) *
         factorial[p.size() - i];
    mask ^= (p.size() - r);

    for (size_t j = i + 1; j < p.size() + 1; j++)
    {
        unsigned const x = p[j - 1] + (p[j - 1] >= r);
        k += (x - __builtin_popcount(mask >> (p.size() - x))) *
             factorial[p.size() - j];
        mask ^= 1 << (p.size() - x);
    }

    return k;
}

void ith_permutation(size_t n, size_t i, vector<unsigned> &p)
{
    assert(factorial.size() >= n);
    assert(i < factorial[n]);
    assert(p.size() >= n);

    // Berechne die Ziffern von i im fakultätsbasierten Zahlensystem.
    for (size_t j = 0; j < n; j++)
    {
        p[j] = i / factorial[n - j - 1];
        i %= factorial[n - j - 1];
    }

    // Addiere zu jeder Ziffer die Anzahl kleiner oder gleicher, links gelegener
    // Ziffern.
    for (size_t j = n - 1; j; j--)
        for (size_t k = j - 1; k < n; k--)
            p[j] += (p[k] <= p[j]);
}

vector<unsigned> ith_permutation(size_t n, size_t i)
{
    vector<unsigned> p(n);
    ith_permutation(n, i, p);
    return p;
}