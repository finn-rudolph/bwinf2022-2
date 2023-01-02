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
    size_t k = 0;
    size_t const lgn = 32 - __builtin_clz(p.size()), m = 1 << lgn;
    unsigned tree[2 * m];
    memset(tree, 0, 2 * m * sizeof(unsigned));

    for (size_t j = 0; j < p.size(); j++)
    {
        size_t z = m + p[j];
        k *= (p.size() - j);
        k += p[j];
        for (size_t l = 0; l < lgn; l++)
        {
            if (z & 1)
                k -= tree[z - 1];
            tree[z]++;
            z >>= 1;
        }
        tree[z]++;
    }

    return k;
}

vector<unsigned> gamma(vector<unsigned> const &p, size_t i)
{
    vector<unsigned> q(p.size() - 1);

    for (size_t j = 0; j < i; j++)
        q[j] = p[i - j - 1] - (p[i - j - 1] > p[i]);
    for (size_t j = i; j < p.size() - 1; j++)
        q[j] = p[j + 1] - (p[j + 1] > p[i]);

    return q;
}

size_t ind_gamma(vector<unsigned> const &p, size_t i)
{
    size_t k = 0;
    size_t const lgn = 32 - __builtin_clz(p.size() - 1), m = 1 << lgn;
    unsigned tree[2 * m];
    memset(tree, 0, 2 * m * sizeof(unsigned));

    for (size_t j = 0; j < i; j++)
    {
        unsigned const x = p[i - j - 1] - (p[i - j - 1] > p[i]);
        size_t z = m + x;
        k *= (p.size() - 1 - j);
        k += x;
        for (size_t l = 0; l < lgn; l++)
        {
            if (z & 1)
                k -= tree[z - 1];
            tree[z]++;
            z >>= 1;
        }
        tree[z]++;
    }

    for (size_t j = i; j < p.size() - 1; j++)
    {
        unsigned const x = p[j + 1] - (p[j + 1] > p[i]);
        size_t z = m + x;
        k *= (p.size() - 1 - j);
        k += x;
        for (size_t l = 0; l < lgn; l++)
        {
            if (z & 1)
                k -= tree[z - 1];
            tree[z]++;
            z >>= 1;
        }
        tree[z]++;
    }

    return k;
}

vector<unsigned> gamma_inv(vector<unsigned> const &p, size_t i, unsigned r)
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

size_t ind_gamma_inv(vector<unsigned> const &p, size_t i, unsigned r)
{
    assert(factorial.size() >= p.size());
    assert(i < p.size() + 1);
    assert(r <= p.size() + 1);

    size_t k = 0;
    uint64_t mask = 0;

    for (size_t j = 0; j < i; j++)
    {
        unsigned const x = p[i - j - 1] + (p[i - j - 1] >= r);
        k += (x - __builtin_popcount(mask >> (p.size() - x))) *
             factorial[p.size() - j];
        mask ^= 1ULL << (p.size() - x);
    }

    k += (r - __builtin_popcount(mask >> (p.size() - r))) *
         factorial[p.size() - i];
    mask ^= 1ULL << (p.size() - r);

    for (size_t j = i + 1; j < p.size() + 1; j++)
    {
        unsigned const x = p[j - 1] + (p[j - 1] >= r);
        k += (x - __builtin_popcount(mask >> (p.size() - x))) *
             factorial[p.size() - j];
        mask ^= 1ULL << (p.size() - x);
    }

    assert(mask == (1ULL << (p.size() + 1)) - 1);
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

vector<unsigned> rev_and_eat(vector<unsigned> const &p, size_t i)
{
    vector<unsigned> q(p.size() - 1);

    copy(p.begin(), p.begin() + i, q.begin());
    reverse(q.begin(), q.begin() + i);
    copy(p.begin() + i + 1, p.end(), q.begin() + i);

    return q;
}
