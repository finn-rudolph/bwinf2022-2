#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

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

void ith_permutation(size_t n, size_t i, vector<unsigned> &p)
{
    // Berechne die Ziffern von i im fakultätsbasierten Zahlensystem.
    for (size_t j = 1; j <= n; j++)
    {
        p[n - j] = i % j;
        i /= j;
    }

    size_t const lgn = 32 - __builtin_clz(p.size()), m = 1 << lgn;
    unsigned tree[2 * m];
    for (size_t l = 0; l <= lgn; l++)
        for (size_t j = 0; j < (1ULL << l); j++)
            tree[(1 << l) + j] = 1 << (lgn - l);

    for (size_t j = 0; j < n; j++)
    {
        size_t z = 1;
        for (size_t l = 0; l < lgn; l++)
        {
            tree[z]--;
            z <<= 1;
            if (p[j] >= tree[z])
                p[j] -= tree[z++];
        }
        tree[z] = 0;
        p[j] = z - m;
    }
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
