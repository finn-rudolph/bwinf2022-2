#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

void update_z(
    size_t k, size_t u, size_t v, uint8_t const *const y, uint8_t *const z,
    size_t i1, size_t i2)
{
    vector<unsigned> p = ith_permutation(k, i1);

    for (size_t i = i1; i < i2; i++)
    {
        z[i] = k;
        z[u - i - 1] = k;

        for (size_t j = 0; j < k; j++)
        {
            size_t const l = ind_gamma(p, j);
            z[i] = min<uint8_t>(z[i], y[l] + 1);
            z[u - i - 1] = min<uint8_t>(z[u - i - 1], y[v - l - 1] + 1);
        }

        next_permutation(p.begin(), p.end());
    }
}

pair<unsigned, size_t> get_max_a(
    size_t n, uint8_t const *const y, size_t i1, size_t i2)
{
    pair<unsigned, size_t> res = {0, -1};
    vector<unsigned> p = ith_permutation(n, i1);

    for (size_t i = i1; i < i2; i++)
    {
        unsigned a = n;

        for (size_t j = 0; j < n; j++)
            a = min<unsigned>(a, y[ind_gamma(p, j)] + 1);

        res = max(res, {a, i});
        next_permutation(p.begin(), p.end());
    }

    return res;
}

int main()
{
    size_t n;
    cin >> n;

    uint8_t *y = (uint8_t *)malloc(sizeof(uint8_t)),
            *z = (uint8_t *)malloc(sizeof(uint8_t));
    y[0] = 0;
    size_t u = 1, v = 1;
    unsigned const n_threads = thread::hardware_concurrency();

    for (size_t k = 2; k < n; k++)
    {
        u *= k;
        z = (uint8_t *)realloc(z, u * sizeof(uint8_t));

        vector<thread> threads;

        for (unsigned i = 0; i < n_threads; i++)
            threads.emplace_back(update_z, k, u, v, y, z,
                                 (u / 2) * i / n_threads,
                                 (u / 2) * (i + 1) / n_threads);

        for (thread &t : threads)
            t.join();

        z[0] = 0;
        swap(y, z);
        v *= k;
    }

    u *= n;
    pair<unsigned, size_t> res = {0, -1};
    vector<future<pair<unsigned, size_t>>> fut;

    for (unsigned i = 0; i < n_threads; i++)
        fut.emplace_back(
            async(get_max_a, n, y, u * i / n_threads, u * (i + 1) / n_threads));

    for (auto &f : fut)
        res = max(res, f.get());

    cout << "P(" << n << ") = " << res.first << '\n';
    cout << "Beispiel mit A(p) = P(" << n << "): ";
    for (unsigned const x : ith_permutation(n, res.second))
        cout << x + 1 << ' ';
    cout << '\n';
}