#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

int main()
{
    size_t n;
    cin >> n;

    uint8_t *y = (uint8_t *)malloc(sizeof(uint8_t)),
            *z = (uint8_t *)malloc(sizeof(uint8_t));
    y[0] = 0;
    size_t u = 1, v = 1;

    for (size_t k = 2; k < n; k++)
    {
        u *= k;
        z = (uint8_t *)realloc(z, u * sizeof(uint8_t));

        vector<unsigned> p(k);
        for (size_t i = 0; i < k; i++)
            p[i] = i;

        for (size_t i = 0; i < u / 2; i++)
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

        z[0] = 0;
        swap(y, z);
        v *= k;
    }

    u *= n;
    unsigned max_a = 0;
    vector<unsigned> p(n);
    for (size_t i = 0; i < n; i++)
        p[i] = i;

    for (size_t i = 0; i < u; i++)
    {
        unsigned a = n;

        for (size_t j = 0; j < n; j++)
            a = min<unsigned>(a, y[ind_gamma(p, j)] + 1);

        max_a = max(max_a, a);
        next_permutation(p.begin(), p.end());
    }

    cout << "P(" << n << ") = " << max_a << '\n';
}