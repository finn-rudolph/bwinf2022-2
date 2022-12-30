#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

unsigned perm_type(
    vector<unsigned> const &p, array<vector<size_t>, 3> const &y)
{
    if (!ind(p))
        return 3;

    unsigned type = 0;

    for (size_t i = 0; i < p.size(); i++)
    {
        size_t j = ind_gamma(p, i);

        if (binary_search(y[0].begin(), y[0].end(), j))
            type = max(type, 0U);
        else if (binary_search(y[1].begin(), y[1].end(), j))
            type = max(type, 1U);
        else if (binary_search(y[2].begin(), y[2].end(), j))
            type = max(type, 2U);
        else
            type = 3;
    }

    return type;
}

int main()
{
    size_t n;
    cin >> n;

    if (n == 1)
    {
        cout << "P(1) = 0\n"
             << "Beispiel für eine Permutation p mit A(p) = 0: 1\n";
        return 0;
    }
    if (n == 2)
    {
        cout << "P(2) = 1\n"
             << "Beispiel für eine Permutation p mit A(p) = 1: 2 1\n";
        return 0;
    }

    precalc_factorial(n);
    array<vector<size_t>, 3> y, z;
    y[0] = {3};
    y[1] = {1, 2, 4, 5};
    y[2] = {0};
    unsigned y_steps = 2, z_steps = 2;

    for (size_t k = 4; k <= n; k++)
    {
        vector<unsigned> p(k - 1);

        for (size_t l : y[0])
        {
            ith_permutation(k - 1, l, p);
            for (size_t i = 0; i < k; i++)
                for (size_t r = 0; r < k; r++)
                {
                    unsigned const t = perm_type(gamma_inv(p, i, r), y);
                    if (t <= 2)
                        z[t].push_back(ind_gamma_inv(p, i, r));
                    if (!t)
                        z_steps = y_steps + 1;
                }
        }

        for (size_t l : y[1])
        {
            ith_permutation(k - 1, l, p);
            for (size_t i = 0; i < k; i++)
                for (size_t r = 0; r < k; r++)
                {
                    unsigned const t = perm_type(gamma_inv(p, i, r), y);
                    if (t <= 2)
                        z[t].push_back(ind_gamma_inv(p, i, r));
                    if (!t)
                        z_steps = y_steps + 1;
                }
        }

        swap(y, z);

        if (y[0].empty())
        {
            swap(y[0], y[1]);
            swap(y[1], y[2]);
        }

        sort(y[0].begin(), y[0].end());
        sort(y[1].begin(), y[1].end());
        sort(y[2].begin(), y[2].end());

        y[0].resize(unique(y[0].begin(), y[0].end()) - y[0].begin());
        y[1].resize(unique(y[1].begin(), y[1].end()) - y[1].begin());
        y[2].resize(unique(y[2].begin(), y[2].end()) - y[2].begin());

        z[0].clear();
        z[1].clear();
        z[2].clear();

        y_steps = z_steps;
    }

    cout << "P(" << n << ") = " << y_steps << '\n'
         << "Beispiel für eine Permutation p mit A(p) = " << y_steps << ": ";
    for (unsigned x : ith_permutation(n, y[0][0]))
        cout << x + 1 << ' ';
    cout << '\n';
}