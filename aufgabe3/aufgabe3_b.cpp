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
        vector<pair<size_t, unsigned>> q;

        for (unsigned a = 0; a < 3; a++)
            for (size_t l : y[a])
            {
                ith_permutation(k - 1, l, p);
                for (size_t i = 0; i < k; i++)
                    for (size_t r = 0; r < k; r++)
                        q.emplace_back(make_pair(ind_gamma_inv(p, i, r), a));
            }

        cout << q.size() << '\n';
        sort(q.begin(), q.end());
        for (size_t i = 0; i < q.size(); i++)
        {
            if (i + k - 1 < q.size() && q[i + k - 1].first == q[i].first)
                z[q[i + k - 1].second].push_back(q[i].first);
        }

        swap(y, z);

        if (!y[2][0])
            y[2].erase(y[2].begin());

        if (y[0].empty())
        {
            swap(y[0], y[1]);
            swap(y[1], y[2]);
        }
        else
            z_steps++;

        if (k == 4)
            y[2].push_back(0);

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