#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

struct node
{
    uint8_t a, c;
};

int main()
{
    size_t n, m;
    cin >> n >> m;

    array<vector<uint8_t>, 2> y;
    y[0] = {0};

    for (size_t k = 2; k <= m; k++)
    {
        y[1].resize(y[0].size() * k);
        vector<unsigned> p(k);

        for (size_t i = 0; i < k; i++)
            p[i] = i;
        size_t j = 0;

        do
        {
            y[1][j] = k;
            for (size_t i = 0; i < k; i++)
                y[1][j] = min<uint8_t>(y[1][j], y[0][ind_gamma(p, i)] + 1);
            j++;
        } while (next_permutation(p.begin(), p.end()));

        y[1][0] = 0;
        swap(y[0], y[1]);
    }

    array<unordered_map<size_t, node>, 2> z;
    for (size_t i = 0; i < y[0].size(); i++)
        z[0][i] = {y[0][i], 0};
    y[0].clear();
    y[1].clear();

    for (size_t k = m; k < n; k++)
    {
        for (auto const &[i, x] : z[0])
        {
            vector<unsigned> p = ith_permutation(k, i);
            for (size_t j = 0; j < k; j++)
                for (size_t r = 0; r < k; r++)
                {
                    auto it = z[1].find(ind(gamma_inv(p, j, r)));
                    if (it == z[1].end())
                        z[1][ind(gamma_inv(p, j, r))] = (node){x.a + 1, 1};
                    else
                    {
                        it->second.a = min<uint8_t>(it->second.a, x.a + 1);
                        it->second.c++;
                    }
                }
        }

        for (auto const &[i, x] : z[1])
            if (i && x.c == k)
                z[0][i] = x;
    }

    size_t max_a = 0, example_i = -1;
    for (auto const &[i, x] : z[0])
    {
        if (x.a > max_a)
        {
            max_a = x.a;
            example_i = i;
        }
    }

    cout << "P(" << n << ") = " << max_a << '\n'
         << "Beispiel für eine Permutation p mit A(p) = " << max_a << ": ";
    for (unsigned x : ith_permutation(n, example_i))
        cout << x << ' ';
    cout << '\n';
}