#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

int main()
{
    size_t n, m;
    cin >> n >> m;

    vector<uint64_t> factorial(n + 1);
    factorial[0] = 1;
    for (size_t i = 1; i <= n; i++)
        factorial[i] = factorial[i - 1] * i;

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

    vector<unordered_map<size_t, pair<uint8_t, uint8_t>>> z(n + 1);
    size_t a_max = 0, example = 0;

    for (size_t i_m = 0; i_m < y[0].size(); i_m++)
    {
        queue<pair<size_t, uint8_t>> q;
        q.push({i_m, y[0][i_m]});

        for (size_t k = m; k < n; k++)
        {
            while (!q.empty())
            {
                auto const [i, a] = q.front();
                q.pop();

                vector<unsigned> const p = ith_permutation(k, i);

                for (size_t j = 0; j < k + 1; j++)
                {
                    vector<unsigned> const s = inv_permutation(gamma_inv(p, j, k));
                    size_t mu = ind(gamma_inv(p, j, 0));

                    for (size_t r = 0; r < k + 1; r++)
                    {
                        auto it = z[k + 1].find(mu);
                        if (it == z[k + 1].end())
                            z[k + 1].emplace(mu, make_pair<uint8_t, uint8_t>(a + 1, 1));
                        else
                        {
                            it->second.first = min<uint8_t>(it->second.first, a + 1);
                            it->second.second++;
                        }

                        if (s[r] < j)
                            mu -= factorial[k - s[r]];
                        else if (s[r] > j)
                            mu += factorial[k - j];
                    }
                }
            }

            auto it = z[k + 1].begin();
            while (it != z[k + 1].end())
            {
                if (it->second.second == k + 1)
                {
                    if (it->first)
                        q.emplace(it->first, it->second.first);
                    it = z[k + 1].erase(it);
                }
                else
                    it++;
            }
        }

        if (!q.empty())
        {
            example = q.front().first;
            a_max = q.front().second;
            break;
        }
    }

    cout << "P(" << n << ") = " << a_max << '\n'
         << "Beispiel für eine Permutation p mit A(p) = " << a_max << ": ";
    for (unsigned x : ith_permutation(n, example))
        cout << x + 1 << ' ';
    cout << '\n';
}