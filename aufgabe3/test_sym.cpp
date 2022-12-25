#include <bits/stdc++.h>
using namespace std;

vector<size_t> factorial;

void calc_factorial(unsigned n)
{
    factorial = vector<size_t>(n + 1);
    factorial[0] = 1;
    for (size_t i = 1; i <= n; i++)
        factorial[i] = factorial[i - 1] * i;
}

// Gibt den Index von p in einer lexikographisch sortierten Liste aller
// Permutationen der Länge |p| zurück.
size_t index(vector<unsigned> p)
{
    size_t x = 0;
    for (size_t i = 0; i < p.size(); i++)
    {
        x += p[i] * factorial[p.size() - i - 1];
        for (size_t j = i + 1; j < p.size(); j++)
            if (p[j] > p[i])
                p[j]--;
    }
    return x;
}

void test_symmetry(unsigned n)
{
    vector<unsigned> p(n);
    for (size_t i = 0; i < n; i++)
        p[i] = i;
    vector<vector<unsigned>> res(factorial[n]);
    size_t k = 0;

    do
    {
        unsigned min_ops = n;
        for (size_t i = 1; i <= n; i++)
        {
            reverse(p.begin(), p.begin() + i);
            for (size_t j = 1; j < n; j++)
                if (p[j] > p[0])
                    p[j]--;

            res[k].push_back(index(vector<unsigned>(p.begin() + 1, p.end())));

            for (size_t j = 1; j < n; j++)
                if (p[j] >= p[0])
                    p[j]++;
            reverse(p.begin(), p.begin() + i);
        }
        k++;

    } while (next_permutation(p.begin(), p.end()));

    for (size_t i = 0; i < factorial[n]; i++)
        for (size_t j = 0; j < n; j++)
            assert(res[i][j] == factorial[n - 1] - res[factorial[n] - i - 1][j] - 1);
}

int main()
{
    size_t n;
    cin >> n;
    calc_factorial(n);

    for (size_t i = 2; i <= n; i++)
        test_symmetry(i);
}