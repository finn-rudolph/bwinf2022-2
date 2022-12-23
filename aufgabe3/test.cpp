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

vector<unsigned> pwue(unsigned n, vector<unsigned> &y)
{
    vector<unsigned> z(factorial[n]);
    vector<unsigned> p(n);
    for (size_t i = 0; i < n; i++)
        p[i] = i;
    vector<unsigned> x(n + 1, 0);
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

            min_ops = min(min_ops, 1 + y[index(vector<unsigned>(p.begin() + 1, p.end()))]);

            for (size_t j = 1; j < n; j++)
                if (p[j] >= p[0])
                    p[j]++;
            reverse(p.begin(), p.begin() + i);
        }
        if (k)
            x[min_ops]++;
        z[k] = min_ops;
        k++;

    } while (next_permutation(p.begin(), p.end()));

    z[0] = 0;
    swap(y, z);

    return x;
}

int main()
{
    size_t n;
    cin >> n;
    calc_factorial(n);

    vector<unsigned> y(1, 0);

    for (size_t i = 2; i <= n; i++)
    {
        vector<unsigned> x = pwue(i, y);
        for (unsigned &z : x)
            cout << z << ' ';
        cout << endl;
    }
}