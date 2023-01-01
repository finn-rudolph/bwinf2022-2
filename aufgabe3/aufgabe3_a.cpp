#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

bool is_increasing(bool incr, unsigned a, unsigned b, unsigned &num_monotone)
{
    if (incr ^ (a < b))
    {
        num_monotone++;
        return !incr;
    }
    return incr;
}

unsigned lbound(vector<unsigned> const &p)
{
    if (p.size() == 1)
        return 0;

    bool incr = p[1] > p[0];
    unsigned num_monotone = 1;

    for (size_t i = 1; i < p.size(); i++)
        incr = is_increasing(incr, p[i - 1], p[i], num_monotone);

    return (num_monotone + 1) / 3;
}

unsigned lbound_gamma(vector<unsigned> const &p, size_t i)
{
    assert(i < p.size());

    if (p.size() <= 2)
        return 0;

    bool incr =
        i >= 2 ? (p[i - 2] > p[i - 1])
               : (i == 1 ? (p[i + 1] > p[i - 1]) : (p[i + 2] > p[i + 1]));
    unsigned num_monotone = 0;

    for (size_t j = 2; j < i; j++)
        incr = is_increasing(incr, p[i - j], p[i - j - 1], num_monotone);

    if (i && i < p.size())
        incr = is_increasing(incr, p[0], p[i + 1], num_monotone);

    for (size_t j = i + 2; j < p.size(); j++)
        incr = is_increasing(incr, p[j - 1], p[j], num_monotone);

    return (num_monotone + 1) / 3;
}

vector<unsigned> shortest_ops(vector<unsigned> const &p)
{
    vector<unsigned> bounds(p.size());
}

int main()
{
    size_t n;
    cin >> n;

    vector<unsigned> p(n);
    for (unsigned &x : p)
    {
        cin >> x;
        x--;
    }
}