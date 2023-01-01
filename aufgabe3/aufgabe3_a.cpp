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

    for (size_t i = 2; i < p.size(); i++)
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
    unsigned num_monotone = 1;

    for (size_t j = 2; j < i; j++)
        incr = is_increasing(incr, p[i - j], p[i - j - 1], num_monotone);

    if (i && i + 1 < p.size())
        incr = is_increasing(incr, p[0], p[i + 1], num_monotone);

    for (size_t j = i + 2; j < p.size(); j++)
        incr = is_increasing(incr, p[j - 1], p[j], num_monotone);

    return (num_monotone + 1) / 3;
}

pair<vector<unsigned>, bool> shortest_ops(
    vector<unsigned> const &p, vector<unordered_set<size_t>> &vis,
    unsigned b = UINT_MAX)
{
    if (!ind(p))
        return {{}, 1};

    if (vis[p.size() - 1].find(ind(p)) != vis[p.size() - 1].end())
        return {{}, 0};

    priority_queue<
        pair<unsigned, size_t>, vector<pair<unsigned, size_t>>,
        greater<pair<unsigned, size_t>>>
        q;
    for (size_t i = 0; i < p.size(); i++)
        q.push({lbound_gamma(p, i), i});

    vector<unsigned> ops;
    bool found_better = 0;

    while (!q.empty() && q.top().first < b)
    {
        auto [res, found] = shortest_ops(gamma(p, q.top().second), vis, b - 1);

        if (found && res.size() + 1 < b)
        {
            b = res.size() + 1;
            ops = res;
            ops.push_back(q.top().second);
            found_better = 1;
        }

        q.pop();
    }

    vis[p.size() - 1].insert(ind(p));
    return {ops, found_better};
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

    precalc_factorial(n);
    vector<unordered_set<size_t>> vis(n);
    vector<unsigned> ops = shortest_ops(p, vis).first;

    if (!ops.empty())
    {
        cout << ops.size() << " Operationen nötig. Dafür hinter folgenden"
                              " Indizes wenden:\n\n";
        cout << "Index |  p\n";

        for (auto it = ops.rbegin(); it != ops.rend(); it++)
        {
            cout << left << setw(6) << *it << "|  ";
            for (unsigned const &x : p)
                cout << left << setw(4) << x + 1;
            p = rev_and_eat(p, *it);
            cout << '\n';
        }

        cout << "      |  ";
        for (unsigned const &x : p)
            cout << left << setw(4) << x + 1;
    }
    else
        cout << "Der Stapel ist bereits sortiert.";
    cout << '\n';
}