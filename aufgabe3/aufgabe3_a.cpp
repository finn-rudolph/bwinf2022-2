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

pair<vector<unsigned>, bool> shortest_ops_dfs_r(
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

    vector<unsigned> res;
    bool found_better = 0;

    while (!q.empty() && q.top().first < b)
    {
        auto [ops, found] =
            shortest_ops_dfs_r(gamma(p, q.top().second), vis, b - 1);

        if (found && ops.size() + 1 < b)
        {
            b = ops.size() + 1;
            res = ops;
            res.push_back(q.top().second);
            found_better = 1;
        }

        q.pop();
    }

    vis[p.size() - 1].insert(ind(p));
    return {res, found_better};
}

vector<unsigned> shortest_ops_dfs(vector<unsigned> const &p)
{
    vector<unordered_set<size_t>> vis(p.size());
    vector<unsigned> res = shortest_ops_dfs_r(p, vis).first;
    reverse(res.begin(), res.end());
    return res;
}

struct node
{
    size_t n, i;
    unsigned b;
    vector<unsigned> ops;

    bool operator<(node const &x) const
    {
        if (b == x.b)
            return n > x.n;
        return b > x.b;
    }
};

vector<unsigned> shortest_ops_bfs(vector<unsigned> const &p)
{
    vector<unordered_set<size_t>> vis(p.size());
    vis[p.size() - 1].insert(ind(p));

    priority_queue<node> q;
    q.push({p.size(), ind(p), lbound(p), {}});

    unsigned b = p.size();
    vector<unsigned> s(p.size());
    vector<unsigned> res;

    while (!q.empty() && q.top().b < b)
    {
        node const x = q.top();
        q.pop();

        if (!x.i)
        {
            if (x.ops.size() < b)
            {
                res = x.ops;
                b = res.size();
            }
            continue;
        }

        s.resize(x.n);
        ith_permutation(x.n, x.i, s);

        for (size_t i = 0; i < x.n; i++)
        {
            unsigned const l = lbound_gamma(s, i);
            if (l + x.ops.size() < b)
            {
                node y = {x.n - 1, ind_gamma(s, i), l + 1, x.ops};
                y.ops.push_back(i);

                if (vis[y.n - 1].find(y.i) == vis[y.n - 1].end())
                {
                    q.push(y);
                    vis[y.n - 1].insert(y.i);
                }
            }
        }
    }

    return res;
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
    vector<unsigned> ops = shortest_ops_bfs(p);

    if (!ops.empty())
    {
        cout << ops.size() << " Operationen nötig. Dafür hinter folgenden"
                              " Indizes wenden:\n\n";
        cout << "Index |  p\n";

        for (auto it = ops.cbegin(); it != ops.cend(); it++)
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