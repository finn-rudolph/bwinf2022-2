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

unsigned get_lbound(vector<unsigned> const &p)
{
    if (p.size() == 1)
        return 0;

    bool incr = p[1] > p[0];
    unsigned num_monotone = 1;

    for (size_t i = 2; i < p.size(); i++)
        incr = is_increasing(incr, p[i - 1], p[i], num_monotone);

    return (num_monotone + 1) / 3;
}

unsigned get_lbound_gamma(vector<unsigned> const &p, size_t i)
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

// Bestimmt, ob eine Folge an maximal ubound gamma-Operationen zum Sortieren
// von p existiert (zweites Element). Falls ja, wird diese im ersten Element
// zurückgegeben.
pair<vector<unsigned>, bool> shortest_op_dfs_r(
    vector<unsigned> const &p, vector<unordered_set<size_t>> &vis,
    unsigned ubound = UINT_MAX)
{
    if (!ind(p)) // Identische Permutation erreicht.
        return {{}, 1};

    if (vis[p.size() - 1].find(ind(p)) != vis[p.size() - 1].end())
        return {{}, 0};

    priority_queue<
        pair<unsigned, size_t>, vector<pair<unsigned, size_t>>,
        greater<pair<unsigned, size_t>>>
        q;
    for (size_t i = 0; i < p.size(); i++)
        q.push({get_lbound_gamma(p, i), i});

    vector<unsigned> res;
    bool found_solution = 0;

    // Gehe die durch eine gamma-Operation erreichbaren Permutationen
    // aufsteigend nach unterer Schranke durch und bestimme rekursiv die
    // kürzeste Operationenfolge.
    while (!q.empty() && q.top().first < ubound)
    {
        auto [op, found] =
            shortest_op_dfs_r(gamma(p, q.top().second), vis, ubound - 1);

        if (found && op.size() + 1 < ubound)
        {
            // Neue kürzeste Folge gefunden.
            ubound = op.size() + 1;
            res = op;
            res.push_back(q.top().second);
            found_solution = 1;
        }

        q.pop();
    }

    vis[p.size() - 1].insert(ind(p));
    return {res, found_solution};
}

vector<unsigned> shortest_op_dfs(vector<unsigned> const &p)
{
    vector<unordered_set<size_t>> vis(p.size());
    vector<unsigned> res = shortest_op_dfs_r(p, vis).first;
    reverse(res.begin(), res.end());
    return res;
}

struct node
{
    size_t n, i;
    unsigned lbound;
    vector<unsigned> op;

    bool operator<(node const &x) const
    {
        if (lbound == x.lbound)
            return n > x.n;
        return lbound > x.lbound;
    }
};

vector<unsigned> shortest_op_bfs(vector<unsigned> const &p)
{
    // Speichert für jede Permutationsgröße die Indizes besuchter Permutationen.
    vector<unordered_set<size_t>> vis(p.size());
    vis[p.size() - 1].insert(ind(p));

    priority_queue<node> q;
    q.push({p.size(), ind(p), get_lbound(p), {}});

    unsigned ubound = p.size(); // aktuelle Oberschranke
    vector<unsigned> res;

    while (!q.empty() && q.top().lbound < ubound)
    {
        node const x = q.top();
        q.pop();

        if (!x.i)
        {
            // Eine identische Permutation wurde gefunden.
            if (x.op.size() < ubound)
            {
                res = x.op;
                ubound = res.size();
            }
            continue;
        }

        vector<unsigned> s = ith_permutation(x.n, x.i);

        // Füge jede durch eine gamma-Operation erreichbare Permutation zur
        // Warteschlange hinzu, die das Ergebnis noch verbessern kann.
        for (size_t i = 0; i < x.n; i++)
        {
            unsigned const lbound = x.op.size() + get_lbound_gamma(s, i) + 1;
            if (lbound < ubound)
            {
                node y = {x.n - 1, ind_gamma(s, i), lbound, x.op};
                y.op.push_back(i);

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

pair<vector<unsigned>, bool> shortest_op_bf_r(
    vector<unsigned> const &p, vector<unordered_set<size_t>> &vis)
{
    if (!ind(p))
        return {{}, 1};

    if (vis[p.size() - 1].find(ind(p)) != vis[p.size() - 1].end())
        return {{}, 0};

    vector<unsigned> res;
    bool found_solution = 0;

    for (size_t i = 0; i < p.size(); i++)
    {
        auto const [op, found] = shortest_op_bf_r(gamma(p, i), vis);

        if (found && (op.size() + 1 < res.size() || !found_solution))
        {
            res = op;
            res.push_back(i);
            found_solution = 1;
        }
    }

    vis[p.size() - 1].insert(ind(p));
    return {res, found_solution};
}

// Findet die kürzeste Folge an gamma-Operationen, um p in eine identische
// Permutation umzuformen, durch Austesten aller möglichen Folgen.
vector<unsigned> shortest_op_bf(vector<unsigned> const &p)
{
    vector<unordered_set<size_t>> vis(p.size());
    vector<unsigned> res = shortest_op_bf_r(p, vis).first;
    reverse(res.begin(), res.end());
    return res;
}

int main(int argc, char *argv[])
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
    vector<unsigned> op;

    if (argc == 2 && !strcmp(argv[1], "--dfs"))
        op = shortest_op_dfs(p);
    else if (argc == 2 && !strcmp(argv[1], "--bf"))
        op = shortest_op_bf(p);
    else
        op = shortest_op_bfs(p);

    if (!op.empty())
    {
        cout << op.size() << " Operationen nötig. Dafür hinter folgenden"
                             " Indizes wenden:\n\n";
        cout << "Index |  p\n";

        for (auto it = op.cbegin(); it != op.cend(); it++)
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