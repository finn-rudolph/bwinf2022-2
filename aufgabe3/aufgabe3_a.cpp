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

// Bestimmt die untere Schranke nach der Anwendung von gamma_i auf p.
unsigned get_lbound_gamma(vector<unsigned> const &p, size_t i)
{
    assert(i < p.size());

    if (p.size() <= 2)
        return 0;

    bool incr =
        i >= 2 ? (p[i - 2] > p[i - 1])
               : (i == 1 ? (p[i + 1] > p[i - 1]) : (p[i + 2] > p[i + 1]));
    unsigned num_monotone = 1;

    for (size_t j = 2; j < i; j++) // umgekehrte Elemente vor i
        incr = is_increasing(incr, p[i - j], p[i - j - 1], num_monotone);

    if (i && i + 1 < p.size()) // neu benachbarte Elemente (p[0], p[i + 1])
        incr = is_increasing(incr, p[0], p[i + 1], num_monotone);

    for (size_t j = i + 2; j < p.size(); j++) // Elemente nach i
        incr = is_increasing(incr, p[j - 1], p[j], num_monotone);

    return (num_monotone + 1) / 3;
}

// Bestimmt rekursiv die kürzestmögliche Folge and gamma-Operationen zum
// Sortieren von p. Falls keine kürzere als ubound existiert, wird das zweite
// Element der Rückgabe auf 0 gesetzt.
pair<vector<unsigned>, bool> min_operations_bnb_r(
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
    bool found_better = 0;

    // Gehe die Nachfolgerknoten (Permutationen, die durch eine gamma-Operation
    // erreichbar sind) aufsteigend nach unterer Schranke durch und bestimme
    // rekursiv die kürzeste Operationenfolge.
    while (!q.empty() && q.top().first < ubound)
    {
        auto const [lbound, i] = q.top();
        q.pop();

        auto const [op, found] =
            min_operations_bnb_r(gamma(p, i), vis, ubound - 1);

        if (found && op.size() + 1 < ubound)
        {
            // Neue kürzeste Folge gefunden.
            ubound = op.size() + 1;
            res = op;
            res.push_back(i);
            found_better = 1;
        }
    }

    vis[p.size() - 1].insert(ind(p));
    return {res, found_better};
}

vector<unsigned> min_operations_bnb(vector<unsigned> const &p)
{
    vector<unordered_set<size_t>> vis(p.size());
    vector<unsigned> res = min_operations_bnb_r(p, vis).first;
    reverse(res.begin(), res.end());
    return res;
}

// Läuft den durch pre gegebenen Suchbaum hoch und sammelt die gamma-Operationen
// auf dem Pfad ein.
vector<unsigned> reconstruct_op(
    vector<unordered_map<size_t, size_t>> const &pre, size_t m, size_t i)
{
    vector<unsigned> op;

    while (m < pre.size())
    {
        vector<unsigned> const s = ith_permutation(m + 1, pre[m - 1].at(i));

        // Suche nach der gamma-Operation, die den Vorgänger (p) in den
        // Nachfolger (i-te Permutation der Länge n) umwandelt.
        for (size_t j = 0; j < m + 1; j++)
            if (ind_gamma(s, j) == i)
            {
                op.push_back(j);
                break;
            }

        i = pre[m - 1].at(i);
        m++;
    }

    reverse(op.begin(), op.end());
    return op;
}

// Speichert Index, Länge, untere Schranke.
typedef tuple<size_t, unsigned, unsigned> node;

inline bool node_compare(node const &x, node const &y)
{
    if (get<2>(x) == get<2>(y))
    {
        return get<1>(x) > get<1>(y);
    }
    return get<2>(x) > get<2>(y);
}

// Gibt die kürzestmögliche Folge an gamma-Operationen zurück. Wie im A*-
// Algorithmus werden die Blätter der Suche in einer Prioritätswarteschlange
// gespeichert und aufsteigend nach unterer Schranke abgearbeitet.
vector<unsigned> min_operations_astar(vector<unsigned> const &p)
{
    // Speichert für jede Permutationsgröße die Indizes besuchter Permutationen
    // und deren Vorgänger im Suchbaum.
    vector<unordered_map<size_t, size_t>> pre(p.size());
    pre[p.size() - 1][ind(p)] = SIZE_MAX;

    priority_queue<node, vector<node>, decltype(&node_compare)>
        q(&node_compare);
    q.push({ind(p), (unsigned)p.size(), get_lbound(p)});

    unsigned ubound = p.size(); // aktuelle Oberschranke
    size_t res_n = SIZE_MAX;

    while (!q.empty() && get<2>(q.top()) < ubound)
    {
        auto const [i, m, lbound] = q.top();
        q.pop();

        if (!i)
        {
            // Eine identische Permutation wurde gefunden.
            if (p.size() - m < ubound)
            {
                res_n = m;
                ubound = p.size() - m;
            }
            continue;
        }

        vector<unsigned> const s = ith_permutation(m, i);

        // Füge jede durch eine gamma-Operation erreichbare Permutation zur
        // Warteschlange hinzu, die das Ergebnis noch verbessern kann.
        for (size_t j = 0; j < m; j++)
        {
            node const y = {ind_gamma(s, j), m - 1,
                            p.size() - m + get_lbound_gamma(s, j) + 1};

            if (get<2>(y) < ubound &&
                pre[m - 2].find(get<0>(y)) == pre[m - 2].end())
            {
                q.push(y);
                pre[m - 2][get<0>(y)] = i;
            }
        }
    }

    return reconstruct_op(pre, res_n, 0);
}

// Findet die kürzeste Folge an gamma-Operationen, um p in eine identische
// Permutation umzuformen, durch Austesten aller möglichen Folgen.
vector<unsigned> min_operations_bfs(vector<unsigned> const &p)
{
    vector<unordered_map<size_t, size_t>> pre(p.size());
    pre[p.size() - 1][ind(p)] = SIZE_MAX;

    // Speichert Länge und Index jedes Blatts im Suchbaum.
    queue<pair<size_t, size_t>> q;
    q.push({p.size(), ind(p)});

    while (!q.empty())
    {
        auto const [n, i] = q.front();
        q.pop();

        if (!i)
            return reconstruct_op(pre, n, i);

        vector<unsigned> const s = ith_permutation(n, i);

        for (size_t j = 0; j < n; j++)
        {
            pair<size_t, size_t> const y = {n - 1, ind_gamma(s, j)};
            if (pre[y.first - 1].find(y.second) == pre[y.first - 1].end())
            {
                pre[y.first - 1][y.second] = i;
                q.push(y);
            }
        }
    }

    exit(EXIT_FAILURE); // Sollte nie erreicht werden.
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

    vector<unsigned> op;

    if (argc == 2 && !strcmp(argv[1], "--bnb"))
        op = min_operations_bnb(p);
    else if (argc == 2 && !strcmp(argv[1], "--bfs"))
        op = min_operations_bfs(p);
    else
        op = min_operations_astar(p);

    if (!op.empty())
    {
        cout << op.size() << " Operationen nötig. Dafür hinter folgenden"
                             " Indizes wenden:\n";
        for (auto it = op.cbegin(); it != --op.cend(); it++)
            cout << *it << ", ";
        cout << *(--op.cend()) << "\n\n"
             << "Index |  p\n";

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