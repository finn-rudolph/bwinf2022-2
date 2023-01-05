#include <bits/stdc++.h>
#include "util.hpp"
#include "aufgabe3_a.hpp"
using namespace std;

// Findet die kuerzeste Folge an gamma-Operationen, um p in eine identische
// Permutation umzuformen, durch Austesten aller moeglichen Folgen.
vector<unsigned> min_operations_bfs(vector<unsigned> const &p)
{
    vector<unordered_map<size_t, size_t>> pre(p.size());
    pre[p.size() - 1][ind(p)] = SIZE_MAX;

    // Speichert Laenge und Index jedes Blatts im Suchbaum.
    queue<pair<size_t, size_t>> q;
    q.push({p.size(), ind(p)});

    while (!q.empty())
    {
        auto const [m, i] = q.front();
        q.pop();

        if (!i)
            return reconstruct_operations(pre, m, i);

        vector<unsigned> const s = ith_permutation(m, i);

        for (size_t j = 0; j < m; j++)
        {
            pair<size_t, size_t> const y = {m - 1, ind_gamma(s, j)};
            if (pre[y.first - 1].find(y.second) == pre[y.first - 1].end())
            {
                pre[y.first - 1][y.second] = i;
                q.push(y);
            }
        }
    }

    exit(EXIT_FAILURE); // Sollte nie erreicht werden.
}

// Laeuft den durch pre gegebenen Suchbaum hoch und sammelt die gamma-
// Operationen auf dem Pfad ein.
vector<unsigned> reconstruct_operations(
    vector<unordered_map<size_t, size_t>> const &pre, size_t m, size_t i)
{
    vector<unsigned> t;

    while (m < pre.size())
    {
        vector<unsigned> const s = ith_permutation(m + 1, pre[m - 1].at(i));

        // Suche nach der gamma-Operation, die den Vorgaenger (p) in den
        // Nachfolger (i-te Permutation der Laenge n) umwandelt.
        for (size_t j = 0; j < m + 1; j++)
            if (ind_gamma(s, j) == i)
            {
                t.push_back(j);
                break;
            }

        i = pre[m - 1].at(i);
        m++;
    }

    reverse(t.begin(), t.end());
    return t;
}

bool is_increasing(bool incr, unsigned a, unsigned b, unsigned &x)
{
    if (incr ^ (a < b))
    {
        x++;
        return !incr;
    }
    return incr;
}

unsigned get_lbound(vector<unsigned> const &p)
{
    if (p.size() == 1)
        return 0;

    unsigned x = 1; // Anzahl monotoner Teilstrings
    bool incr = p[1] > p[0];

    for (size_t i = 2; i < p.size(); i++)
        incr = is_increasing(incr, p[i - 1], p[i], x);

    return (x + 1) / 3;
}

// Bestimmt die untere Schranke nach der Anwendung von gamma_i auf p.
unsigned get_lbound_gamma(vector<unsigned> const &p, size_t i)
{
    assert(i < p.size());

    if (p.size() <= 2)
        return 0;

    unsigned x = 1;
    bool incr =
        i >= 2 ? (p[i - 2] > p[i - 1])
               : (i == 1 ? (p[i + 1] > p[i - 1]) : (p[i + 2] > p[i + 1]));

    for (size_t j = 2; j < i; j++) // umgekehrte Elemente vor i
        incr = is_increasing(incr, p[i - j], p[i - j - 1], x);

    if (i && i + 1 < p.size()) // neu benachbarte Elemente (p[0], p[i + 1])
        incr = is_increasing(incr, p[0], p[i + 1], x);

    for (size_t j = i + 2; j < p.size(); j++) // Elemente nach i
        incr = is_increasing(incr, p[j - 1], p[j], x);

    return (x + 1) / 3;
}

// Speichert Index, Laenge, untere Schranke.
typedef tuple<size_t, unsigned, unsigned> node;

inline bool node_compare(node const &x, node const &y)
{
    if (get<2>(x) == get<2>(y))
    {
        return get<1>(x) > get<1>(y);
    }
    return get<2>(x) > get<2>(y);
}

// Gibt die kuerzestmoegliche Folge an gamma-Operationen zurueck. Wie im A*-
// Algorithmus werden die Blaetter der Suche in einer Prioritaetswarteschlange
// gespeichert und aufsteigend nach unterer Schranke abgearbeitet.
vector<unsigned> min_operations_astar(vector<unsigned> const &p)
{
    priority_queue<node, vector<node>, decltype(&node_compare)>
        q(&node_compare);
    q.push({ind(p), (unsigned)p.size(), get_lbound(p)});

    // Speichert fuer jede Permutationslaenge die Indizes besuchter
    // Permutationen und deren Vorgaenger im Suchbaum.
    vector<unordered_map<size_t, size_t>> pre(p.size());
    pre[p.size() - 1][ind(p)] = SIZE_MAX;

    unsigned ubound = p.size(); // aktuelle Oberschranke
    size_t n_res = SIZE_MAX;

    while (!q.empty() && get<2>(q.top()) < ubound)
    {
        auto const [i, m, lbound] = q.top();
        q.pop();

        if (!i)
        {
            // Eine identische Permutation wurde gefunden.
            if (p.size() - m < ubound)
            {
                n_res = m;
                ubound = p.size() - m;
            }
            continue;
        }

        vector<unsigned> const s = ith_permutation(m, i);

        // Fuege jede durch eine gamma-Operation erreichbare Permutation zur
        // Warteschlange hinzu, die das Ergebnis noch verbessern kann.
        for (size_t j = 0; j < m; j++)
        {
            node const y = {ind_gamma(s, j), m - 1,
                            p.size() - m + get_lbound_gamma(s, j) + 1};

            if (pre[m - 2].find(get<0>(y)) == pre[m - 2].end() &&
                get<2>(y) < ubound)
            {
                pre[m - 2][get<0>(y)] = i;
                q.push(y);
            }
        }
    }

    return reconstruct_operations(pre, n_res, 0);
}

// Bestimmt rekursiv die kuerzestmoegliche Folge and gamma-Operationen zum
// Sortieren von p. Falls keine kuerzere als ubound existiert, wird das zweite
// Element der Rueckgabe auf 0 gesetzt.
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
    // rekursiv die kuerzeste Operationenfolge.
    while (!q.empty() && q.top().first < ubound)
    {
        auto const [lbound, i] = q.top();
        q.pop();

        auto const [op, found] =
            min_operations_bnb_r(gamma(p, i), vis, ubound - 1);

        if (found && op.size() + 1 < ubound)
        {
            // Neue kuerzeste Folge gefunden.
            ubound = op.size() + 1;
            res = op;
            res.push_back(i);
            found_better = 1;
        }
    }

    vis[p.size() - 1].insert(ind(p));
    return {res, found_better};
}

// Sollte zum Bestimmen von A(p) aufgerufen werden. Stellt vis fuer die
// rekursive Funktion bereit und bringt die Operationen in die richtige
// Reihenfolge.
vector<unsigned> min_operations_bnb(vector<unsigned> const &p)
{
    vector<unordered_set<size_t>> vis(p.size());
    vector<unsigned> res = min_operations_bnb_r(p, vis).first;
    reverse(res.begin(), res.end());
    return res;
}

void print_operations(vector<unsigned> p, vector<unsigned> const &op)
{
    if (!op.empty())
    {
        cout << op.size() << " Operationen notwendig. Dazu hinter folgenden"
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

    print_operations(p, op);
}