#include <bits/stdc++.h>
#include "aufgabe3_a.hpp"
#include "util.hpp"
using namespace std;

struct Node // Ein Knoten im Suchbaum. Enthält Index, Länge und untere Schranke.
{           // Bei Vergleich werden untere Schranke und Länge verglichen.
    size_t index;
    unsigned length, lbound;

    Node(size_t index_, unsigned length_, unsigned lbound_)
    {
        index = index_, length = length_, lbound = lbound_;
    }

    // Umgekehrter Vergleich wegen absteigender Sortierung in der priority_queue
    bool operator<(Node const &x) const
    {
        return lbound == x.lbound ? length > x.length : lbound > x.lbound;
    }
};

// Findet die kürzeste Folge an gamma-Operationen, um p in eine identische
// Permutation umzuformen durch Austesten aller möglichen Operationen.
vector<unsigned> min_operations_bfs(vector<unsigned> const &p)
{
    vector<unordered_map<size_t, size_t>> pre(p.size());
    pre[p.size() - 1][ind(p)] = SIZE_MAX;

    queue<Node> q;
    q.emplace(ind(p), p.size(), 0); // Untere Schranke wird nicht verwendet.

    while (!q.empty())
    {
        auto const [index, length, _] = q.front();
        q.pop();

        if (!index) // Identische Permutation gefunden
            return reconstruct_operations(pre, length, index);

        vector<unsigned> const s = ith_permutation(length, index);

        for (size_t j = 0; j < length; j++) // Wende alle möglichen gamma_i-
        {                                   // Operationen (0 <= i < length) an.
            Node const y(ind_gamma(s, j), length - 1, 0);
            if (pre[y.length - 1].find(y.index) == pre[y.length - 1].end())
            {
                pre[y.length - 1][y.index] = index; // Vorgänger im Suchbaum.
                q.push(y);
            }
        }
    }

    exit(EXIT_FAILURE); // Sollte nie erreicht werden.
}

// Gibt zurück, ob a < b und erhäht x um 1, wenn sich das Steigungsverhalten
// ändert.
bool is_increasing(bool incr, unsigned a, unsigned b, unsigned &x)
{
    if (incr ^ (a < b))
    {
        x++;
        return !incr;
    }
    return incr;
}

// Bestimmt eine untere Schrankte für A(p).
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

// Bestimmt eine untere Schranke für A(gamma_i p).
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

// Gibt die kürzestmögliche Folge an gamma-Operationen zurück. Wie im A*-
// Algorithmus werden die Blätter der Suche in einer Prioritätswarteschlange
// gespeichert und aufsteigend nach unterer Schranke abgearbeitet.
vector<unsigned> min_operations_astar(vector<unsigned> const &p)
{
    priority_queue<Node> q;
    q.emplace(ind(p), (unsigned)p.size(), get_lbound(p));

    // Speichert für jede Permutationslänge die Indizes besuchter
    // Permutationen und deren Vorgänger.
    vector<unordered_map<size_t, size_t>> pre(p.size());
    pre[p.size() - 1][ind(p)] = SIZE_MAX;

    unsigned ubound = p.size(); // aktuelle Oberschranke
    size_t n_res = SIZE_MAX;

    while (!q.empty() && q.top().lbound < ubound)
    {
        auto const [index, length, lbound] = q.top();
        q.pop();

        if (!index)
        {
            // Eine identische Permutation wurde gefunden.
            if (p.size() - length < ubound)
            {
                n_res = length;
                ubound = p.size() - length;
            }
            continue;
        }

        vector<unsigned> const s = ith_permutation(length, index);

        // Füge jede durch eine gamma-Operation erreichbare Permutation zur
        // Warteschlange hinzu, die das Ergebnis noch verbessern kann.
        for (size_t j = 0; j < length; j++)
        {
            Node const y(ind_gamma(s, j), length - 1,
                         p.size() - length + get_lbound_gamma(s, j) + 1);

            if (pre[length - 2].find(y.index) == pre[length - 2].end() &&
                y.lbound < ubound)
            {
                pre[length - 2][y.index] = index;
                q.push(y);
            }
        }
    }

    return reconstruct_operations(pre, n_res, 0);
}

// Bestimmt rekursiv die kürzestmögliche Folge an gamma-Operationen zum
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

    priority_queue<Node> q;
    for (size_t i = 0; i < p.size(); i++)
        q.emplace(i, p.size() - 1, get_lbound_gamma(p, i));

    vector<unsigned> res;
    bool found_better = 0;

    // Gehe die Nachfolgerknoten aufsteigend nach unterer Schranke durch und
    // bestimme rekursiv die kürzeste Operationenfolge.
    while (!q.empty() && q.top().lbound < ubound)
    {
        auto const [index, lbound, length] = q.top();
        q.pop();

        auto const [op, found] = // Rekursionsschritt
            min_operations_bnb_r(gamma(p, index), vis, ubound - 1);

        if (found && op.size() + 1 < ubound)
        {
            // Neue kürzeste Folge gefunden.
            ubound = op.size() + 1;
            res = op;
            res.push_back(index);
            found_better = 1;
        }
    }

    vis[p.size() - 1].insert(ind(p));
    return {res, found_better};
}

// Stellt vis für die rekursive Funktion min_operations_bnb_r bereit und bringt
// die Operationen in die richtige Reihenfolge.
vector<unsigned> min_operations_bnb(vector<unsigned> const &p)
{
    vector<unordered_set<size_t>> vis(p.size());
    vector<unsigned> res = min_operations_bnb_r(p, vis).first;
    reverse(res.begin(), res.end());
    return res;
}

// Läuft den durch pre gegebenen Suchbaum hoch und sammelt die gamma-
// Operationen auf dem Pfad ein. m ist die Länge der anfänglichen Permutation,
// i ihr Index.
vector<unsigned> reconstruct_operations(
    vector<unordered_map<size_t, size_t>> const &pre, size_t m, size_t i)
{
    vector<unsigned> t;

    while (m < pre.size())
    {
        vector<unsigned> const s = ith_permutation(m + 1, pre[m - 1].at(i));

        // Suche nach der gamma-Operation, die den Vorgänger (p) in den
        // Nachfolger (i-te Permutation der Länge n) umwandelt.
        for (size_t j = 0; j < m + 1; j++)
            if (ind_gamma(s, j) == i)
            {
                t.push_back(j);
                break;
            }

        i = pre[m - 1].at(i); // Speichere den Index der vorherigen Permutation
        m++;                  // und erhöhe die Permutationslänge um 1.
    }

    reverse(t.begin(), t.end()); // Die Operationen wurden umgekehrt eingefügt.
    return t;
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
            p = reverse_and_eat(p, *it);
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