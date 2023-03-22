#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

constexpr size_t TwoOptIterationLimit = 1000000;

vector<complex<double>> randomized_obtuse_path(vector<complex<double>> const &z)
{
    size_t const n = z.size();
    vector<size_t> path;
    list<size_t> unvisited;
    bool front_is_dead_end, back_is_dead_end;
    size_t no_added_node;

    auto restart_search = [&]()
    {
        front_is_dead_end = back_is_dead_end = 0; // Setze alle Datenstrukturen
        unvisited.clear();                        // zurück.
        path.clear();
        no_added_node = 0;
        srand(time(0));
        path.push_back(rand() % n);    // Wähle einen zufälligen Startknoten.
        for (size_t i = 0; i < n; i++) // Fülle unvisited mit allen Knoten
            if (i != path.front())     // außer dem Startknoten.
                unvisited.push_back(i);
    };

    restart_search();

    while (path.size() < n)
    {
        if (no_added_node > n / 2) // Zu große Anzahl aufeinanderfolgener
            restart_search();      // Iterationen ohne Hinzufügen eines Knotens.

        // u: letzter Knoten, v: vorletzter Knoten (wenn existent)
        // w: Iterator in unvisited zum neu hinzugefügten Knoten
        size_t u = path.back(),
               v = path.size() >= 2 ? *++path.rbegin() : SIZE_MAX,
               candidates = 0;
        list<size_t>::iterator w = unvisited.end();

        for (auto it = unvisited.begin(); it != unvisited.end(); it++)
            if (path.size() < 2 || dot_product(z[u] - z[v], z[*it] - z[u]) >= 0)
            {
                candidates++;
                if (!(rand() % candidates)) // wahr mit Wahrscheinlichkeit
                    w = it;                 // 1 / candidates.
            }
        if (candidates)
        {
            path.push_back(*w); // Erweitere den Pfad um w.
            unvisited.erase(w);
            no_added_node = 0;
            continue;
        }

        // Der Pfad kann von u aus nicht mehr erweitert werden, da alle mit
        // Abbiegewinkel <= pi / 2 schon besucht wurden. Das Ende des Pfads wird
        // als Sackgasse markiert.
        back_is_dead_end = 1;
        no_added_node++; // In dieser Iteration wurde kein Knoten hinzugefügt.

        if (front_is_dead_end && back_is_dead_end)
        {
            size_t candidates = 0;
            // it: Iterator zum aktuell betrachteten Knoten
            // w: Iterator zum Knoten, nach dem der Pfad aufgebrochen wird.
            auto it = path.rbegin() + 2, w = path.rend();

            while (it != path.rend()) // Überprüfe die Winkelbeschränkungen und
            {                         // ziehe it als Kandidaten in Betracht.
                if (dot_product(z[u] - z[v], z[*it] - z[u]) >= 0 &&
                    (it + 1 == path.rend() ||
                     dot_product(z[*it] - z[u], z[*(it + 1)] - z[*it]) >= 0))
                {
                    candidates++;
                    if (!(rand() % candidates)) // wahr mit Wahrscheinlichkeit
                        w = it;                 // 1 / candidates.
                }
                it++;
            }

            if (candidates) // Breche den Pfad nach w auf und verbinde u mit w.
            {
                reverse(path.rbegin(), w); // Kehre das Suffix bis w um.
                back_is_dead_end = 0;
            }
            else
                reverse(path.begin(), path.end());
        }
        else // Versuche den Pfad am anderen Ende zu erweitern.
        {
            reverse(path.begin(), path.end());
            swap(front_is_dead_end, back_is_dead_end);
        }
    }

    vector<complex<double>> point_order;
    for (size_t i = 0; i < n; i++)
        point_order.push_back(z[path[i]]);
    return point_order;
}

// Gibt einen mit der 2-opt Heuristik optimierten, ohne
// Abbiegewinkel > pi / 2 zu erzeugen.
vector<complex<double>> optimize_path(vector<complex<double>> const &path)
{
    size_t const n = path.size();
    vector<array<size_t, 2>> nodes(n);
    queue<pair<size_t, size_t>> q;
    for (size_t i = 0; i + 1 < n; i++)
    { // Füge alle Knoten in die Warteschlange ein.
        q.emplace(i, i + 1), q.emplace(i + 1, i);
        nodes[i][1] = i + 1, nodes[i + 1][0] = i;
    }
    nodes[0][0] = nodes[n - 1][1] = SIZE_MAX;

    size_t iteration_count = 0;
    while (!q.empty() && iteration_count < TwoOptIterationLimit)
    {
        auto const [v, w] = q.front();
        q.pop();
        if (nodes[v][0] != w && nodes[v][1] != w)
            continue;
        bool const direction = nodes[w][0] == v;
        iteration_count++;

        // Die aktuell bearbeitete Kante ist {v, w}. u kommt vor v, x nach w.
        // Als mögliche Tauschpartner werden nur Kanten in Richtung von x in
        // Betracht gezogen.
        size_t const u = nodes[v][!direction], x = nodes[w][direction];
        size_t a = w, b = x;

        while (b != SIZE_MAX && nodes[b][direction] != SIZE_MAX)
        {
            size_t c = nodes[b][direction], d = nodes[c][direction];

            // Überprüfe, ob die Tour durch Erstetzen von {v, w}, {b, c} durch
            // {v, b}, {w, c} verkürzt wird und die 4 neuen Abbiegewinkel (uvb,
            // vba, xwc, wcd) alle <= pi / 2 sind.
            if ((u >= n ||
                 dot_product(path[v] - path[u], path[b] - path[v]) >= 0) &&
                dot_product(path[b] - path[v], path[a] - path[b]) >= 0 &&
                dot_product(path[w] - path[x], path[c] - path[w]) >= 0 &&
                (d >= n ||
                 dot_product(path[c] - path[w], path[d] - path[c]) >= 0) &&
                abs(path[v] - path[w]) + abs(path[b] - path[c]) >
                    abs(path[v] - path[b]) + abs(path[w] - path[c]))
            {
                // Kehre den Teil des Pfads von x bis a um.
                for (size_t i = x; i != b; i = nodes[i][!direction])
                    swap(nodes[i][0], nodes[i][1]);
                // Entferne {v, w}, {b, c} und füge {v, b}, {w, c} ein.
                nodes[v][direction] = b;
                nodes[b][direction] = a;
                nodes[b][!direction] = v;
                nodes[w][!direction] = x;
                nodes[w][direction] = c;
                nodes[c][!direction] = w;
                // Die neu eingefügten Kanten können erneut mit anderen
                // vertauscht werden, daher werden sie zu q hinzugefügt.
                q.emplace(v, b);
                q.emplace(b, v);
                q.emplace(w, c);
                q.emplace(c, w);
                break;
            }
            a = nodes[a][direction]; // Gehe zur nächsten Kante im Pfad.
            b = nodes[b][direction];
        }
    }

    vector<complex<double>> new_path;
    size_t start;
    bool direction;
    for (size_t i = 0; i < n; i++) // Suche nach einem Knoten mit Grad 1.
        if (nodes[i][0] == SIZE_MAX || nodes[i][1] == SIZE_MAX)
        {
            start = i;
            direction = nodes[i][1] != SIZE_MAX;
            break;
        }
    for (size_t i = start; i != SIZE_MAX; i = nodes[i][direction])
        new_path.push_back(path[i]);
    return new_path;
}

int main()
{
    vector<complex<double>> z;
    double x, y;
    while (scanf("%lf %lf", &x, &y) == 2)
        z.emplace_back(x, y);

    vector<complex<double>> path = randomized_obtuse_path(z);
    cerr << "Zulässige Tour mit Länge " << path_length(path) << " gefunden.\n"
         << "Starte 2-opt...\n";
    path = optimize_path(path);
    cerr << "Tourlänge nach Optimierung: " << path_length(path) << '\n';

    cout << setprecision(6) << fixed
         << "Tourlänge: " << path_length(path) << '\n';
    for (complex<double> const &u : path)
        cout << u.real() << ' ' << u.imag() << '\n';
}