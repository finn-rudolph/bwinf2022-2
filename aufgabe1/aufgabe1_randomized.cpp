#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

constexpr size_t TwoOptIterationLimit = 1000000;

double path_length(vector<complex<double>> const &z)
{
    double length = 0.0;
    for (size_t i = 1; i < z.size(); i++)
        length += abs(z[i] - z[i - 1]);
    return length;
}

// Optimiert den gegebenen Hamiltonpfad in-place mit der 2-opt Heuristik, ohne
// Abbiegewinkel > pi / 2 zu erzeugen.
void optimize_path(vector<complex<double>> &path)
{
    size_t const n = path.size();
    array<vector<bool>, 2> processed;
    processed[0] = vector<bool>(n, 0);
    processed[1] = vector<bool>(n, 0);
    queue<pair<size_t, ptrdiff_t>> q;
    for (size_t i = 0; i + 1 < n; i++)
        q.emplace(i, -1), q.emplace(i + 1, 1);

    size_t iteration_count = 0;
    while (!q.empty() && iteration_count < TwoOptIterationLimit)
    {
        auto const [w, direction] = q.front();
        q.pop();
        if (processed[direction][w])
            continue;
        processed[(direction + 1) / 2][w] = 1;
        iteration_count++;

        size_t const v = w - direction, u = v - direction, x = w + direction;
        size_t a = w, b = x;

        while (b && b + 1 < n)
        {
            size_t c = b + direction, d = c + direction;

            // Überprüfe, ob die 4 neuen Abbiegewinkel (uvb, vba, xwc, wcd) alle
            // <= pi / 2 sind.
            if ((u >= n ||
                 dot_product(path[v] - path[u], path[b] - path[v]) >= 0) &&
                dot_product(path[b] - path[v], path[a] - path[b]) >= 0 &&
                dot_product(path[w] - path[x], path[c] - path[w]) >= 0 &&
                (d >= n ||
                 dot_product(path[c] - path[w], path[d] - path[c]) >= 0) &&
                abs(path[v] - path[w]) + abs(path[b] - path[c]) >
                    abs(path[v] - path[b]) + abs(path[w] - path[c]))
            {
                reverse(path.begin() + min(w, b), path.begin() + max(w, b) + 1);
                q.emplace(v, -direction); // Die neu eingefügten Kanten können
                q.emplace(b, direction);  // erneut mit anderen vertauscht
                q.emplace(w, -direction); // werden.
                q.emplace(c, direction);
                processed[(-direction + 1) / 2][v] = 0;
                processed[(direction + 1) / 2][b] = 0;
                processed[(-direction + 1) / 2][w] = 0;
                processed[(direction + 1) / 2][c] = 0;
                break;
            }
            a += direction; // Gehe zur nächsten Kante im Pfad.
            b += direction;
        }
    }
}

vector<complex<double>> obtuse_path(vector<complex<double>> const &z)
{
    size_t const n = z.size();
    vector<size_t> path;
    list<size_t> unvisited;
    bool front_is_dead_end, back_is_dead_end;
    size_t no_added_node;

    auto restart_search = [&]()
    {
        front_is_dead_end = back_is_dead_end = 0;
        unvisited.clear();
        path.clear();
        srand(time(0));
        path.push_back(rand() % n); // Wähle einen zufälligen Startknoten.
        for (size_t i = 0; i < n; i++)
            if (i != path.front())
                unvisited.push_back(i);
        no_added_node = 0;
    };

    restart_search();

    while (path.size() < n)
    {
        if (no_added_node > n / 2) // Zu große Anzahl aufeinanderfolgener
            restart_search();      // Iterationen ohne Hinzufügen eines Knotens.

        // u: erster Knoten, v: zweiter Knoten (wenn existent)
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
            auto it = path.rbegin() + 2, w = path.rend();

            while (it != path.rend())
            {
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

            if (w != path.rend())
            {
                reverse(path.rbegin(), w);
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

int main()
{
    vector<complex<double>> z;
    double x, y;
    while (scanf("%lf %lf", &x, &y) == 2)
        z.emplace_back(x, y);

    vector<complex<double>> path = obtuse_path(z);
    cerr << "Zulässige Tour mit Länge " << path_length(path) << " gefunden.\n"
         << "Starte 2-opt...\n";
    optimize_path(path);
    cerr << "Tourlänge nach Optimierung: " << path_length(path) << '\n';

    cout << setprecision(6) << fixed
         << "Tourlänge: " << path_length(path) << '\n';
    for (complex<double> const &u : path)
        cout << u.real() << ' ' << u.imag() << '\n';
}