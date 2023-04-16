#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

vector<complex<double>> randomized_obtuse_path(vector<complex<double>> const &z)
{
    size_t const n = z.size(), sqrtn = sqrt(n);
    deque<size_t> path;
    list<size_t> unvisited;
    bool front_is_dead_end, back_is_dead_end, extending_back;
    size_t no_added_node;

    auto restart_search = [&]()
    {
        front_is_dead_end = back_is_dead_end = 0; // Setze alle Datenstrukturen
        extending_back = 0;                       // zurück.
        unvisited.clear();
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
        if (no_added_node > sqrtn) // Zu große Anzahl aufeinanderfolgener
            restart_search();      // Iterationen ohne Hinzufügen eines Knotens.

        // u: letzter / erster Knoten, v: vorletzter / zweiter Knoten
        // w: Iterator in unvisited zum neu hinzugefügten Knoten
        size_t u, v, candidates = 0;
        list<size_t>::iterator w = unvisited.end();

        if (extending_back)
            u = path.back(), v = path.size() >= 2 ? *++path.rbegin() : SIZE_MAX;
        else
            u = path.front(), v = path.size() >= 2 ? *++path.begin() : SIZE_MAX;

        for (auto it = unvisited.begin(); it != unvisited.end(); it++)
            if (path.size() < 2 || dot_product(z[u] - z[v], z[*it] - z[u]) >= 0)
            {
                candidates++;
                if (!(rand() % candidates)) // wahr mit Wahrscheinlichkeit
                    w = it;                 // 1 / candidates.
                if (candidates > sqrtn)
                    break;
            }
        if (candidates)
        {
            if (extending_back)
                path.push_back(*w); // Erweitere den Pfad um w.
            else
                path.push_front(*w);
            unvisited.erase(w);
            no_added_node = 0;
            continue;
        }

        // Der Pfad kann von u aus nicht mehr erweitert werden. Das aktuell
        // behandelte Ende des Pfads wird als Sackgasse markiert.
        (extending_back ? back_is_dead_end : front_is_dead_end) = 1;
        no_added_node++; // In dieser Iteration wurde kein Knoten hinzugefügt.

        if (front_is_dead_end && back_is_dead_end)
        {
            // i: Index des aktuell betrachteten Knotens
            // w: Index des Knotens, nach dem der Pfad aufgebrochen wird.
            size_t i = extending_back ? path.size() - 3 : 2, w = SIZE_MAX,
                   candidates = 0;

            while (i < path.size())
            {
                // j: Nachbar von path[i], mit dem beim Einfügen der Kante
                //    {u, path[i]} ein Abbiegewinkel entstehen würde.
                size_t j = extending_back ? i - 1 : i + 1;
                if (dot_product(z[u] - z[v], z[path[i]] - z[u]) >= 0 &&
                    (j >= path.size() ||
                     dot_product(z[path[i]] - z[u], z[path[j]] - z[path[i]]) >= 0))
                {
                    candidates++;
                    if (!(rand() % candidates)) // wahr mit Wahrscheinlichkeit
                        w = i;                  // 1 / candidates.
                    if (candidates > sqrtn)
                        break;
                }
                i += extending_back ? -1 : 1;
            }

            if (candidates) // Breche den Pfad nach w auf und verbinde u mit w.
            {
                if (extending_back) // Kehre das Suffix strikt nach w um.
                    reverse(path.begin() + w + 1, path.end());
                else // Kehre das Präfix strikt vor w um.
                    reverse(path.begin(), path.begin() + w);
                (extending_back ? back_is_dead_end : front_is_dead_end) = 0;
            }
            else
                extending_back = !extending_back;
        }
        else // Versuche den Pfad am anderen Ende zu erweitern.
            extending_back = !extending_back;
    }

    vector<complex<double>> point_order;
    for (size_t i = 0; i < n; i++)
        point_order.push_back(z[path[i]]);
    return point_order;
}

// Optimiert den gegebenen Pfad mit der 2-opt Heuristik, ohne die Beschränkung
// von Abbiegewinkeln zu verletzen. Mit time_limit wird ein Zeitlimit in s
// gesetzt.
vector<complex<double>> optimize_path(
    vector<complex<double>> const &path, double time_limit)
{
    auto start_time = chrono::system_clock::now();

    size_t const n = path.size();
    vector<array<size_t, 2>> nodes(n);
    queue<pair<size_t, size_t>> q;
    for (size_t i = 0; i + 1 < n; i++)
    { // Füge alle Knoten in die Warteschlange ein.
        q.emplace(i, i + 1), q.emplace(i + 1, i);
        nodes[i][1] = i + 1, nodes[i + 1][0] = i;
    }
    nodes[0][0] = nodes[n - 1][1] = SIZE_MAX;

    while (!q.empty() &&
           chrono::duration<double>(chrono::system_clock::now() - start_time)
                   .count() < time_limit)
    {
        auto const [v, w] = q.front();
        q.pop();
        if (nodes[v][0] != w && nodes[v][1] != w) // Überprüfe, ob es die Kante
            continue;                             // noch gibt.
        bool const direction = nodes[w][0] == v;

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
            if ((u == SIZE_MAX ||
                 dot_product(path[v] - path[u], path[b] - path[v]) >= 0) &&
                dot_product(path[b] - path[v], path[a] - path[b]) >= 0 &&
                dot_product(path[w] - path[x], path[c] - path[w]) >= 0 &&
                (d == SIZE_MAX ||
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

    vector<complex<double>> new_path; // nue
    size_t start = SIZE_MAX;
    bool direction = 0;
    for (size_t i = 0; i < n; i++) // Suche nach einem Knoten mit Grad 1.
        if (nodes[i][0] == SIZE_MAX || nodes[i][1] == SIZE_MAX)
        {
            start = i; // Startknoten gefunden.
            direction = nodes[i][1] != SIZE_MAX;
            break;
        }
    for (size_t i = start; i != SIZE_MAX; i = nodes[i][direction])
        new_path.push_back(path[i]);
    return new_path;
}

int main(int argc, char **argv)
{
    vector<complex<double>> z;
    double x, y;
    while (scanf("%lf %lf", &x, &y) == 2)
        z.emplace_back(x, y);

    auto start_time = chrono::system_clock::now();

    vector<complex<double>> path = randomized_obtuse_path(z);
    auto duration = chrono::duration_cast<chrono::duration<double>>(
        chrono::system_clock::now() - start_time);

    cerr << setprecision(6) << fixed << "Zulässige Tour mit Länge "
         << path_length(path) << " nach " << duration.count()
         << " s gefunden.\nStarte 2-opt..." << endl;

    double two_opt_time_limit = DBL_MAX;
    for (int i = 1; i + 1 < argc; i++)
        if (!strcmp(argv[i], "--2-opt-time-limit"))
            two_opt_time_limit = stod(argv[i + 1]);

    path = optimize_path(path, two_opt_time_limit);
    duration = chrono::duration_cast<chrono::duration<double>>(
        chrono::system_clock::now() - start_time);
    cerr << "Tourlänge nach Optimierung: " << path_length(path) << '\n'
         << "Laufzeit: " << duration.count() << " s" << endl;

    cout << setprecision(6) << fixed
         << "Tourlänge: " << path_length(path) << '\n';
    for (complex<double> const &u : path)
        cout << u.real() << ' ' << u.imag() << '\n';
}