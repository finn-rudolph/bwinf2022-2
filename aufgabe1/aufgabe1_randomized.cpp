#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

constexpr size_t TwoOptIterationLimit = 400000;

struct Node
{
    Node *adj[2]; // 0: Vorgänger, 1: Nachfolger.
    complex<double> p;

    Node() { adj[0] = adj[1] = nullptr; }
};

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
    vector<Node> nodes(path.size());
    for (size_t i = 0; i < path.size(); i++)
        nodes[i].p = path[i];
    for (size_t i = 0; i + 1 < nodes.size(); i++)
        nodes[i].adj[1] = &nodes[i + 1], nodes[i + 1].adj[0] = &nodes[i];

    queue<pair<Node *, bool>> q;
    for (size_t i = 0; i + 1 < nodes.size(); i++)
        q.emplace(&nodes[i], 0), q.emplace(&nodes[i + 1], 1);

    size_t iteration_count = 0;
    while (!q.empty() && iteration_count < TwoOptIterationLimit)
    {
        auto const [w, direction] = q.front();
        q.pop();
        if (!w->adj[!direction])
            continue;

        iteration_count++;
        Node *v = w->adj[!direction], *u = v->adj[!direction],
             *x = w->adj[direction], *a = w, *b = x;

        while (b && b->adj[direction])
        {
            Node *c = b->adj[direction], *d = c->adj[direction];

            // Überprüfe, ob die 4 neuen Abbiegewinkel (uvb, vba, xwc, wcd) alle
            // <= pi / 2 sind.
            if ((!u || dot_product(v->p - u->p, b->p - v->p) >= 0) &&
                dot_product(b->p - v->p, a->p - b->p) >= 0 &&
                dot_product(w->p - x->p, c->p - w->p) >= 0 &&
                (!d || dot_product(c->p - w->p, d->p - c->p) >= 0) &&
                abs(v->p - w->p) + abs(b->p - c->p) >
                    abs(v->p - b->p) + abs(w->p - c->p))
            {
                Node *z = x;   // Vertausche die Nachbarn aller Knoten von x
                while (z != b) // bis a.
                {
                    swap(z->adj[0], z->adj[1]);
                    z = z->adj[!direction];
                }
                v->adj[direction] = b; // Entferne die Kanten {v, w}, {b, c}
                b->adj[direction] = a; // und füge {v, b}, {w, c} ein.
                b->adj[!direction] = v;
                w->adj[!direction] = x;
                w->adj[direction] = c;
                c->adj[!direction] = w;
                q.emplace(v, !direction); // Die neu eingefügten Kanten können
                q.emplace(b, direction);  // erneut mit anderen vertauscht
                q.emplace(w, !direction); // werden.
                q.emplace(c, direction);
                break;
            }
            a = a->adj[direction]; // Gehe zur nächsten Kante im Pfad.
            b = b->adj[direction];
        }
    }

    Node *x = nullptr;
    bool direction;
    for (size_t i = 0; i < nodes.size(); i++)     // Find den Anfang des Pfads
        if (!nodes[i].adj[0] || !nodes[i].adj[1]) // (Knoten mit Grad 1).
        {
            x = &nodes[i];
            direction = nodes[i].adj[1]; // Lege die Laufrichtung fest.
            break;
        }
    path.clear();
    while (x) // Schreibe die neue Knotenfolge in path und berechne die neue
    {         // Länge.
        path.push_back(x->p);
        x = x->adj[direction];
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
        size_t u = path.back(), v = path.size() >= 2 ? *++path.rbegin() : SIZE_MAX,
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
    optimize_path(path);

    cout << setprecision(6) << fixed << "Tourlänge: " << path_length(path) << '\n';
    for (complex<double> const &u : path)
        cout << u.real() << ' ' << u.imag() << '\n';
}