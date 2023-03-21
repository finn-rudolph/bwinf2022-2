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

// Optimiert den gegebenen Hamiltonpfad in-place mit der 2-opt Heuristik, ohne
// Abbiegewinkel > pi / 2 zu erzeugen. Gibt die neue Länge zurück.
double optimize_path(vector<complex<double>> &path)
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
    double new_length = 0.0;
    while (x) // Schreibe die neue Knotenfolge in path und berechne die neue
    {         // Länge.
        path.push_back(x->p);
        x = x->adj[direction];
        if (x)
            new_length += abs(x->p - path.back());
    }
    return new_length;
}

pair<vector<complex<double>>, double> obtuse_path(
    vector<complex<double>> const &z)
{
    size_t const n = z.size();
    deque<size_t> path;
    vector<bool> visited(n);
    bool front_is_dead_end, back_is_dead_end, extending_back;
    size_t no_added_node;

    auto restart_search = [&]()
    {
        front_is_dead_end = back_is_dead_end = extending_back = 0;
        fill(visited.begin(), visited.end(), 0);
        path.clear();
        srand(time(0));
        path.push_back(rand() % n);
        visited[path.front()] = 1;
        no_added_node = 0;
    };

    restart_search();

    while (path.size() < n)
    {
        if (no_added_node > n / 2) // Zu große Anzahl aufeinanderfolgener
            restart_search();      // Iterationen ohne Hinzufügen eines Knotens.

        size_t const u = extending_back ? path.back() : path.front();

        vector<size_t> successors;
        for (size_t v = 0; v < n; v++)
            if (!visited[v])
            {
                bool can_extend = path.size() < 2;
                if (!can_extend) // Überprüfe, ob eine Erweiterung des Pfads mit
                {                // v zu einem Abbiegewinkel <= pi / 2 führt.
                    size_t const pre = extending_back ? *++path.crbegin()
                                                      : *++path.cbegin();
                    can_extend = dot_product(z[u] - z[pre], z[v] - z[u]) >= 0;
                }
                if (can_extend)
                    successors.push_back(v);
            }
        if (!successors.empty()) // Wähle einen zufälligen Knoten aus, um den
        {                        // Pfad zu erweitern.
            size_t const v = successors[rand() % successors.size()];
            if (extending_back)
                path.push_back(v);
            else
                path.push_front(v);
            visited[v] = 1;
            no_added_node = 0;
            continue;
        }

        // Der Pfad kann von u aus nicht mehr erweitert werden, da alle mit
        // Abbiegewinkel <= pi / 2 schon besucht wurden. Das Ende des Pfads mit
        // u wird als Sackgasse markiert.
        (extending_back ? back_is_dead_end : front_is_dead_end) = 1;
        no_added_node++;

        if (front_is_dead_end && back_is_dead_end)
        {
            // Laufe den Pfad vom Ende zurück und finde alle Suffixe, die
            // umgekehrt werden können.
            vector<deque<size_t>::reverse_iterator> rev_points;
            auto it = path.rbegin() + 2;
            size_t const last = *path.rbegin(), scn_last = *++path.rbegin();
            while (it != path.rend())
            {
                if (dot_product(z[last] - z[scn_last], z[*it] - z[last]) >= 0 &&
                    (it + 1 == path.rend() ||
                     dot_product(z[*it] - z[last], z[*(it + 1)] - z[*it]) >= 0))
                {
                    rev_points.push_back(it);
                }
                it++;
            }
            if (!rev_points.empty()) // Drehe ein zufälliges Suffix des Pfades
            {                        // um (von den möglichen Suffixen).
                reverse(path.rbegin(), rev_points[rand() % rev_points.size()]);
                back_is_dead_end = 0;
            }
            else // Kein Suffix kann umgekehrt werden -> Versuche ein Präfix.
                reverse(path.begin(), path.end());
        }
        else // Versuche den Pfad am anderen Ende zu erweitern.
            extending_back = !extending_back;
    }

    vector<complex<double>> point_order = {z[path[0]]};
    for (size_t i = 1; i < n; i++)
        point_order.push_back(z[path[i]]);
    double length = optimize_path(point_order);
    return {point_order, length};
}

int main()
{
    vector<complex<double>> z;
    double x, y;
    while (scanf("%lf %lf", &x, &y) == 2)
        z.emplace_back(x, y);

    auto const [tour, length] = obtuse_path(z);

    cout << setprecision(6) << fixed << "Tourlänge: " << length << '\n';
    for (complex<double> const &u : tour)
        cout << u.real() << ' ' << u.imag() << '\n';
}