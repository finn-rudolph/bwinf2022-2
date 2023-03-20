#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

struct Node
{
    Node *adj[2];
    complex<double> p;

    Node() { adj[0] = adj[1] = nullptr; }
};

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

    while (!q.empty())
    {
        auto const [w, direction] = q.front();
        q.pop();
        if (!w->adj[!direction])
            continue;
        Node *v = w->adj[!direction], *u = v->adj[!direction],
             *x = w->adj[direction], *a = w, *b = x;

        while (b && b->adj[direction])
        {
            Node *c = b->adj[direction], *d = c->adj[direction];

            if ((!u || dot_product(v->p - u->p, b->p - v->p) >= 0) &&
                dot_product(b->p - v->p, a->p - b->p) >= 0 &&
                dot_product(w->p - x->p, c->p - w->p) >= 0 &&
                (!d || dot_product(c->p - w->p, d->p - c->p) >= 0) &&
                abs(v->p - w->p) + abs(b->p - c->p) >
                    abs(v->p - b->p) + abs(w->p - c->p))
            {
                Node *z = x;
                while (z != b)
                {
                    swap(z->adj[0], z->adj[1]);
                    z = z->adj[!direction];
                }
                v->adj[direction] = b;
                b->adj[direction] = a;
                b->adj[!direction] = v;
                w->adj[!direction] = x;
                w->adj[direction] = c;
                c->adj[!direction] = w;
                q.emplace(v, !direction);
                q.emplace(b, direction);
                q.emplace(w, !direction);
                q.emplace(c, direction);
                break;
            }
            a = a->adj[direction];
            b = b->adj[direction];
        }
    }

    Node *x = nullptr;
    bool direction;
    for (size_t i = 0; i < nodes.size(); i++)
        if (!nodes[i].adj[0] || !nodes[i].adj[1])
        {
            x = &nodes[i];
            direction = nodes[i].adj[1];
            break;
        }
    path.clear();
    double new_length = 0.0;
    while (x)
    {
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
        if (no_added_node > n / 2)
            restart_search();

        size_t const u = extending_back ? path.back() : path.front();

        vector<size_t> successors;
        for (size_t v = 0; v < n; v++)
            if (!visited[v])
            {
                bool can_extend = path.size() < 2;
                if (!can_extend)
                {
                    size_t const pre = extending_back ? *++path.crbegin()
                                                      : *++path.cbegin();
                    can_extend = dot_product(z[u] - z[pre], z[v] - z[u]) >= 0;
                }
                if (can_extend)
                    successors.push_back(v);
            }
        if (!successors.empty())
        {
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
            // Laufe den Pfad vom Ende zurück und versuche ein Suffix
            // umzukehren.
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
            if (!rev_points.empty())
            {
                reverse(path.rbegin(), rev_points[rand() % rev_points.size()]);
                back_is_dead_end = 0;
            }
            else
            {
                reverse(path.begin(), path.end());
            }
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