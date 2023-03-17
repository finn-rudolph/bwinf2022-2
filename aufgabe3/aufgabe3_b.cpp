#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

// Schreibt A(p) für jede Permutation der Länge k mit Index im Intervall
// [i1, i2) in z. In y muss A(q) für jede Permutation q der Laenge k - 1 stehen.
// Des Weiteren gilt u = k! und v = (k - 1)!.
void update_z(
    size_t k, size_t u, size_t v, uint8_t const *const y, uint8_t *const z,
    size_t i1, size_t i2)
{
    vector<unsigned> p = ith_permutation(k, i1);

    for (size_t i = i1; i < i2; i++)
    {
        z[i] = k;         // Durch die Symmetrie des Pancake-Graphen kann ein
        z[u - i - 1] = k; // symmetrisches Paar gleichzeitig behandelt werden.

        for (size_t j = 0; j < k; j++)
        {
            size_t const l = ind_gamma(p, j);
            z[i] = min<uint8_t>(z[i], y[l] + 1);
            z[u - i - 1] = min<uint8_t>(z[u - i - 1], y[v - l - 1] + 1);
        }

        next_permutation(p.begin(), p.end());
    }
}

// Findet ein maximales A(p) unter den Permutationen der Länge n mit Index
// zwischen i1 und i2. Zurückgegeben wird A(p) und der Index von p.
pair<unsigned, size_t> get_max_a(
    size_t n, size_t u, size_t v, uint8_t const *const y, size_t i1, size_t i2)
{
    pair<unsigned, size_t> res = {0, -1};
    vector<unsigned> p = ith_permutation(n, i1);

    for (size_t i = i1; i < i2; i++)
    {
        unsigned a1 = n, a2 = n;

        for (size_t j = 0; j < n; j++)
        {
            size_t const l = ind_gamma(p, j);
            a1 = min<unsigned>(a1, y[l] + 1);
            a2 = min<unsigned>(a2, y[v - l - 1] + 1);
        }

        res = max(res, max(make_pair(a1, i), make_pair(a2, u - i - 1)));
        next_permutation(p.begin(), p.end());
    }

    return res;
}

pair<unsigned, size_t> pwue(size_t n)
{
    // Arrays zum Speichern von A(p) aller Permutationen einer bestimmten Länge.
    // In der folgenden for-Schleife gilt y[i] = A(p), wenn ind(p) = i, für jede
    // Permutation p der Länge k - 1. In z wird dann A(p) für alle Permutationen
    // von Länge k geschrieben.
    uint8_t *y = (uint8_t *)malloc(sizeof *y),
            *z = (uint8_t *)malloc(sizeof *z);
    y[0] = 0;
    size_t u = 1, v = 1; // In der for-Schleife gilt u = k!, v = (k - 1)!
    unsigned const n_threads = thread::hardware_concurrency();

    for (size_t k = 2; k < n; k++)
    {
        u *= k;
        z = (uint8_t *)realloc(z, u * sizeof *z);
        vector<thread> threads;

        // Die k! Permutationen werden ausgeglichen unter den n_threads Threads
        // aufgeteilt.
        for (unsigned i = 0; i < n_threads; i++)
            threads.emplace_back(update_z, k, u, v, y, z,
                                 (u / 2) * i / n_threads,
                                 (u / 2) * (i + 1) / n_threads);

        for (thread &t : threads)
            t.join();

        z[0] = 0;
        swap(y, z);
        v *= k;
    }

    free(z);
    u *= n;
    pair<unsigned, size_t> res = {0, -1};
    vector<future<pair<unsigned, size_t>>> fut;

    // Für Länge n ist es nicht mehr nötig, A(p) zu speichern, daher wird
    // get_max_a verwendet.
    for (unsigned i = 0; i < n_threads; i++)
        fut.emplace_back(async(get_max_a, n, u, v, y, (u / 2) * i / n_threads,
                               (u / 2) * (i + 1) / n_threads));

    for (auto &f : fut)
        res = max(res, f.get());
    free(y);

    return res;
}

int main()
{
    size_t n;
    cin >> n;

    pair<unsigned, size_t> const res = pwue(n);

    cout << "P(" << n << ") = " << res.first << '\n';
    cout << "Beispiel mit A(p) = P(" << n << "): ";
    for (unsigned const x : ith_permutation(n, res.second))
        cout << x + 1 << ' ';
    cout << '\n';
}