#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

// Schreibt A(p) für jede Permutation der Länge k mit Index zwischen i1 und i2
// in z. In y muss A(q) für jede Permutation q der Länge k - 1 stehen.
// u = k!, v = (k - 1)!
void update_z(
    size_t k, size_t u, size_t v, uint8_t const *const y, uint8_t *const z,
    size_t i1, size_t i2)
{
    vector<unsigned> p = ith_permutation(k, i1);

    for (size_t i = i1; i < i2; i++)
    {
        z[i] = k;
        z[u - i - 1] = k;

        for (size_t j = 0; j < k; j++)
        {
            size_t const l = ind_gamma(p, j);
            z[i] = min<uint8_t>(z[i], y[l] + 1);
            // Hier wird die Symmetrie des Pancake-Graphen ausgenutzt.
            z[u - i - 1] = min<uint8_t>(z[u - i - 1], y[v - l - 1] + 1);
        }

        next_permutation(p.begin(), p.end());
    }
}

// Findet ein maximales A(p) unter den Permutationen der Länge n mit Index
// zwischen i1 und i2. Zurückgegeben wird A(p) und der Index von p.
pair<unsigned, size_t> get_max_a(
    size_t n, uint8_t const *const y, size_t i1, size_t i2)
{
    pair<unsigned, size_t> res = {0, -1};
    vector<unsigned> p = ith_permutation(n, i1);

    for (size_t i = i1; i < i2; i++)
    {
        unsigned a = n;

        for (size_t j = 0; j < n; j++)
            a = min<unsigned>(a, y[ind_gamma(p, j)] + 1);

        res = max(res, {a, i});
        next_permutation(p.begin(), p.end());
    }

    return res;
}

int main()
{
    size_t n;
    cin >> n;

    // Arrays zum Speichern von A(p) aller Permutationen einer bestimmten Länge.
    // y[i] = A(p), wenn ind(p) = i. In z werden dann die A(p) der
    // nächstgrößeren Permutationen hineingeschrieben.
    uint8_t *y = (uint8_t *)malloc(sizeof(uint8_t)),
            *z = (uint8_t *)malloc(sizeof(uint8_t));
    y[0] = 0;
    size_t u = 1, v = 1; // In der for-Schleife gilt u = k!, v = (k - 1)!
    unsigned const n_threads = thread::hardware_concurrency();

    for (size_t k = 2; k < n; k++)
    {
        u *= k;
        z = (uint8_t *)realloc(z, u * sizeof(uint8_t));
        vector<thread> threads;

        // Die k! Permutationen werden gleichmäßig unter den n_threads Threads
        // aufgeteilt.
        for (unsigned i = 0; i < n_threads; i++)
            threads.emplace_back(update_z, k, u, v, y, z,
                                 (u / 2) * i / n_threads,
                                 (u / 2) * (i + 1) / n_threads);

        for (thread &t : threads) // Warte, bis alle Threads beendet sind.
            t.join();

        z[0] = 0;
        swap(y, z);
        v *= k;
    }

    u *= n;
    pair<unsigned, size_t> res = {0, -1};
    vector<future<pair<unsigned, size_t>>> fut;

    // Für n ist es nicht mehr nötig, A(p) zu speichern, daher wird get_max_a
    // verwendet.
    for (unsigned i = 0; i < n_threads; i++)
        fut.emplace_back(
            async(get_max_a, n, y, u * i / n_threads, u * (i + 1) / n_threads));

    for (auto &f : fut)
        res = max(res, f.get());

    cout << "P(" << n << ") = " << res.first << '\n';
    cout << "Beispiel mit A(p) = P(" << n << "): ";
    for (unsigned const x : ith_permutation(n, res.second))
        cout << x + 1 << ' ';
    cout << '\n';
}