#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

// Schreibt A(q) für jede Permutation der Länge k mit Index im Intervall
// [i1, i2) in z. In y muss A(p) für jede Permutation p der Laenge k - 1 stehen.
void update_z(
    size_t k, uint8_t const *const y, uint8_t *const z, size_t i1, size_t i2)
{
    vector<unsigned> p = ith_permutation(k, i1);
    size_t const u = factorial(k), v = factorial(k - 1);

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
    size_t n, uint8_t const *const y, size_t i1, size_t i2)
{
    pair<unsigned, size_t> res = {0, -1};
    vector<unsigned> p = ith_permutation(n, i1);
    size_t const u = factorial(n), v = factorial(n - 1);

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
    unsigned const num_threads = thread::hardware_concurrency();

    for (size_t k = 2; k < n; k++)
    {
        size_t const k_factorial = factorial(k);
        z = (uint8_t *)realloc(z, k_factorial * sizeof *z);
        vector<thread> threads;

        // Die k! Permutationen werden ausgeglichen unter den num_threads
        // Threads aufgeteilt.
        for (unsigned i = 0; i < num_threads; i++)
            threads.emplace_back(update_z, k, y, z,
                                 (k_factorial / 2) * i / num_threads,
                                 (k_factorial / 2) * (i + 1) / num_threads);

        for (thread &t : threads) // Warte, bis alle Threads fertig sind.
            t.join();

        z[0] = 0;
        swap(y, z);
    }

    free(z);
    pair<unsigned, size_t> res = {0, -1}; // P(n), Index der Beispielpermutation.
    vector<future<pair<unsigned, size_t>>> fut;

    // Für Länge n ist es nicht mehr nötig, A(p) zu speichern, daher wird
    // get_max_a als Threadfunktion verwendet. Die future-Objekte erlauben es
    // auf die Rückgabe der Threads zuzugreifen.
    size_t const n_factorial = factorial(n);
    for (unsigned i = 0; i < num_threads; i++)
        fut.emplace_back(async(get_max_a, n, y,
                               (n_factorial / 2) * i / num_threads,
                               (n_factorial / 2) * (i + 1) / num_threads));

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