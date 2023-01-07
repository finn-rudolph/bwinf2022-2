#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

// Berechnet den Index von p in einer lexikographisch sortierten Liste aller
// |p|! Permutationen der Laenge |p|.
size_t ind(vector<unsigned> const &p)
{
    size_t k = 0; // Der Index von p.
    size_t const lgn = 32 - __builtin_clz(p.size()), m = 1 << lgn;
    unsigned tree[2 * m];
    memset(tree, 0, 2 * m * sizeof(unsigned));

    for (size_t j = 0; j < p.size(); j++)
    {
        size_t z = m + p[j]; // Aktueller Knoten.
        k *= (p.size() - j);
        k += p[j];
        // Angefangen beim p[j]-ten Blatt wird der Segmentbaum nach oben
        // durchlaufen.
        for (size_t l = 0; l < lgn; l++)
        {
            // Die Anzahl links gelegener, kleinerer Elemente wird abgezogen.
            if (z & 1)
                k -= tree[z - 1];
            tree[z]++;
            z >>= 1;
        }
        tree[z]++;
    }

    return k;
}

vector<unsigned> gamma(vector<unsigned> const &p, size_t i)
{
    vector<unsigned> s(p.size() - 1);

    for (size_t j = 0; j < i; j++)
        s[j] = p[i - j - 1] - (p[i - j - 1] > p[i]);
    for (size_t j = i; j < p.size() - 1; j++)
        s[j] = p[j + 1] - (p[j + 1] > p[i]);

    return s;
}

// Berechnet den Index von gamma_i p.
size_t ind_gamma(vector<unsigned> const &p, size_t i)
{
    size_t k = 0;
    size_t const lgn = 32 - __builtin_clz(p.size() - 1), m = 1 << lgn;
    unsigned tree[2 * m];
    memset(tree, 0, 2 * m * sizeof(unsigned));

    // Die Elemente vor i werden in umgekehrter Reihenfolge bearbeitet. Das j-te
    // Element gamma_i p ist das (i - j - 1)-te Element von p, fuer j < i.
    for (size_t j = 0; j < i; j++)
    {
        unsigned const x = p[i - j - 1] - (p[i - j - 1] > p[i]);
        size_t z = m + x;
        k *= (p.size() - 1 - j);
        k += x;
        for (size_t l = 0; l < lgn; l++)
        {
            if (z & 1)
                k -= tree[z - 1];
            tree[z]++;
            z >>= 1;
        }
        tree[z]++;
    }

    for (size_t j = i; j < p.size() - 1; j++)
    {
        unsigned const x = p[j + 1] - (p[j + 1] > p[i]);
        size_t z = m + x;
        k *= (p.size() - 1 - j);
        k += x;
        for (size_t l = 0; l < lgn; l++)
        {
            if (z & 1)
                k -= tree[z - 1];
            tree[z]++;
            z >>= 1;
        }
        tree[z]++;
    }

    return k;
}

// Schreibt die Ziffern von i im fakultaetsbasierten Zahlensystem in s.
void calc_factorial_digits(size_t i, vector<unsigned> &s)
{
    for (size_t j = 1; j <= s.size(); j++)
    {
        s[s.size() - j] = i % j;
        i /= j;
    }
}

vector<unsigned> ith_permutation(size_t n, size_t i)
{
    vector<unsigned> p(n);
    calc_factorial_digits(i, p);

    size_t const lgn = 32 - __builtin_clz(p.size()), m = 1 << lgn;
    unsigned tree[2 * m];
    for (size_t l = 0; l <= lgn; l++) // Initialisiere den Baum mit Einsen.
        for (size_t j = 0; j < (1ULL << l); j++)
            tree[(1 << l) + j] = 1 << (lgn - l);

    for (size_t j = 0; j < n; j++)
    {
        size_t z = 1;
        for (size_t l = 0; l < lgn; l++)
        {
            tree[z]--;
            z <<= 1;
            // Wenn nach rechts gegangen wird, muss die Anahl benoetigter,
            // kleinerer Elemente entsprechend verringert werden.
            if (p[j] >= tree[z])
                p[j] -= tree[z++];
        }
        tree[z] = 0;
        p[j] = z - m;
    }

    return p;
}

// Die eigentliche Wende-und-Ess-Operation.
vector<unsigned> rev_and_eat(vector<unsigned> const &p, size_t i)
{
    vector<unsigned> s(p.size() - 1);

    copy(p.begin(), p.begin() + i, s.begin());
    reverse(s.begin(), s.begin() + i);
    copy(p.begin() + i + 1, p.end(), s.begin() + i);

    return s;
}
