#include <bits/stdc++.h>
#include "util.hpp"
using namespace std;

uint64_t factorial(uint64_t n) { return !n ? 1 : n * factorial(n - 1); }

// Berechnet den Index von p in einer lexikographisch sortierten Liste aller
// |p|! Permutationen der Laenge |p|.
uint64_t ind(vector<unsigned> const &p)
{
    uint64_t k = 0; // Der Index von p.
    size_t const lgn = countl_zero<size_t>(0) - countl_zero(p.size()),
                 m = 1U << lgn;
    unsigned tree[2 * m]; // Summensegmentbaum
    memset(tree, 0, sizeof tree);

    for (size_t j = 0; j < p.size(); j++)
    {
        size_t z = m + p[j]; // Aktueller Knoten.
        k *= (p.size() - j);
        k += p[j];
        // Angefangen beim p[j]-ten Blatt wird der Segmentbaum nach oben
        // durchlaufen.
        for (size_t l = 0; l < lgn; l++)
        {
            if (z & 1)            // Wenn bei einem rechten Nachfolger, ziehe die
                k -= tree[z - 1]; // Zahl links gelegener, kleinerer Elemente ab.
            tree[z]++;
            z >>= 1; // Gehe zum Elternknoten.
        }
        tree[z]++;
    }

    return k;
}

vector<unsigned> gamma(vector<unsigned> const &p, unsigned i)
{
    vector<unsigned> s(p.size() - 1);

    for (size_t j = 0; j < i; j++) // Verringere Elemente > p[i] um 1.
        s[j] = p[i - j - 1] - (p[i - j - 1] > p[i]);
    for (size_t j = i; j < p.size() - 1; j++)
        s[j] = p[j + 1] - (p[j + 1] > p[i]);

    return s;
}

// Berechnet mu(gamma_i p).
uint64_t ind_gamma(vector<unsigned> const &p, unsigned i)
{
    uint64_t k = 0;
    size_t const lgn = countl_zero<size_t>(0) - countl_zero(p.size() - 1),
                 m = 1U << lgn;
    unsigned tree[2 * m]; // Summensegmentbaum
    memset(tree, 0, sizeof tree);

    // Die Elemente vor i werden in umgekehrter Reihenfolge bearbeitet. Das j-te
    // Element gamma_i p ist das (i - j - 1)-te Element von p, für j < i.
    // Daneben wird jedes Element > p[i] von gamma_i um 1 verringert.
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

    for (size_t j = i; j < p.size() - 1; j++) // Elemente nach i.
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

// Schreibt die Ziffern von i im fakultätsbasierten Zahlensystem in digits.
void calc_factorial_digits(uint64_t i, vector<unsigned> &digits)
{
    for (size_t j = 1; j <= digits.size(); j++)
    {
        digits[digits.size() - j] = i % j;
        i /= j;
    }
}

vector<unsigned> ith_permutation(unsigned n, uint64_t i)
{
    vector<unsigned> p(n);       // Verwende p zunächst als Speicher für die
    calc_factorial_digits(i, p); // Ziffern im fankultätsbasierten Zahlensystem.

    size_t const lgn = countl_zero<size_t>(0) - countl_zero(p.size()),
                 m = 1 << lgn;
    unsigned tree[2 * m];
    for (size_t l = 0; l <= lgn; l++) // Initialisiere den Baum mit Einsen.
        for (size_t j = 0; j < (1U << l); j++)
            tree[(1 << l) + j] = 1 << (lgn - l);

    for (size_t j = 0; j < n; j++)
    {
        size_t z = 1; // Index des aktuellen Knotens im Segmentbaum
        for (size_t l = 0; l < lgn; l++)
        {
            tree[z]--; // p[j] ist nun vorhanden -> setzte seinen Wert auf 0.
            z <<= 1;
            // Wenn nach rechts gegangen wird, muss die Anahl benötigter,
            // kleinerer Elemente um die Zahl kleinerer Elemente im linken
            // Teilbaum verringert werden.
            if (p[j] >= tree[z])
                p[j] -= tree[z++];
        }
        tree[z] = 0;
        p[j] = z - m;
    }

    return p;
}

// Die eigentliche Wende-und-Ess-Operation.
vector<unsigned> reverse_and_eat(vector<unsigned> const &p, unsigned i)
{
    vector<unsigned> s(p.size() - 1);

    copy(p.begin(), p.begin() + i, s.begin());
    reverse(s.begin(), s.begin() + i);
    copy(p.begin() + i + 1, p.end(), s.begin() + i);

    return s;
}
