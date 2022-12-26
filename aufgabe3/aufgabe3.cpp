#include <iostream>
#include <vector>
#include <algorithm>

std::vector<uint8_t> y;
std::vector<size_t> factorial;

void precalc_factorial(size_t n)
{
    factorial = std::vector<size_t>(n + 1);
    factorial[0] = 1;
    for (size_t i = 1; i <= n; i++)
        factorial[i] = factorial[i - 1] * i;
}

inline size_t index(std::vector<unsigned> const &p)
{
    uint32_t mask = 0;
    size_t k = 0;

    for (size_t j = 0; j < p.size(); j++)
    {
        k += (p[j] - __builtin_popcount(mask << (32 - p[j]))) *
             factorial[p.size() - j - 1];
        mask ^= 1 << p[j];
    }

    return k;
}

inline size_t gamma_index(std::vector<unsigned> const &p, size_t i)
{
    uint32_t mask = 0;
    size_t k = 0;

    // Um die Permutation in richtiger Reihenfolge von links nach rechts
    // abzuarbeiten, wird der umgekehrte Teil umgekehrt durchgegangen.
    for (size_t j = i - 1; j < p.size(); j--)
    {
        unsigned const x = p[j] - (p[j] > p[i]);
        k += (x - __builtin_popcount(mask >> (p.size() - x - 2))) *
             factorial[p.size() - (i - j) - 1];
        mask ^= 1 << (p.size() - x - 2);
    }

    for (size_t j = i + 1; j < p.size(); j++)
    {
        unsigned const x = p[j] - (p[j] > p[i]);
        k += (x - __builtin_popcount(mask >> (p.size() - x - 2))) *
             factorial[p.size() - j - 1];
        mask ^= 1 << (p.size() - x - 2);
    }

    return k;
}

int main()
{
    size_t n;
    std::cin >> n;

    std::vector<unsigned> p(n);
    for (unsigned &x : p)
    {
        std::cin >> x;
        x--;
    }

    precalc_factorial(n);
    y = std::vector<uint8_t>(1, 0);
    std::vector<uint8_t> z;

    for (size_t i = 2; i <= std::min<size_t>(n / 2, 12); i++)
    {
        z.resize(factorial[i]);
        std::vector<unsigned> s(i);
        for (size_t j = 0; j < i; j++)
            s[j] = j;

        for (size_t j = 0; j < factorial[i] / 2; j++)
        {
            unsigned min_ops1 = i, min_ops2 = i;
            for (size_t k = 0; k < i; k++)
            {
                size_t l = gamma_index(s, k);
                min_ops1 = std::min(min_ops1, y[l] + 1U);
                min_ops2 = std::min(min_ops2, y[factorial[i - 1] - l - 1] + 1U);
            }
            z[j] = min_ops1;
            z[factorial[i] - j - 1] = min_ops2;
            std::next_permutation(s.begin(), s.end());
        }

        z[0] = 0;
        swap(y, z);
    }
}
