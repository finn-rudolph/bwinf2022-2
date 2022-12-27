#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <memory.h>

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

__uint128_t *factorial;

void precalc_factorial(size_t n)
{
    factorial = malloc((n + 1) * sizeof(__uint128_t));
    factorial[0] = 1;
    for (size_t i = 1; i <= n; i++)
        factorial[i] = factorial[i - 1] * i;
}

size_t p_index(size_t n, unsigned const *const p)
{
    uint32_t mask = 0; // Bit i ist 1, wenn n - i - 1 bereits aufgetreten ist.
    size_t k = 0;

    for (size_t j = 0; j < n; j++)
    {
        k += (p[j] - __builtin_popcount(mask >> (n - p[j] - 1))) *
             factorial[n - j - 1];
        mask ^= 1 << (n - p[j] - 1);
    }

    return k;
}

size_t p_index_gamma(size_t n, unsigned const *const p, size_t i)
{
    uint32_t mask = 0;
    size_t k = 0;

    // Um die Permutation in richtiger Reihenfolge von links nach rechts
    // abzuarbeiten, wird der umgekehrte Teil umgekehrt durchgegangen.
    for (size_t j = i - 1; j < n; j--)
    {
        unsigned const x = p[j] - (p[j] > p[i]);
        k += (x - __builtin_popcount(mask >> (n - x - 2))) *
             factorial[n - (i - j) - 1];
        mask ^= 1 << (n - x - 2);
    }

    for (size_t j = i + 1; j < n; j++)
    {
        unsigned const x = p[j] - (p[j] > p[i]);
        k += (x - __builtin_popcount(mask >> (n - x - 2))) *
             factorial[n - j - 1];
        mask ^= 1 << (n - x - 2);
    }

    return k;
}

void ith_permutation(size_t n, size_t i, unsigned *const p)
{
    for (size_t j = 0; j < n; j++)
    {
        p[i] = i / factorial[n - j - 1];
        i %= factorial[n - j - 1];
    }

    for (size_t j = n - 1; j < n; j--)
        for (size_t k = j - 1; k < n; k--)
            if (p[k] <= p[j])
                p[j]++;
}

void next_permutation(size_t n, unsigned *const p)
{
    if (n == 1)
        return;

    size_t m = n - 2, k = n - 1;
    while (m && p[m] > p[m + 1])
        m--;
    while (p[m] > p[k])
        k--;
    unsigned x = p[k];
    p[k] = p[m];
    p[m] = x;

    m++;
    k = n - 1;
    while (m < k)
    {
        unsigned x = p[k];
        p[k] = p[m];
        p[m] = x;
        m++, k--;
    }
}

typedef struct node node;
struct node
{
    size_t i, parent1, parent2;
};

size_t make_unique(size_t n, node *const nodes)
{
    if (n == 1)
        return;

    size_t i = 0, j = 0;
    while (++j < n)
    {
        if (nodes[i].i != nodes[j].i && nodes[j].i != -1)
            nodes[++i] = nodes[j];
    }
    return i + 1;
}

int node_cmp(void const *const a, void const *const b)
{
    if (((node *)a)->i > ((node *)b)->i)
        return 1;
    else if (((node *)a)->i < ((node *)b)->i)
        return -1;
    return 0;
}

unsigned *reconstruct_path(size_t n, node *z)
{
    unsigned *path = malloc(n * sizeof(unsigned));
    return path;
}

bool next_rec(size_t n, size_t l1, size_t *l2, node *y, node *z)
{
    z = malloc(n * l1 * sizeof(node));
    unsigned s[n];
    *l2 = 0;

    for (size_t i = 0; i < l1; i++)
    {
        ith_permutation(n, y[i].i, s);

        for (size_t j = 0; j < n; j++)
        {
            z[*l2] = (node){.i = p_index_gamma(n, s, j),
                            .parent1 = z[i].i,
                            .parent2 = factorial[n] - y[i].i - 1};
            if (y[i].parent2 == -1)
                z[*l2].parent2 = -1;
            (*l2)++;
        }
    }

    qsort(z, *l2, sizeof(node), node_cmp);
    if (!z[0].i)
        return 1;

    *l2 = make_unique(l2, z);
    size_t i = 1, j = *l2 - 1;

    while (i < j)
    {
        if (z[i].i < factorial[n - 1] - z[j].i - 1)
            i++;
        else if (z[i].i > factorial[n - 1] - z[j].i - 1)
            j--;
        else
        {
            z[j].i = -1;
            j--, i++;
        }
    }

    *l2 = make_unique(l2, z);
    return 0;
}

int main()
{
    size_t n;
    scanf("%zu", &n);

    unsigned p[n];
    for (size_t i = 0; i < n; i++)
    {
        scanf("%u", p + i);
        p[i]--;
    }

    precalc_factorial(n);

    node *z[n];
    z[n - 1] = malloc(sizeof(node));
    z[n - 1][0] = (node){.i = p_index(n, p), .parent1 = -1, .parent2 = -1};
    size_t l[n];
    l[n - 1] = 1;

    free(factorial);
}
