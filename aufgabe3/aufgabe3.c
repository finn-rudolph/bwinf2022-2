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

size_t opt_precomputation_size(size_t n)
{
    size_t a = 1, b = min(n - 1, 12);
    while (a < b)
    {
        size_t mid = (a + b + 1) / 2;
        if (factorial[mid] > factorial[n] / factorial[mid])
            b = mid - 1;
        else
            a = mid;
    }
    return a;
}

void precompute_next(size_t n, uint8_t const *const y, uint8_t *const z)
{
    unsigned s[n];
    for (size_t j = 0; j < n; j++)
        s[j] = j;

    for (size_t j = 0; j < factorial[n] / 2; j++)
    {
        unsigned min_ops1 = n, min_ops2 = n;
        for (size_t k = 0; k < n; k++)
        {
            size_t l = p_index_gamma(n, s, k);
            min_ops1 = min(min_ops1, y[l] + 1U);
            min_ops2 = min(min_ops2, y[factorial[n - 1] - l - 1] + 1U);
        }
        z[j] = min_ops1;
        z[factorial[n] - j - 1] = min_ops2;
        next_permutation(n, s);
    }

    z[0] = 0;
}

typedef struct node node;
struct node
{
    size_t i, parent1, parent2;
};

size_t binary_search(size_t n, node *const nodes, size_t i)
{
    size_t a = 0, b = n;
    while (a < b)
    {
        size_t mid = (a + b) / 2;
        if (nodes[mid].i < i)
            a = mid + 1;
        else
            b = mid;
    }
    return a;
}

size_t make_unique(size_t n, node *const nodes)
{
    size_t i = 0, j = 1;
    while (j < n)
    {
        if (nodes[i].i == nodes[j].i)
            j++;
        else
            nodes[++i] = nodes[j++];
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

uint8_t *optimal_gamma_seq(
    size_t n, unsigned const *const p, uint8_t const *const *const y)
{
    node *nodes[n];
    nodes[n - 1] = malloc(sizeof(node));
    nodes[n - 1][0] = (node){.i = p_index(n, p), .parent1 = -1, .parent2 = -1};
    size_t len[n];
    len[n - 1] = 1;
    unsigned s[n];

    for (size_t m = n - 1; m < n; m--)
    {
        size_t const l = len[m - 1] = len[m] * m;
        node *z = nodes[m - 1] = malloc(l * sizeof(node));
        size_t curr_len = 0;

        for (size_t i = 0; i < len[m]; i++)
        {
            node *x = nodes[m] + i;

            if (!x->i)
            {
            }
            else if (x->parent2 != -1 && x->i == factorial[m + 1] - 1)
            {
            }

            ith_permutation(m, x->i, s);

            for (size_t j = 0; j < m + 1; j++)
            {
                z[curr_len] = (node){.i = p_index_gamma(m + 1, s, j),
                                     x->i,
                                     factorial[m + 1] - x->i - 1};
                if (x->parent2 == -1)
                    z[curr_len].parent2 = -1;
                curr_len++;
            }
        }

        qsort(z, l, sizeof(node *), node_cmp);

        // Entferne einen Knoten symmetrischer Paare.
        for (size_t i = 0; i < l; i++)
        {
            size_t const sym = binary_search(l, z, factorial[m] - z[i].i - 1);

            if (z[sym].i == factorial[m] - z[i].i - 1)
            {
                z[i].parent2 = nodes[sym]->parent1;
                for (size_t j = sym;
                     j < l && z[j].i == factorial[m] - z[i].i - 1; j++)
                    z[j].i = z[sym - 1].i;
            }
        }

        len[m - 1] = make_unique(l, z);
    }
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
    size_t opt_size = opt_precomputation_size(n);
    uint8_t *y[opt_size];
    for (size_t i = 0; i < opt_size; i++)
        y[i] = malloc(factorial[i + 1] * sizeof(uint8_t));
    y[0][0] = 0;

    for (size_t i = 2; i <= opt_precomputation_size(n); i++)
        precompute_next(i, y[i - 2], y[i - 1]);

    for (size_t i = 0; i < opt_size; i++)
        free(y[i]);
    free(factorial);
}
