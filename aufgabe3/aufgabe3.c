#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

typedef struct node node;
struct node
{
    size_t i, parent;
};

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
    uint64_t mask = 0; // Bit i ist 1, wenn n - i - 1 bereits aufgetreten ist.
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
    uint64_t mask = 0;
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
    // Berechne die Ziffern von i im fakultätsbasierten Zahlensystem.
    for (size_t j = 0; j < n; j++)
    {
        p[j] = i / factorial[n - j - 1];
        i %= factorial[n - j - 1];
    }

    // Addiere zu jeder Ziffer die Anzahl kleiner oder gleicher, links gelegener
    // Ziffern.
    for (size_t j = n - 1; j; j--)
        for (size_t k = j - 1; k < n; k--)
            p[j] += (p[k] <= p[j]);
}

size_t make_unique(size_t n, node *const nodes)
{
    size_t i = 0, j = 0;
    while (++j < n)
    {
        if (nodes[i].i != nodes[j].i)
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

unsigned *reconstruct_path(
    size_t n, size_t m, size_t const *const l, node *const *const tree)
{
    size_t i = 0;
    unsigned *path = malloc((n - m) * sizeof(unsigned));

    while (i != -1)
    {
        // Suche den akutellen Knoten mit Binärsuche.
        size_t a = 0, b = l[m - 1] - 1;

        while (a < b)
        {
            size_t mid = (a + b) / 2;
            if (tree[m - 1][mid].i < i)
                a = mid + 1;
            else
                b = mid;
        }

        i = tree[m - 1][a].parent;

        // Finde die Wendeoperation, die benutzt wurde, um zum vorherigen
        // Knoten zu gelangen.
        unsigned p[m + 1];
        ith_permutation(m + 1, tree[m - 1][a].parent, p);
        for (size_t j = 0; j < m + 1; j++)
        {
            if (p_index_gamma(m + 1, p, j) == tree[m - 1][a].i)
            {
                path[n - m - 1] = j;
                break;
            }
        }

        m++;
    }

    assert(m == n + 1);
    return path;
}

bool recurse(size_t n, size_t l1, size_t *l2, node *y, node **z)
{
    *z = malloc(n * l1 * sizeof(node));
    unsigned s[n];
    *l2 = 0;

    for (size_t i = 0; i < l1; i++)
    {
        ith_permutation(n, y[i].i, s);

        for (size_t j = 0; j < n; j++)
            (*z)[(*l2)++] =
                (node){.i = p_index_gamma(n, s, j), .parent = y[i].i};
    }

    qsort(*z, *l2, sizeof(node), node_cmp);
    *l2 = make_unique(*l2, *z);
    if (!(*z)[0].i)
        return 1;

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

    node *tree[n];
    memset(tree, 0, n * sizeof(node *));
    tree[n - 1] = malloc(sizeof(node));
    tree[n - 1][0] = (node){.i = p_index(n, p), .parent = -1};
    size_t l[n];
    l[n - 1] = 1;

    size_t m = n;
    while (!recurse(m, l[m - 1], l + m - 2, tree[m - 1], tree + m - 2))
        m--;
    m--;

    printf("%zu Operationen erforderlich. "
           "Indizes, nach denen gewendet werden muss (0-indexiert):\n",
           n - m);
    unsigned *path = reconstruct_path(n, m, l, tree);
    for (size_t i = 0; i < n - m; i++)
        printf("%u ", path[i]);
    putchar('\n');

    free(factorial);
    for (size_t i = 0; i < n; i++)
        free(tree[i]);
}
