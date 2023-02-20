#include <stdio.h>
#include <complex.h>
#include <tgmath.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <glpk.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define min(x, y) ((x < y) ? (x) : (y))

double dot_product(complex double const x, complex double const y)
{
    return fma(creal(x), creal(y), cimag(x) * cimag(y));
}

size_t nchoose2(size_t const n)
{
    return n * (n - 1) / 2;
}

size_t edge_index(size_t n, size_t u, size_t v)
{
    return nchoose2(n) - nchoose2(n - min(u, v)) + max(u, v) + 1;
}

void add_angle_constraints(
    glp_prob *ip, size_t n, complex double const *const z)
{
    int ind[3];
    double val[3];
    val[1] = val[2] = 1;

    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
            for (size_t k = j + 1; k < n; k++)
                if (dot_product(z[j] - z[i], z[k] - z[j]) < 0)
                {
                    size_t const i0 = glp_add_rows(ip, 1);
                    glp_set_row_bnds(ip, i0, GLP_DB, 0, 1);
                    ind[1] = edge_index(n, i, j);
                    ind[2] = edge_index(n, j, k);
                    glp_set_mat_row(ip, i0, 2, ind, val);
                }
}

void add_degree_constraints(glp_prob *ip, size_t n)
{
    size_t const i0 = glp_add_rows(ip, n);
    for (size_t i = 0; i < n; i++)
        glp_set_row_bnds(ip, i0 + i, GLP_DB, 1, 2);

    int *ind = malloc(n * sizeof *ind);
    double *val = malloc(n * sizeof *val);

    for (size_t i = 1; i < n; i++)
        val[i] = 1;

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
            if (i != j)
                ind[j + 1 - (j > i)] = edge_index(n, i, j);
        glp_set_mat_row(ip, i0 + i, n - 1, ind, val);
    }

    free(ind);
    free(val);
}

void add_connectivity_constraint(glp_prob *ip, size_t n)
{
    size_t const i0 = glp_add_rows(ip, 1);
    glp_set_row_bnds(ip, i0, GLP_FX, n - 1, n - 1);
    int *ind = malloc((nchoose2(n) + 1) * sizeof *ind);
    double *val = malloc((nchoose2(n) + 1) * sizeof *val);

    for (size_t i = 1; i < nchoose2(n) + 1; i++)
        ind[i] = i, val[i] = 1;

    glp_set_mat_row(ip, i0, nchoose2(n), ind, val);
    free(ind);
    free(val);
}

void add_subtour_elimination_constraint(
    glp_prob *ip, size_t n, size_t m, size_t const *const subtour)
{
    size_t const i0 = glp_add_rows(ip, 1);
    glp_set_row_bnds(ip, i0, GLP_DB, 0, m - 1);

    int *ind = malloc((nchoose2(m) + 1) * sizeof *ind);
    double *val = malloc((nchoose2(m) + 1) * sizeof *val);

    for (size_t i = 1; i < nchoose2(m) + 1; i++)
        val[i] = 1;

    for (size_t i = 0; i < m; i++)
        for (size_t j = i + 1; j < m; j++)
            ind[edge_index(m, i, j)] = edge_index(n, subtour[i], subtour[j]);

    glp_set_mat_row(ip, i0, nchoose2(m), ind, val);
    free(ind);
    free(val);
}

int main()
{
    complex double *z = malloc(sizeof *z);
    size_t n = 0, capacity = 1;
    double x, y;
    while (scanf("%lf %lf", &x, &y) == 2)
    {
        if (n == capacity)
            capacity <<= 1, z = realloc(z, capacity * sizeof *z);
        z[n++] = x + y * I;
    }

    glp_prob *ip = glp_create_prob();
    glp_set_obj_dir(ip, GLP_MIN);
    glp_add_cols(ip, nchoose2(n));

    add_angle_constraints(ip, n, z);
    add_degree_constraints(ip, n);
    add_connectivity_constraint(ip, n);
    for (size_t i = 1; i < nchoose2(n) + 1; i++)
        glp_set_col_kind(ip, i, GLP_BV);

    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
            glp_set_obj_coef(ip, edge_index(n, i, j), cabs(z[i] - z[j]));

    glp_iocp parameters;
    glp_init_iocp(&parameters);
    parameters.presolve = GLP_ON;
    parameters.ps_heur = GLP_ON;

    bool solution_found = 0;
    bool *visited = malloc(n * sizeof *visited);
    size_t *successor = malloc(n * sizeof *successor),
           *cycle = malloc(n * sizeof *cycle);

    do
    {
        glp_intopt(ip, &parameters);

        memset(visited, 0, n * sizeof *visited);
        for (size_t i = 0; i < n; i++)
            successor[i] = SIZE_MAX;

        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                if (glp_mip_col_val(ip, i * n + j + 1) > 0.1)
                    successor[i] = j;

        solution_found = 1;
        for (size_t i = 0; i < n; i++)
        {
            if (!visited[i])
            {
                size_t cycle_len = 0, u = i;
                do
                {
                    cycle[cycle_len++] = u;
                    visited[u] = 1;
                    u = successor[u];
                } while (u != SIZE_MAX && u != i);

                if (u == i)
                {
                    add_subtour_elimination_constraint(ip, n, cycle_len, cycle);
                    solution_found = 0;
                }
            }
        }

    } while (!solution_found && glp_mip_status(ip) != GLP_NOFEAS);

    if (glp_mip_status(ip) == GLP_NOFEAS)
    {
        printf("Keine Tour mit Abbiegewinkeln kleiner gleich 90 Grad möglich.\n");
        return 0;
    }

    printf("Gesamtlänge: %lf\n", glp_mip_obj_val(ip));

    size_t u;
    for (size_t i = 0; i < n; i++)
        if (glp_mip_row_val(ip, 3 * n + i + 1) < 0.1)
            u = i;

    for (size_t i = 0; i < n; i++)
        printf("%zu: %zu\n", i, successor[i]);

    while (u != SIZE_MAX)
    {
        printf("%lf %lf\n", creal(z[u]), cimag(z[u]));
        u = successor[u];
    }

    glp_delete_prob(ip);
    free(z);
    free(visited);
    free(successor);
    free(cycle);
}