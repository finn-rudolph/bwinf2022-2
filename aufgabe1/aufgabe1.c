#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
#include <stdint.h>
#include <glpk.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void add_angle_constraints(glp_prob *ip, size_t n, complex double const *const z)
{
    for (size_t i = 0; i < n; i++)
        glp_set_row_bnds(ip, i + 1, GLP_LO, -M_PI / 2, 0.0);
    for (size_t i = 0; i < n; i++)
        glp_set_row_bnds(ip, n + i + 1, GLP_UP, 0.0, 2.5 * M_PI);

    int *ind = malloc(2 * n * sizeof *ind);
    double *val = malloc(2 * n * sizeof *val);

    for (size_t i = 0; i < n; i++)
    {
        ind[2 * n - 1] = n * n + i + 1;
        val[2 * n - 1] = 2 * M_PI;

        size_t k = 1;
        for (size_t j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            ind[k] = i * n + j + 1;
            val[k] = carg(z[j] - z[i]) - M_PI / 2;
            ind[k + 1] = j * n + i + 1;
            val[k + 1] = -carg(z[i] - z[j]) - M_PI / 2;
            k += 2;
        }
        glp_set_mat_row(ip, i + 1, 2 * n - 1, ind, val);

        k = 1;
        for (size_t j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            val[k] = carg(z[j] - z[i]) + M_PI / 2;
            val[k + 1] = -carg(z[i] - z[j]) + M_PI / 2;
            k += 2;
        }
        glp_set_mat_row(ip, n + i + 1, 2 * n - 1, ind, val);
    }

    free(ind);
    free(val);
}

void add_degree_constraints(glp_prob *ip, size_t n)
{
    for (size_t i = 0; i < 2 * n; i++)
        glp_set_row_bnds(ip, 2 * n + i + 1, GLP_DB, 0.0, 1.0);
    for (size_t i = 0; i < n; i++)
        glp_set_row_bnds(ip, 4 * n + i + 1, GLP_DB, 1.0, 2.0);
    for (size_t i = 0; i < n; i++)
        glp_set_row_bnds(ip, 5 * n + i + 1, GLP_FX, 0.0, 0.0);

    int *ind = malloc((2 * n - 1) * sizeof *ind);
    double *val = malloc((2 * n - 1) * sizeof *val);

    for (size_t i = 0; i < 2 * n - 1; i++)
        val[i] = 1.0;

    for (size_t i = 0; i < n; i++)
    {
        size_t k = 1;
        for (size_t j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            ind[k] = i * n + j + 1;
            ind[n - 1 + k] = j * n + i + 1;
            k++;
        }
        glp_set_mat_row(ip, 2 * n + i + 1, n - 1, ind, val);
        glp_set_mat_row(ip, 3 * n + i + 1, n - 1, ind + n - 1, val + n - 1);
        glp_set_mat_row(ip, 4 * n + i + 1, 2 * n - 2, ind, val);
    }

    for (size_t i = 0; i < n; i++)
    {
        ind[1] = i * n + i + 1;
        glp_set_mat_row(ip, 5 * n + i + 1, 1, ind, val);
    }

    free(ind);
    free(val);
}

void add_connectivity_constraint(glp_prob *ip, size_t n)
{
    glp_set_row_bnds(ip, 6 * n + 1, GLP_FX, n - 1, n - 1);
    int *ind = malloc((n * n + 1) * sizeof *ind);
    double *val = malloc((n * n + 1) * sizeof *val);

    for (size_t i = 1; i <= n * n; i++)
    {
        ind[i] = i;
        val[i] = 1.0;
    }

    glp_set_mat_row(ip, 6 * n + 1, n * n, ind, val);
    free(ind);
    free(val);
}

void add_subtour_elimination_constraint(
    glp_prob *ip, size_t n, size_t m, size_t const *const subtour)
{
    glp_add_rows(ip, 1);
    glp_set_row_bnds(ip, glp_get_num_rows(ip), GLP_DB, 0.0, m - 1);

    int *ind = malloc((m * m + 1) * sizeof *ind);
    double *val = malloc((m * m + 1) * sizeof *val);

    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < m; j++)
            ind[i * m + j + 1] = subtour[i] * n + subtour[j] + 1;

    for (size_t i = 1; i <= m * m; i++)
        val[i] = 1.0;

    glp_set_mat_row(ip, glp_get_num_rows(ip), m * m, ind, val);
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
    glp_add_rows(ip, 6 * n + 1);
    glp_add_cols(ip, n * n + n);

    add_angle_constraints(ip, n, z);
    add_degree_constraints(ip, n);
    add_connectivity_constraint(ip, n);
    for (size_t i = 0; i < n * n + n; i++)
        glp_set_col_kind(ip, i + 1, GLP_BV);

    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            glp_set_obj_coef(ip, i * n + j + 1, cabs(z[i] - z[j]));

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