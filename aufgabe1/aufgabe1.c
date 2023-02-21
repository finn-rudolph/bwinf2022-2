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

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))

bool is_acute(complex double a, complex double b, complex double c)
{
    return fma(creal(b), creal(c),
               fma(creal(a), creal(b),
                   fma(cimag(b), cimag(c), cimag(a) * cimag(b)))) -
               fma(creal(b), creal(b),
                   fma(creal(a), creal(c),
                       fma(cimag(b), cimag(b), cimag(a) * cimag(c)))) <
           0.0;
}

size_t edge_index(size_t n, size_t u, size_t v)
{
    return n * u + v + 1;
}

void ban_triple(glp_prob *ip, size_t n, size_t i, size_t j, size_t k)
{
    int ind[3];
    double val[3];
    size_t const i0 = glp_add_rows(ip, 1);
    glp_set_row_bnds(ip, i0, GLP_DB, 0, 1);
    ind[1] = edge_index(n, i, j);
    ind[2] = edge_index(n, j, k);
    val[1] = val[2] = 1;
    glp_set_mat_row(ip, i0, 2, ind, val);
}

void add_angle_constraints(
    glp_prob *ip, size_t n, complex double const *const z)
{

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i + 1; j < n; j++)
        {
            for (size_t k = j + 1; k < n; k++)
            {
                if (is_acute(z[i], z[j], z[k]) && is_acute(z[k], z[i], z[j]) &&
                    is_acute(z[j], z[k], z[i]))
                {
                    // Im Dreieck ijk gibt es keinen Winkel >= pi / 2, es kann
                    // also maximal eine Kante daraus vorkommen.
                    int ind[7];
                    double val[7];
                    size_t const i0 = glp_add_rows(ip, 1);
                    glp_set_row_bnds(ip, i0, GLP_DB, 0, 1);
                    ind[1] = edge_index(n, i, j);
                    ind[2] = edge_index(n, j, i);
                    ind[3] = edge_index(n, j, k);
                    ind[4] = edge_index(n, k, j);
                    ind[5] = edge_index(n, k, i);
                    ind[6] = edge_index(n, i, k);
                    val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = 1;
                    glp_set_mat_row(ip, i0, 6, ind, val);
                }
                else
                {
                    // Füge die Beschränkung hinzu, dass von jedem Kantenpaar
                    // mit Innenwinkel < pi / 2 maximal eine Kante vorkommt.
                    if (is_acute(z[i], z[j], z[k]))
                        ban_triple(ip, n, i, j, k), ban_triple(ip, n, k, j, i);
                    if (is_acute(z[k], z[i], z[j]))
                        ban_triple(ip, n, k, i, j), ban_triple(ip, n, j, i, k);
                    if (is_acute(z[j], z[k], z[i]))
                        ban_triple(ip, n, j, k, i), ban_triple(ip, n, i, k, j);
                }
            }
        }
    }
}

void add_degree_constraints(glp_prob *ip, size_t n)
{
    size_t const i0 = glp_add_rows(ip, 3 * n);
    for (size_t i = 0; i < 2 * n; i++)
        glp_set_row_bnds(ip, i0 + i, GLP_DB, 0, 1);

    int *ind = malloc((n + 1) * sizeof *ind);
    double *val = malloc((n + 1) * sizeof *val);

    for (size_t i = 1; i < n + 1; i++)
        val[i] = 1;

    for (size_t i = 0; i < n; i++)
    {
        // Schreibe den Index von x_ij für j != i in ind, sodass der Koeffizient
        // jeder anliegenden Kante 1 ist.
        for (size_t j = 0; j < n; j++)
            ind[j + 1] = edge_index(n, i, j);
        glp_set_mat_row(ip, i0 + i, n, ind, val);
        for (size_t j = 0; j < n; j++)
            ind[j + 1] = edge_index(n, j, i);
        glp_set_mat_row(ip, i0 + n + i, n, ind, val);
    }

    for (size_t i = 0; i < n; i++)
    {
        ind[1] = edge_index(n, i, i);
        glp_set_row_bnds(ip, i0 + 2 * n + i, GLP_FX, 0, 0);
        glp_set_mat_row(ip, i0 + 2 * n + i, 1, ind, val);
    }

    free(ind);
    free(val);
}

void add_connectivity_constraint(glp_prob *ip, size_t n)
{
    size_t const i0 = glp_add_rows(ip, 1);
    glp_set_row_bnds(ip, i0, GLP_FX, n - 1, n - 1);
    int *ind = malloc((n * n + 1) * sizeof *ind);
    double *val = malloc((n * n + 1) * sizeof *val);

    // Summiere die Werte der Variablen aller Kanten.
    for (size_t i = 1; i < n * n + 1; i++)
        ind[i] = i, val[i] = 1;

    glp_set_mat_row(ip, i0, n * n, ind, val);
    free(ind);
    free(val);
}

void add_subtour_elimination_constraints(glp_prob *ip, size_t n)
{
    size_t j0 = glp_add_cols(ip, n * n),
           i0 = glp_add_rows(ip, n * n);

    for (size_t i = 0; i < n * n; i++)
    {
        glp_set_col_kind(ip, j0 + i, GLP_CV);
        glp_set_col_bnds(ip, j0 + i, GLP_DB, 0, n - 1);
    }

    int *ind = malloc(2 * n * sizeof *ind);
    double *val = malloc(2 * n * sizeof *val);
    val[1] = 1, val[2] = -((int64_t)n - 1);

    for (size_t i = 0; i < n * n; i++)
    {
        ind[1] = j0 + i, ind[2] = j0 + i - n * n;
        glp_set_row_bnds(ip, i0 + i, GLP_DB, -((int64_t)n - 1), 0);
        glp_set_mat_row(ip, i0 + i, 2, ind, val);
    }

    i0 = glp_add_rows(ip, n - 1);
    for (size_t i = 1; i < n + 1; i++)
        val[i] = 1;
    for (size_t i = n + 1; i < 2 * n; i++)
        val[i] = -1;

    for (size_t i = 1; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
            if (i != j)
                ind[j + (j < i)] = j0 + edge_index(n, j, i);
        for (size_t j = 1; j < n; j++)
            if (i != j)
                ind[n + j - 1 - (j > i)] = j0 + edge_index(n, i, j);
        glp_set_mat_row(ip, i0 + i - 1, 2 * n - 3, ind, val);
    }

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
    glp_add_cols(ip, n * n);

    // Setze alle Variablen binär.
    for (size_t i = 1; i < n * n + 1; i++)
        glp_set_col_kind(ip, i, GLP_BV);

    // Die Länge jeder verwendeten Kante wird zur Kostenfunktion addiert.
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            glp_set_obj_coef(ip, edge_index(n, i, j), cabs(z[i] - z[j]));

    add_angle_constraints(ip, n, z);
    add_degree_constraints(ip, n);
    add_connectivity_constraint(ip, n);
    add_subtour_elimination_constraints(ip, n);

    glp_iocp parameters;
    glp_init_iocp(&parameters);
    parameters.presolve = GLP_ON;
    parameters.ps_heur = GLP_ON;

    glp_intopt(ip, &parameters);

    if (glp_mip_status(ip) == GLP_NOFEAS)
    {
        printf("Keine Tour mit Abbiegewinkeln von höchstens 90 Grad möglich.\n");
        return 0;
    }

    printf("Gesamtlänge: %lf\n", glp_mip_obj_val(ip));

    size_t *successor = malloc(n * sizeof *successor);
    bool *has_predecessor = calloc(n, sizeof *has_predecessor);
    for (size_t i = 0; i < n; i++)
        successor[i] = SIZE_MAX;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (glp_mip_col_val(ip, edge_index(n, i, j)))
                successor[i] = j, has_predecessor[j] = 0;

    size_t u;
    for (size_t i = 0; i < n; i++)
        if (!has_predecessor[i])
            u = i;
    while (u != SIZE_MAX)
    {
        printf("%lf %lf\n", creal(z[u]), cimag(z[u]));
        u = successor[u];
    }

    glp_delete_prob(ip);
    free(z);
    free(successor);
    free(has_predecessor);
}