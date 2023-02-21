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

size_t nchoose2(size_t const n)
{
    return n * (n - 1) / 2;
}

size_t edge_index(size_t n, size_t u, size_t v)
{
    return nchoose2(n) - nchoose2(n - min(u, v)) + max(u, v) - min(u, v);
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
                    int ind[4];
                    double val[4];
                    size_t const i0 = glp_add_rows(ip, 1);
                    glp_set_row_bnds(ip, i0, GLP_DB, 0, 1);
                    ind[1] = edge_index(n, i, j);
                    ind[2] = edge_index(n, j, k);
                    ind[3] = edge_index(n, k, i);
                    val[1] = val[2] = val[3] = 1;
                    glp_set_mat_row(ip, i0, 3, ind, val);
                }
                else
                {
                    // Füge die Beschränkung hinzu, dass von jedem Kantenpaar
                    // mit Innenwinkel < pi / 2 maximal eine Kante vorkommt.
                    if (is_acute(z[i], z[j], z[k]))
                        ban_triple(ip, n, i, j, k);
                    if (is_acute(z[k], z[i], z[j]))
                        ban_triple(ip, n, k, i, j);
                    if (is_acute(z[j], z[k], z[i]))
                        ban_triple(ip, n, j, k, i);
                }
            }
        }
    }
}

void add_degree_constraints(glp_prob *ip, size_t n)
{
    size_t const i0 = glp_add_rows(ip, n);
    for (size_t i = 0; i < n; i++) // Beschränke den Wert jeder Zeile auf [1, 2].
        glp_set_row_bnds(ip, i0 + i, GLP_DB, 1, 2);

    int *ind = malloc(n * sizeof *ind);
    double *val = malloc(n * sizeof *val);

    for (size_t i = 1; i < n; i++)
        val[i] = 1;

    for (size_t i = 0; i < n; i++)
    {
        // Schreibe den Index von x_ij für j != i in ind, sodass der Koeffizient
        // jeder anliegenden Kante 1 ist.
        for (size_t j = 0; j < n; j++)
            if (i != j)
                ind[j + (j < i)] = edge_index(n, i, j);
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

    // Summiere die Werte der Variablen aller Kanten.
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

    // Für jedes unterschiedliche Knotenpaar in subtour, füge den Index der
    // entsprechenden Kante zu ind hinzu.
    for (size_t i = 0; i < m; i++)
        for (size_t j = i + 1; j < m; j++)
            ind[edge_index(m, i, j)] = edge_index(n, subtour[i], subtour[j]);

    glp_set_mat_row(ip, i0, nchoose2(m), ind, val);
    free(ind);
    free(val);
}

// Erstellt die Adjazenzliste mit den von GLPK als Lösung vorgeschlagenen
// Kanten.
void build_graph(glp_prob *ip, size_t n, size_t *const *const graph)
{
    // Da der Grad jedes Knoten nur 1 oder 2 sein kann, werden die Längen der
    // Adjazenzlisten nicht explizit gespeichert, sondern durch SIZE_MAX das
    // Nichtvorhandensein eines zweiten benachbarten Knoten signalisiert.
    for (size_t i = 0; i < n; i++)
        graph[i][0] = graph[i][1] = SIZE_MAX;

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i + 1; j < n; j++)
        {
            if (glp_mip_col_val(ip, edge_index(n, i, j)) > 0.1)
            {
                graph[i][1] = graph[i][0];
                graph[i][0] = j;
                graph[j][1] = graph[j][0];
                graph[j][0] = i;
            }
        }
    }
}

// Gibt die Länge des Zyklus, der u enthält zurück, oder 0, wenn u nicht Teil
// eines Zyklus ist. Die Knoten des Zyklus werden in subtour geschrieben.
size_t find_cycle_len(
    size_t *const *const graph, size_t *const subtour, bool *const visited,
    size_t u)
{
    size_t cycle_len = 0, v = u, last = u;

    do
    {
        subtour[cycle_len++] = v;
        visited[v] = 1;
        size_t const next = graph[v][0] == last ? graph[v][1] : graph[v][0];
        last = v;
        v = next;
    } while (v != SIZE_MAX && v != u);

    return u == v ? cycle_len : 0;
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

    // Setze alle Variablen binär.
    for (size_t i = 1; i < nchoose2(n) + 1; i++)
        glp_set_col_kind(ip, i, GLP_BV);

    // Die Länge jeder verwendeten Kante wird zur Kostenfunktion addiert.
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
            glp_set_obj_coef(ip, edge_index(n, i, j), cabs(z[i] - z[j]));

    glp_iocp parameters;
    glp_init_iocp(&parameters);
    parameters.presolve = GLP_ON;
    parameters.ps_heur = GLP_ON;

    bool solution_found = 0;
    bool *visited = malloc(n * sizeof *visited);
    size_t **graph = malloc(n * sizeof *graph),
           *subtour = malloc(n * sizeof *subtour);
    for (size_t i = 0; i < n; i++)
        graph[i] = malloc(2 * sizeof *graph[i]);

    do
    {
        glp_intopt(ip, &parameters);

        build_graph(ip, n, graph);
        memset(visited, 0, n * sizeof *visited);
        solution_found = 1;

        for (size_t i = 0; i < n; i++)
        {
            if (!visited[i])
            {
                // Besuche jeden Knoten des Graphen und eliminiere ggf. Zyklen.
                size_t const m = find_cycle_len(graph, subtour, visited, i);
                solution_found = solution_found && !m;
                if (m)
                    add_subtour_elimination_constraint(ip, n, m, subtour);
            }
        }

    } while (!solution_found && glp_mip_status(ip) != GLP_NOFEAS);

    if (glp_mip_status(ip) == GLP_NOFEAS)
    {
        printf("Keine Tour mit Abbiegewinkeln von höchstens 90 Grad möglich.\n");
        return 0;
    }

    printf("Gesamtlänge: %lf\n", glp_mip_obj_val(ip));

    size_t u;
    for (size_t i = 0; i < n; i++)
        if (graph[i][1] == SIZE_MAX)
            u = i;
    memset(visited, 0, n * sizeof *visited);
    while (u != SIZE_MAX)
    {
        printf("%lf %lf\n", creal(z[u]), cimag(z[u]));
        visited[u] = 1;
        u = visited[graph[u][0]] ? graph[u][1] : graph[u][0];
    }

    glp_delete_prob(ip);
    free(z);
    free(visited);
    for (size_t i = 0; i < n; i++)
        free(graph[i]);
    free(graph);
    free(subtour);
}