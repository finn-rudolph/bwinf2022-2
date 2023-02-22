#include <bits/stdc++.h>
#include "Highs.h"
using namespace std;

size_t nchoose2(size_t n) { return n * (n - 1) / 2; }

size_t edge_index(size_t n, size_t i, size_t j)
{
    return nchoose2(n) - nchoose2(n - min(i, j)) + max(i, j) - min(i, j) - 1;
}

bool is_acute(complex<double> a, complex<double> b, complex<double> c)
{
    return (b.real() - a.real()) * (c.real() - b.real()) +
               (b.imag() - a.imag()) * (c.imag() - b.imag()) <
           0.0;
}

void ban_triangle(HighsModel &model, size_t n, size_t i, size_t j, size_t k)
{
    model.lp_.a_matrix_.start_.push_back(model.lp_.a_matrix_.index_.size());
    model.lp_.row_lower_.push_back(0);
    model.lp_.row_upper_.push_back(1);
    model.lp_.a_matrix_.index_.push_back(edge_index(n, i, j));
    model.lp_.a_matrix_.index_.push_back(edge_index(n, j, k));
    model.lp_.a_matrix_.index_.push_back(edge_index(n, k, i));
    model.lp_.a_matrix_.value_.push_back(1);
    model.lp_.a_matrix_.value_.push_back(1);
    model.lp_.a_matrix_.value_.push_back(1);
}

void ban_triple(HighsModel &model, size_t n, size_t i, size_t j, size_t k)
{
    model.lp_.a_matrix_.start_.push_back(model.lp_.a_matrix_.index_.size());
    model.lp_.row_lower_.push_back(0);
    model.lp_.row_upper_.push_back(1);
    model.lp_.a_matrix_.index_.push_back(edge_index(n, i, j));
    model.lp_.a_matrix_.index_.push_back(edge_index(n, j, k));
    model.lp_.a_matrix_.value_.push_back(1);
    model.lp_.a_matrix_.value_.push_back(1);
}

void add_angle_constraints(HighsModel &model, vector<complex<double>> &z)
{
    size_t x = model.lp_.row_lower_.size();
    for (size_t i = 0; i < z.size(); i++)
    {
        for (size_t j = i + 1; j < z.size(); j++)
        {
            for (size_t k = j + 1; k < z.size(); k++)
            {
                if (is_acute(z[i], z[j], z[k]) && is_acute(z[k], z[i], z[j]) &&
                    is_acute(z[j], z[k], z[i]))
                {
                    ban_triangle(model, z.size(), i, j, k);
                }
                else
                {
                    if (is_acute(z[i], z[j], z[k]))
                        ban_triple(model, z.size(), i, j, k);
                    if (is_acute(z[k], z[i], z[j]))
                        ban_triple(model, z.size(), k, i, j);
                    if (is_acute(z[j], z[k], z[i]))
                        ban_triple(model, z.size(), j, k, i);
                }
            }
        }
    }
}

void add_degree_constraints(HighsModel &model, size_t n)
{
    for (size_t i = 0; i < n; i++)
    {
        model.lp_.a_matrix_.start_.push_back(model.lp_.a_matrix_.index_.size());
        model.lp_.row_lower_.push_back(1);
        model.lp_.row_upper_.push_back(2);
        for (size_t j = 0; j < n; j++)
        {
            if (i != j)
            {
                model.lp_.a_matrix_.index_.push_back(edge_index(n, i, j));
                model.lp_.a_matrix_.value_.push_back(1);
            }
        }
    }

    model.lp_.a_matrix_.start_.push_back(model.lp_.a_matrix_.index_.size());
    model.lp_.row_lower_.push_back(n - 1);
    model.lp_.row_upper_.push_back(n - 1);
    for (size_t i = 0; i < nchoose2(n); i++)
    {
        model.lp_.a_matrix_.index_.push_back(i);
        model.lp_.a_matrix_.value_.push_back(1);
    }
}

void add_subtour_elimination_constraint(
    Highs &highs, size_t n, vector<size_t> const &tour)
{
    HighsInt *ind = (HighsInt *)malloc(nchoose2(tour.size()) * sizeof *ind);
    double *val = (double *)malloc(nchoose2(tour.size()) * sizeof *val);
    size_t k = 0;

    for (size_t i = 0; i < tour.size(); i++)
    {
        for (size_t j = i + 1; j < tour.size(); j++)
        {
            ind[k] = edge_index(n, tour[i], tour[j]);
            val[k] = 1;
            k++;
        }
    }

    highs.addRow(0, tour.size() - 1, nchoose2(tour.size()), ind, val);
    free(ind);
    free(val);
}

bool check_for_subtours(Highs &highs, size_t n)
{
    vector<vector<unsigned>> graph(n);
    HighsSolution const &solution = highs.getSolution();

    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
            if (solution.col_value[edge_index(n, i, j)] > 0.5)
                graph[i].push_back(j), graph[j].push_back(i);

    vector<bool> visited(n, 0);
    bool has_subtours = 0;

    for (size_t i = 0; i < n; i++)
    {
        if (!visited[i])
        {
            vector<size_t> subtour;
            size_t u = i, last = SIZE_MAX;

            do
            {
                visited[u] = 1;
                subtour.push_back(u);
                size_t next = SIZE_MAX;
                for (size_t v : graph[u])
                    if (v != last)
                        next = v;
                last = u;
                u = next;
            } while (u != i && u != SIZE_MAX);

            if (u == i)
            {
                has_subtours = 1;
                add_subtour_elimination_constraint(highs, n, subtour);
            }
        }
    }

    return has_subtours;
}

int main()
{
    vector<complex<double>> z;
    double x, y;
    while (scanf("%lf %lf", &x, &y) == 2)
        z.emplace_back(x, y);
    size_t const n = z.size();

    HighsModel model;
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    model.lp_.num_col_ = nchoose2(n);

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i + 1; j < n; j++)
        {
            model.lp_.col_cost_.push_back(abs(z[i] - z[j]));
            model.lp_.integrality_.push_back(HighsVarType::kInteger);
            model.lp_.col_lower_.push_back(0);
            model.lp_.col_upper_.push_back(1);
        }
    }

    add_angle_constraints(model, z);
    add_degree_constraints(model, n);
    model.lp_.a_matrix_.start_.push_back(model.lp_.a_matrix_.index_.size());
    model.lp_.num_row_ = model.lp_.row_lower_.size();

    Highs highs;
    HighsStatus status;
    status = highs.passModel(model);
    assert(status == HighsStatus::kOk);
    bool has_subtours = 1;

    while (has_subtours && highs.getModelStatus() != HighsModelStatus::kInfeasible)
    {
        status = highs.run();
        assert(status == HighsStatus::kOk);
        has_subtours = check_for_subtours(highs, n);
    }

    if (highs.getModelStatus() == HighsModelStatus::kInfeasible)
    {
        cout << "Keine Tour mit Abbiegewinkeln von höchstens 90 Grad möglich.\n";
    }
    else
    {
        cout << "Tourlänge: " << highs.getObjectiveValue() << '\n';

        HighsSolution const &solution = highs.getSolution();
        vector<vector<size_t>> graph(n);

        for (size_t i = 0; i < n; i++)
            for (size_t j = i + 1; j < n; j++)
                if (solution.col_value[edge_index(n, i, j)] > 0.5)
                    graph[i].push_back(j), graph[j].push_back(i);

        size_t u, last = SIZE_MAX;
        for (size_t i = 0; i < n; i++)
            if (graph[i].size() == 1)
                u = i;

        do
        {
            cout << z[u].real() << ' ' << z[u].imag() << '\n';
            size_t next = SIZE_MAX;
            for (size_t v : graph[u])
                if (v != last)
                    next = v;
            last = u;
            u = next;
        } while (u != SIZE_MAX);
    }
}