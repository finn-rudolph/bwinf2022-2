#include <bits/stdc++.h>
#include "util.hpp"
#include "Highs.h"
using namespace std;

// Gibt den Index der zur Kante {i, j} zugehörigen Variable zurück.
size_t edge_index(size_t n, size_t i, size_t j)
{
    return nchoose2(n) - nchoose2(n - min(i, j)) + max(i, j) - min(i, j) - 1;
}

// Fügt für jedes Tripel i, j, k (i != j != k, i < k) die Bedingung hinzu,
// dass die Kanten ij und jk nicht gleichzeitig verwendet werden dürfen, wenn
// der Betrag ihres Außenwinkels > pi / 2 ist.
void add_angle_constraints(HighsModel &model, vector<complex<double>> const &z)
{
    for (size_t j = 0; j < z.size(); j++)
    {
        for (size_t i = 0; i < z.size(); i++)
        {
            if (i == j)
                continue;
            for (size_t k = i + 1; k < z.size(); k++)
            {
                if (k != j && dot_product(z[j] - z[i], z[k] - z[j]) < 0.0)
                {
                    HighsLp &lp = model.lp_;
                    lp.a_matrix_.index_.push_back(edge_index(z.size(), i, j));
                    lp.a_matrix_.index_.push_back(edge_index(z.size(), j, k));
                    lp.a_matrix_.value_.push_back(1); // Der Koeffizient jeder
                    lp.a_matrix_.value_.push_back(1); // Kante ist 1.
                    lp.row_lower_.push_back(0);
                    lp.row_upper_.push_back(1);
                    //
                    lp.a_matrix_.start_.push_back(lp.a_matrix_.index_.size());
                }
            }
        }
    }
}

// Schränkt den Grad jedes Knoten auf 1 oder 2 ein.
void add_degree_constraints(HighsModel &model, size_t n)
{
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            if (i != j) // Die Kante zu jedem Knoten != i wird mit Koeffizient
            {           // 1 zur aktuellen Zeile hinzugefügt.
                model.lp_.a_matrix_.index_.push_back(edge_index(n, i, j));
                model.lp_.a_matrix_.value_.push_back(1);
            }
        }
        model.lp_.row_lower_.push_back(1); // Setze Unter- und Oberschranke der
        model.lp_.row_upper_.push_back(2); // Zeile.
        model.lp_.a_matrix_.start_.push_back(model.lp_.a_matrix_.index_.size());
    }
}

// Schränkt die Anzahl verwendeter Kanten auf genau n - 1 ein.
void add_num_edges_constraint(HighsModel &model, size_t n)
{
    for (size_t i = 0; i < nchoose2(n); i++) // Iteriere über alle Kanten.
    {
        model.lp_.a_matrix_.index_.push_back(i);
        model.lp_.a_matrix_.value_.push_back(1);
    }
    model.lp_.row_lower_.push_back(n - 1); // Fixiere den Wert der neuen Zeile
    model.lp_.row_upper_.push_back(n - 1); // auf genau n - 1.
    model.lp_.a_matrix_.start_.push_back(model.lp_.a_matrix_.index_.size());
}

// Fügt einen SEC für die Knoten in tour ein.
void add_subtour_elimination_constraint(
    Highs &highs, size_t n, vector<size_t> const &tour)
{
    // Arrays für die neuen Spaltenindizes und Werte.
    HighsInt *ind = (HighsInt *)malloc(nchoose2(tour.size()) * sizeof *ind);
    double *val = (double *)malloc(nchoose2(tour.size()) * sizeof *val);
    size_t k = 0; // Anzahl bereits in ind bzw. val eingefügter Elemente

    // Setze den Koeffizienten jeder Kante zwischen Knoten der Subtour auf 1.
    for (size_t i = 0; i < tour.size(); i++)
    {
        for (size_t j = i + 1; j < tour.size(); j++)
        {
            ind[k] = edge_index(n, tour[i], tour[j]);
            val[k] = 1;
            k++;
        }
    }

    // Verwendet Highs::addRow(double lower, double upper, HighsInt num_new_nz,
    //                         const HighsInt *indices, const double *values)
    HighsStatus status;
    status = highs.addRow(0, tour.size() - 1, nchoose2(tour.size()), ind, val);
    assert(status == HighsStatus::kOk);
    free(ind);
    free(val);
}

// Gibt den als Lösung gefundenen Graphen in Form einer Adjazenzliste zurück.
vector<vector<size_t>> build_graph(Highs const &highs, size_t n)
{
    HighsSolution const &solution = highs.getSolution();
    vector<vector<size_t>> graph(n); // Adjazenzliste

    for (size_t i = 0; i < n; i++)         // Überprüfe für jede Kante, ob der
        for (size_t j = i + 1; j < n; j++) // Wert ihrer Variablen 1 ist.
            if (solution.col_value[edge_index(n, i, j)] > 0.5)
                graph[i].push_back(j), graph[j].push_back(i);

    return graph;
}

// Gibt zurück, ob in der Lösung Subtouren existieren.
bool check_for_subtours(Highs &highs, size_t n)
{
    vector<vector<size_t>> graph = build_graph(highs, n);
    vector<bool> visited(n, 0);
    bool has_subtours = 0;

    for (size_t i = 0; i < n; i++)
    {
        if (!visited[i])
        {
            vector<size_t> subtour;
            size_t j = i, last = SIZE_MAX; // aktueller und vorheriger Knoten

            do
            {
                visited[j] = 1;
                subtour.push_back(j);
                size_t next = SIZE_MAX;
                for (size_t k : graph[j]) // Da der Grad jedes Knoten <= 2 ist,
                    if (k != last)        // wurde jeder Knoten unterschiedlich
                        next = k;         // zum letzten noch nicht besucht oder
                last = j;                 // ein Zyklus gefunden.
                j = next;
            } while (j != i && j != SIZE_MAX);

            if (j == i) // Zyklus gefunden.
            {
                has_subtours = 1;
                add_subtour_elimination_constraint(highs, n, subtour);
            }
        }
    }

    return has_subtours;
}

// Gibt den kürzesten Hamiltonpfad zurück, auf dem der Abbiegewinkel jedes
// benachbarten Kantenpaares <= pi / 2 ist. Existiert kein solcher Pfad, ist der
// zurückgegebene Vektor leer.
vector<complex<double>> shortest_obtuse_path(vector<complex<double>> const &z)
{
    size_t const n = z.size();

    HighsModel model; // Objekt, in dem das IP spezifiziert wird.
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    model.lp_.a_matrix_.start_ = {0}; // Die erste Zeile beginnt bei Index 0.
    model.lp_.num_col_ = nchoose2(n);

    // Füge die Länge jeder Kante als ihren Koeffizienten zur Kostenfunktion
    // hinzu und beschränke ihre Variable auf 0 oder 1.
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
        {
            model.lp_.col_cost_.push_back(abs(z[i] - z[j]));
            model.lp_.integrality_.push_back(HighsVarType::kInteger);
            model.lp_.col_lower_.push_back(0);
            model.lp_.col_upper_.push_back(1);
        }

    add_angle_constraints(model, z);
    add_degree_constraints(model, n);
    add_num_edges_constraint(model, n);
    model.lp_.num_row_ = model.lp_.row_lower_.size();

    Highs highs;
    HighsStatus status;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    status = highs.passModel(model); // Übergebe model an highs.
    assert(status == HighsStatus::kOk);
    bool has_subtours = 1;

    while (has_subtours && model_status != HighsModelStatus::kInfeasible)
    {
        status = highs.run(); // Löse das ganzzahlige lineare Programm.
        assert(status == HighsStatus::kOk);
        model_status = highs.getModelStatus();
        has_subtours = check_for_subtours(highs, n); // Füge SECs ein.
    }

    if (model_status == HighsModelStatus::kInfeasible) // Keine Tour möglich.
        return vector<complex<double>>();

    vector<vector<size_t>> graph = build_graph(highs, n);
    vector<complex<double>> path;

    size_t j = SIZE_MAX, last = SIZE_MAX; // aktueller und vorheriger Knoten
    for (size_t i = 0; i < n; i++)
        if (graph[i].size() == 1) // Ein Knoten mit Grad 1 muss Anfang des Pfads
            j = i;                // sein.

    while (j != SIZE_MAX)
    {
        path.push_back(z[j]);
        size_t next = SIZE_MAX;
        for (size_t k : graph[j]) // Wähle den der maximal zwei benachbarten
            if (k != last)        // Knoten als Nachfolger, der nicht Vorgänger
                next = k;         // ist.
        last = j;
        j = next;
    }

    return path;
}

int main()
{
    vector<complex<double>> z;
    double x, y;
    while (scanf("%lf %lf", &x, &y) == 2)
        z.emplace_back(x, y);

    vector<complex<double>> const path = shortest_obtuse_path(z);

    if (path.empty())
    {
        cout << "Keine Tour mit maximalem Abbiegewinkel von 90 Grad möglich.\n";
    }
    else
    {
        cout << setprecision(6) << fixed
             << "Tourlänge: " << path_length(path) << '\n';
        for (complex<double> const &u : path)
            cout << u.real() << ' ' << u.imag() << '\n';
    }
}