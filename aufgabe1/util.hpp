#include <bits/stdc++.h>
using namespace std;

template <typename T>
T nchoose2(T n) { return n * (n - 1) / 2; }

template <typename T>
T dot_product(complex<T> const &a, complex<T> const &b)
{
    return (a * conj(b)).real();
}