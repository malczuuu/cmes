#include "main.hpp"
#include <iomanip>
#include <iostream>

const double A = 1.0;
const double F = -2.0;

const double ALPHA0 = 0.0;
const double BETA0 = 1.0;
const double GAMMA0 = 1.0;

const double ALPHA1 = 1.0;
const double BETA1 = 0.0;
const double GAMMA1 = 2.0;

const int X0 = 0.0;
const int X1 = 1.0;

const int ELEMENTS = 4;

using namespace std;

int main(int argc, char* argv[])
{
    shared_ptr<func1d_t> a_func = make_shared<polynomial_t>(A);
    shared_ptr<func1d_t> b_func = make_shared<polynomial_t>();
    shared_ptr<func1d_t> c_func = make_shared<polynomial_t>();
    shared_ptr<func1d_t> f_func = make_shared<polynomial_t>(F);

    fem_t solver;

    solver.a_func(a_func);
    solver.b_func(b_func);
    solver.c_func(c_func);
    solver.f_func(f_func);

    solver.alpha0(ALPHA0);
    solver.beta0(BETA0);
    solver.gamma0(GAMMA0);

    solver.alpha1(ALPHA1);
    solver.beta1(BETA1);
    solver.gamma1(GAMMA1);

    solver.x0(X0);
    solver.x1(X1);
    solver.elements(ELEMENTS);

    vector<array<double, 2>> points = solver.solve();

    for (int i = 0; i < points.size(); ++i) {
        cout << setprecision(8) << fixed << points[i][0] << ";" << points[i][1] << endl;
    }

    return 0;
}
