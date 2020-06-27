#include "main.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    shared_ptr<func1d_t> a_func = make_shared<polynomial_t>(1.0);
    shared_ptr<func1d_t> f_func = make_shared<polynomial_t>(-2.0);

    double alpha0 = 0.0;
    double beta0 = 1.0;
    double gamma0 = 1.0;

    double alpha1 = 1.0;
    double beta1 = 0.0;
    double gamma1 = 2.0;

    int elements = 3;
    double x0 = 0.0;
    double x1 = 1.0;

    mes_solver_t solver;

    solver.a_func(a_func);
    solver.f_func(f_func);

    solver.alpha0(alpha0);
    solver.beta0(beta0);
    solver.gamma0(gamma0);

    solver.alpha1(alpha1);
    solver.beta1(beta1);
    solver.gamma1(gamma1);

    solver.x0(x0);
    solver.x1(x1);
    solver.elements(elements);

    solver.solve();

    return 0;
}
