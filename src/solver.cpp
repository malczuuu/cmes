#include "solver.hpp"
#include <iostream>

using namespace std;

static double integral(const func1d_t& func, double x0, double x1)
{
    double result = 0.0;
    double diff = (x1 - x0) / 1000.0;
    double xa = x0;
    double xb = xa + diff;
    for (int i = 0; i < 1000; ++i) {
        result += (func.func(xa) + func.func(xb)) * 0.5 * diff;
        xa += diff;
        xb += diff;
    }
    return result;
}

point2d_t::point2d_t(double x, double y)
    : _x(x)
    , _y(y)
{
}

point2d_t::point2d_t(const point2d_t& copy)
        : _x(copy._x)
        , _y(copy._y)
{
}

double point2d_t::x() const { return _x; }

double point2d_t::y() const { return _y; }

string point2d_t::str() { return "(" + to_string(_x) + ", " + to_string(_y) + ")"; }

mes_solver_t::mes_solver_t()
    : _a_func(make_shared<polynomial_t>(1.0))
    , _f_func(make_shared<polynomial_t>(-2.0))
    , _alpha0(0.0)
    , _beta0(1.0)
    , _gamma0(1.0)
    , _alpha1(1.0)
    , _beta1(0.0)
    , _gamma1(2.0)
    , _x0(0.0)
    , _x1(1.0)
    , _elements(3)
{
}

vector<point2d_t> mes_solver_t::solve() const
{
    int nodes = this->nodes();
    double el_size = this->element_size();

    vector<point2d_t> points(nodes);
    vector<shared_ptr<func1d_t>> funcs(nodes);
    vector<vector<double>> stiffness(nodes);
    for (int i = 0; i < nodes; ++i) {
        stiffness[i].resize(nodes);
    }
    vector<double> load(nodes);

    double xa = _x0;
    double xb = xa + el_size;

    for (int i = 0; i < _elements; ++i) {
        double ref_shape1_param0 = 1.0;
        double ref_shape1_param1 = -1.0;
        shared_ptr<func1d_t> ref_shape1 = make_shared<polynomial_t>(ref_shape1_param0, ref_shape1_param1);
        shared_ptr<func1d_t> ref_shape1_dx = make_shared<polynomial_t>(ref_shape1_param1);

        double ref_shape2_param0 = 0.0;
        double ref_shape2_param1 = 1.0;
        shared_ptr<func1d_t> ref_shape2 = make_shared<polynomial_t>(ref_shape2_param0, ref_shape2_param1);
        shared_ptr<func1d_t> ref_shape2_dx = make_shared<polynomial_t>(ref_shape2_param1);

        double ref_trans_param0 = xa;
        double ref_trans_param1 = xb - xa;
        shared_ptr<func1d_t> ref_trans = make_shared<polynomial_t>(ref_trans_param0, ref_trans_param1);
        shared_ptr<func1d_t> ref_trans_dx = make_shared<polynomial_t>(ref_trans_param1);

        double ref_trans_inv_param0 = 1.0 / ref_trans_param1;
        shared_ptr<func1d_t> ref_trans_inv_dx = make_shared<polynomial_t>(ref_trans_inv_param0);
        shared_ptr<func1d_t> shape1_dx = make_shared<polynomial_t>(-1.0 * ref_trans_inv_param0);
        shared_ptr<func1d_t> shape2_dx = make_shared<polynomial_t>(1.0 * ref_trans_inv_param0);

        shared_ptr<func1d_t> stiffness11func = make_shared<product_t>(_a_func, shape1_dx, shape1_dx, ref_trans_dx);
        shared_ptr<func1d_t> stiffness12func = make_shared<product_t>(_a_func, shape1_dx, shape2_dx, ref_trans_dx);
        shared_ptr<func1d_t> stiffness22func = make_shared<product_t>(_a_func, shape2_dx, shape2_dx, ref_trans_dx);

        stiffness[i][i] += integral(*stiffness11func, 0.0, 1.0);
        stiffness[i][i + 1] += integral(*stiffness12func, 0.0, 1.0);
        stiffness[i + 1][i] += integral(*stiffness12func, 0.0, 1.0);
        stiffness[i + 1][i + 1] += integral(*stiffness22func, 0.0, 1.0);

        shared_ptr<func1d_t> load1fun = make_shared<product_t>(_f_func, ref_shape1, ref_trans_dx);
        shared_ptr<func1d_t> load2fun = make_shared<product_t>(_f_func, ref_shape2, ref_trans_dx);
        shared_ptr<func1d_t> ref_trans_inv = make_shared<polynomial_t>(-xa / (xb - xa), 1.0 / (xb - xa));
        shared_ptr<func1d_t> shape1 = make_shared<composition_t>(ref_shape1, ref_trans_inv);
        shared_ptr<func1d_t> shape2 = make_shared<composition_t>(ref_shape2, ref_trans_inv);

        load[i] += integral(*load1fun, 0.0, 1.0);
        load[i + 1] += integral(*load2fun, 0.0, 1.0);

        funcs[i] = shape1;
        funcs[i + 1] = shape1;

        xa += el_size;
        xb += el_size;

        if (i == 0 && _alpha0 != 0.0) {
            stiffness[i][i] += _a_func->func(_x0) * _beta0 / _alpha0 * shape1->func(_x0) * shape1->func(_x0);
            load[i] += _a_func->func(_x0) * _gamma0 / _alpha0 * shape1->func(_x0);
        } else if (i == _elements - 1 && _alpha1 != 0.0) {
            stiffness[i + 1][i + 1] += _a_func->func(_x1) * _beta1 / _alpha1 * shape2->func(_x1) * shape2->func(_x1);
            load[i + 1] += _a_func->func(_x1) * _gamma1 / _alpha1 * shape2->func(_x1);
        }
    }

    for (int i = 0; i < nodes; ++i) {

        for (int j = 0; j < nodes; ++j) {
            cout << stiffness[i][j] << " ";
        }
        cout << "; " << load[i] << endl;
    }

    return points;
}

void mes_solver_t::a_func(shared_ptr<func1d_t> value) { _a_func = value; }

void mes_solver_t::f_func(shared_ptr<func1d_t> value) { _f_func = value; }

void mes_solver_t::alpha0(double value) { _alpha0 = value; }

void mes_solver_t::beta0(double value) { _beta0 = value; }

void mes_solver_t::gamma0(double value) { _gamma0 = value; }

void mes_solver_t::alpha1(double value) { _alpha1 = value; }

void mes_solver_t::beta1(double value) { _beta1 = value; }

void mes_solver_t::gamma1(double value) { _gamma1 = value; }

void mes_solver_t::x0(double value) { _x0 = value; }

void mes_solver_t::x1(double value) { _x1 = value; }

void mes_solver_t::elements(int value) { _elements = value; }

int mes_solver_t::nodes() const { return _elements + 1; };

double mes_solver_t::element_size() const { return (_x1 - _x0) / _elements; };
