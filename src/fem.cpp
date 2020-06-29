#include "fem.hpp"
#include "equsys.hpp"

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

fem_t::fem_t()
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

vector<array<double, 2>> fem_t::solve() const
{
    int nodes = this->nodes();
    double el_size = this->element_size();

    vector<array<double, 2>> points(nodes);
    vector<shared_ptr<func1d_t>> funcs(nodes);

    equsys_t equsys(nodes);

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

        equsys.inc(i, i, integral(*stiffness11func, 0.0, 1.0));
        equsys.inc(i, i + 1, integral(*stiffness12func, 0.0, 1.0));
        equsys.inc(i + 1, i, integral(*stiffness12func, 0.0, 1.0));
        equsys.inc(i + 1, i + 1, integral(*stiffness22func, 0.0, 1.0));

        shared_ptr<func1d_t> load1fun = make_shared<product_t>(_f_func, ref_shape1, ref_trans_dx);
        shared_ptr<func1d_t> load2fun = make_shared<product_t>(_f_func, ref_shape2, ref_trans_dx);
        shared_ptr<func1d_t> ref_trans_inv = make_shared<polynomial_t>(-xa / (xb - xa), 1.0 / (xb - xa));
        shared_ptr<func1d_t> shape1 = make_shared<composition_t>(ref_shape1, ref_trans_inv);
        shared_ptr<func1d_t> shape2 = make_shared<composition_t>(ref_shape2, ref_trans_inv);

        equsys.inc(i, integral(*load1fun, 0.0, 1.0));
        equsys.inc(i + 1, integral(*load2fun, 0.0, 1.0));

        funcs[i] = shape1;
        funcs[i + 1] = shape1;

        points[i][0] = xa;
        points[i + 1][0] = xb;
        xa += el_size;
        xb += el_size;

        if (i == 0 && _alpha0 != 0.0) {
            equsys.inc(i, i, _a_func->func(_x0) * _beta0 / _alpha0 * shape1->func(_x0) * shape1->func(_x0));
            equsys.inc(i, _a_func->func(_x0) * _gamma0 / _alpha0 * shape1->func(_x0));
        } else if (i == _elements - 1 && _alpha1 != 0.0) {
            equsys.inc(i + 1, i + 1, _a_func->func(_x1) * _beta1 / _alpha1 * shape2->func(_x1) * shape2->func(_x1));
            equsys.inc(i + 1, _a_func->func(_x1) * _gamma1 / _alpha1 * shape2->func(_x1));
        }
    }

    if (_alpha0 == 0.0) {
        equsys.set(0, 0, 1);
        for (int i = 1; i < nodes; ++i) {
            equsys.set(0, i, 0);
        }
        equsys.set(0, _gamma0 / _beta0);
    }
    if (_alpha1 == 0.0) {
        int last = nodes = 1;
        equsys.set(last, last, 1);
        for (int i = 0; i < last; ++i) {
            equsys.set(last, i, 0);
        }
        equsys.set(last, _gamma1 / _beta1);
    }

    vector<double> results = equsys.solve();

    for (int i = 0; i < nodes; ++i) {
        points[i][1] = results[i];
    }

    return points;
}

void fem_t::a_func(shared_ptr<func1d_t> value) { _a_func = value; }

void fem_t::b_func(shared_ptr<func1d_t> value) { _b_func = value; }

void fem_t::c_func(shared_ptr<func1d_t> value) { _c_func = value; }

void fem_t::f_func(shared_ptr<func1d_t> value) { _f_func = value; }

void fem_t::alpha0(double value) { _alpha0 = value; }

void fem_t::beta0(double value) { _beta0 = value; }

void fem_t::gamma0(double value) { _gamma0 = value; }

void fem_t::alpha1(double value) { _alpha1 = value; }

void fem_t::beta1(double value) { _beta1 = value; }

void fem_t::gamma1(double value) { _gamma1 = value; }

void fem_t::x0(double value) { _x0 = value; }

void fem_t::x1(double value) { _x1 = value; }

void fem_t::elements(int value) { _elements = value; }

int fem_t::nodes() const { return _elements + 1; }

double fem_t::element_size() const { return (_x1 - _x0) / _elements; }
