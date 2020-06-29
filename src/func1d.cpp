#include "func1d.hpp"

using namespace std;

static double simple_pow(double x, double pow)
{
    if (pow == 0) {
        return 1.0;
    }
    if (pow < 0) {
        pow = -pow;
        x = 1.0 / x;
    }
    double result = 1.0;
    for (int i = 0; i < pow; ++i) {
        result *= x;
    }
    return result;
}

polynomial_t::polynomial_t()
    : _params(1)
{
    _params[0] = 0.0;
}

polynomial_t::polynomial_t(double param0)
    : _params(1)
{
    _params[0] = param0;
}

polynomial_t::polynomial_t(double param0, double param1)
    : _params(2)
{
    _params[0] = param0;
    _params[1] = param1;
}

polynomial_t::polynomial_t(const polynomial_t& copy)
    : _params(copy._params)
{
}

double polynomial_t::func(double x) const
{
    double result = 0.0;
    if (!_params.empty()) {
        result += _params[0];
        for (int i = 1; i < _params.size(); ++i) {
            result += _params[i] * simple_pow(x, i);
        }
    }
    return result;
}

product_t::product_t(const shared_ptr<func1d_t>& func0)
        : _funcs(1)
{
    _funcs[0] = func0;
}

product_t::product_t(const shared_ptr<func1d_t>& func0, const shared_ptr<func1d_t>& func1)
        : _funcs(2)
{
    _funcs[0] = func0;
    _funcs[1] = func1;
}

product_t::product_t(
        const shared_ptr<func1d_t>& func0, const shared_ptr<func1d_t>& func1, const shared_ptr<func1d_t>& func2)
        : _funcs(3)
{
    _funcs[0] = func0;
    _funcs[1] = func1;
    _funcs[2] = func2;
}

product_t::product_t(const shared_ptr<func1d_t>& func0, const shared_ptr<func1d_t>& func1,
                     const shared_ptr<func1d_t>& func2, const shared_ptr<func1d_t>& func3)
        : _funcs(4)
{
    _funcs[0] = func0;
    _funcs[1] = func1;
    _funcs[2] = func2;
    _funcs[3] = func3;
}

product_t::product_t(const product_t& copy)
        : _funcs(copy._funcs)
{
}

double product_t::func(double x) const
{
    double result = 1.0;
    for (const shared_ptr<func1d_t>& func : _funcs) {
        result *= func->func(x);
    }
    return result;
}

sum_func_t::sum_func_t(const shared_ptr<func1d_t>& func0)
        : _funcs(1)
{
    _funcs[0] = func0;
}

sum_func_t::sum_func_t(const shared_ptr<func1d_t>& func0, const shared_ptr<func1d_t>& func1)
        : _funcs(2)
{
    _funcs[0] = func0;
    _funcs[1] = func1;
}

sum_func_t::sum_func_t(
        const shared_ptr<func1d_t>& func0, const shared_ptr<func1d_t>& func1, const shared_ptr<func1d_t>& func2)
        : _funcs(3)
{
    _funcs[0] = func0;
    _funcs[1] = func1;
    _funcs[2] = func2;
}

sum_func_t::sum_func_t(const shared_ptr<func1d_t>& func0, const shared_ptr<func1d_t>& func1,
                     const shared_ptr<func1d_t>& func2, const shared_ptr<func1d_t>& func3)
        : _funcs(4)
{
    _funcs[0] = func0;
    _funcs[1] = func1;
    _funcs[2] = func2;
    _funcs[3] = func3;
}

sum_func_t::sum_func_t(const sum_func_t& copy)
        : _funcs(copy._funcs)
{
}

double sum_func_t::func(double x) const
{
    double result = 0.0;
    for (const shared_ptr<func1d_t>& func : _funcs) {
        result += func->func(x);
    }
    return result;
}

composition_t::composition_t(const shared_ptr<func1d_t>& outer, const shared_ptr<func1d_t>& inner)
    : _outer(outer)
    , _inner(inner)
{
}

composition_t::composition_t(const composition_t& copy)
    : _outer(copy._outer)
    , _inner(copy._inner)
{
}

double composition_t::func(double x) const { return _outer->func(_inner->func(x)); }
