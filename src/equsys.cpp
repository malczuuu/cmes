#include "equsys.hpp"

using namespace std;

static bool _equals(double a, double b) { return abs(a - b) < 1e-9; }

void equsys_t::_diag_fix(int i)
{
    int index = i;
    for (int j = i; j < size(); ++j) {
        if (abs((double)get(j, i)) > abs((double)get(index, i))) {
            index = j;
        }
    }

    if (index != i) {
        swap(_inds[i], _inds[index]);
    }
}

equsys_t::equsys_t(int size)
    : _size(size)
{
    if (size > 0) {
        _a_matr.resize(size * size);
        _b_vect.resize(size);
        _inds.resize(size);
        for (int i = 0; i < size; ++i) {
            _inds[i] = i;
        }
    }
}

equsys_t::~equsys_t() = default;

const int& equsys_t::size() const { return _size; }

const double& equsys_t::get(int i, int j) const { return _a_matr[_coords(i, j)]; }

const double& equsys_t::get(int i) const { return _b_vect[_coords(i)]; }

void equsys_t::set(int i, int j, double value) { _a_matr[_coords(i, j)] = value; }

void equsys_t::set(int i, double value) { _b_vect[_coords(i)] = value; }

void equsys_t::inc(int i, int j, double value) { set(i, j, get(i, j) + value); }

void equsys_t::inc(int i, double value) { set(i, get(i) + value); }

void equsys_t::_decompose_lu()
{
    for (int i = 0; i < size(); ++i) {
        if (_equals(get(i, i), 0.0)) {
            _diag_fix(i);
            if (_equals(get(i, i), 0.0)) {
                continue;
            }
        }

        for (int j = i + 1; j < size(); ++j) {
            double param = get(j, i) / get(i, i);
            for (int k = i; k < size(); ++k) {
                set(j, k, get(j, k) - get(i, k) * param);
            }
            set(j, i, param);

            set(j, get(j) - get(i) * param);
        }
    }
}

vector<double> equsys_t::solve()
{
    _decompose_lu();
    vector<double> result(size());

    double sum = 0.0;
    vector<double> y_vect(size());
    y_vect[0] = get(0);
    for (int i = 1; i < size(); ++i) {
        for (int j = 0; j < i; ++j) {
            sum += get(i, j) * result[j];
        }
        y_vect[i] = get(i) - sum;
        sum = 0.0;
    }

    result[size() - 1] = y_vect[size() - 1] / get(size() - 1, size() - 1);

    for (int i = size() - 2; i >= 0; --i) {
        for (int j = i + 1; j < size(); ++j) {
            sum += get(i, j) * result[j];
        }
        result[i] = (y_vect[i] - sum) / get(i, i);
        sum = 0.0;
    }

    return result;
}
