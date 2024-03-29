#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "func1d.hpp"
#include <array>

class fem_t {
private:
    std::shared_ptr<func1d_t> _a_func;
    std::shared_ptr<func1d_t> _b_func;
    std::shared_ptr<func1d_t> _c_func;
    std::shared_ptr<func1d_t> _f_func;
    double _alpha0;
    double _beta0;
    double _gamma0;
    double _alpha1;
    double _beta1;
    double _gamma1;
    double _x0;
    double _x1;
    int _elements;

public:
    fem_t();

    std::vector<std::array<double, 2>> solve() const;

    void a_func(std::shared_ptr<func1d_t> value);

    void b_func(std::shared_ptr<func1d_t> value);

    void c_func(std::shared_ptr<func1d_t> value);

    void f_func(std::shared_ptr<func1d_t> value);

    void alpha0(double value);

    void beta0(double value);

    void gamma0(double value);

    void alpha1(double value);

    void beta1(double value);

    void gamma1(double value);

    void x0(double value);

    void x1(double value);

    void elements(int value);

    int nodes() const;

    double element_size() const;
};

#endif
