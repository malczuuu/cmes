#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "func1d.hpp"

class point2d_t {
private:
    double _x;
    double _y;

public:
    point2d_t(double x = 0.0, double y = 0.0);

    point2d_t(const point2d_t& copy);

    double x() const;

    double y() const;

    std::string str();
};

class mes_solver_t {
private:
    std::shared_ptr<func1d_t> _a_func;
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
    mes_solver_t();

    std::vector<point2d_t> solve() const;

    void a_func(std::shared_ptr<func1d_t> value);

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
