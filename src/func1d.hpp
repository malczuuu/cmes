#ifndef FUNC1D_HPP
#define FUNC1D_HPP

#include <memory>
#include <string>
#include <vector>

/**
 * Abstract class for representing mathematical functions as objects.
 */
class func1d_t {
public:
    virtual double func(double x) const = 0;
};

/**
 * Class which represents objects of polynomial functions.
 */
class polynomial_t : public func1d_t {
private:
    std::vector<double> _params;

public:
    polynomial_t();

    /**
     * Constructs a constant-value function.
     */
    explicit polynomial_t(double param0);

    /**
     * Constructs a linear function. Note that order of parameters is inverted, so the function will match the following
     * param1 * x + param0.
     */
    polynomial_t(double param0, double param1);

    double func(double x) const override;
};

/**
 * Class which represents the multiplication product of lots of other functions.
 */
class product_t : public func1d_t {
private:
    std::vector<std::shared_ptr<func1d_t>> _funcs;

public:
    explicit product_t(const std::shared_ptr<func1d_t>& func0);

    product_t(const std::shared_ptr<func1d_t>& func0, const std::shared_ptr<func1d_t>& func1);

    product_t(const std::shared_ptr<func1d_t>& func0, const std::shared_ptr<func1d_t>& func1,
        const std::shared_ptr<func1d_t>& func2);

    product_t(const std::shared_ptr<func1d_t>& func0, const std::shared_ptr<func1d_t>& func1,
        const std::shared_ptr<func1d_t>& func2, const std::shared_ptr<func1d_t>& func3);

    double func(double x) const override;
};

/**
 * Class which represents the composition of two functions outer(inner(x)).
 */
class composition_t : public func1d_t {
private:
    std::shared_ptr<func1d_t> _outer;
    std::shared_ptr<func1d_t> _inner;

public:
    composition_t(const std::shared_ptr<func1d_t>& outer, const std::shared_ptr<func1d_t>& inner);

    double func(double x) const override;
};

#endif
