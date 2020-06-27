#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

using namespace std;

class equsys_t {
private:
    /**
     * Represents A in the system of linear equations matrix form (A * x = b).
     */
    vector<double> _a_matr;

    /**
     * Represents b in the system of linear equations matrix form (A * x = b).
     */
    vector<double> _b_vect;

    /**
     * Since the algorithms reorders rows in matrix we to avoid moving large arrays the program just swaps the indexes.
     */
    vector<int> _inds;

    int _size;

    int _coords(int i, int j) const { return _inds[i] * _size + j; }

    int _coords(int i) const { return _inds[i]; }

    /**
     * Reorders row with the next first row which has nonzero i-diagonal if such row exists. Otherwise does nothing.
     */
    void _diag_fix(int i);

    void _decompose_lu();

public:
    explicit equsys_t(int size);

    virtual ~equsys_t();

    const int& size() const;

    const double& get(int i, int j) const;

    const double& get(int i) const;

    void set(int i, int j, double value);

    void set(int i, double value);

    void inc(int i, int j, double value);

    void inc(int i, double value);

    /**
     *  { Ly = b
     *  { Ux = y
     *
     *  L = [1 0 0 0]    U = [* * * *]
     *      [* 1 0 0]        [0 * * *]
     *      [* * 1 0]        [0 0 * *]
     *      [* * * 1]        [0 0 0 *]
     *  INFO:
     *  When calculating Ly = b, there's no need to divide result by an element on main diagonal - there's always 1.
     *
     *  Note that this matrix changes the internal state of equation system object.
     */
    vector<double> solve();
};
