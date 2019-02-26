#include "stdafx.h"
#include <iostream>
#include <string>
#include <vector>


//          ______     _                             _       _ 
//          | ___ \   | |                           (_)     | |
//          | |_/ ___ | |_   _ _ __   ___  _ __ ___  _  __ _| |
//          |  __/ _ \| | | | | '_ \ / _ \| '_ ` _ \| |/ _` | |
//          | | | (_) | | |_| | | | | (_) | | | | | | | (_| | |
//          \_|  \___/|_|\__, |_| |_|\___/|_| |_| |_|_|\__,_|_|
//                        __/ |                                
//                       |___/                                 

// Simple Polynomial Template class mainly for integer coefficient Polynomials.
// Coefficients are stored in a dense format as a vector of T.
//
// Functionality:
// - constructors: from single coeficient T, vector of coefficients or from first-last iterators.
// - internal data size: size() and Degree() (0 Polynomial has Degree of -1)
// - data access: [] operator (only onst version), begin(), end() iterators for coeficient vector
// - math operators: +; -; *, +=; -=; *=;
// - polynomial composition: & operator;
// - equality checks: operators == and != (also works for comparing Polynomials of Degree 0 and T)
// - writing to an ofstream: << operator (writes the Polynomial in a nice format)
//

using namespace std;

template<typename T>
T pow(T val, size_t power) {
    T res = T(1);
    while (power) {
        if (!(power % 2)) {
            power /= 2;
            val *= val;
        } else {
            --power;
            res *= val;
        }
    }
    return res;
}

template<typename T>
class Polynomial {
private:
    //
    //  DATA
    //
    vector<T> coefficients;
    //
    // Internal functions
    //
    void removeLeadingZeroes() {
        size_t i = 1;
        for (; i <= coefficients.size(); ++i) {
            if (coefficients[coefficients.size() - i] != T(0))
                break;
        }
        --i;
        if (coefficients.size() - i < coefficients.size() && coefficients.size() - i >= 0) {
            coefficients.resize(coefficients.size() - i);
            coefficients.shrink_to_fit();
        }
    }
    // Sign = 0 is "+="; Sign = 1 is "-="
    void Add(const Polynomial<T>& other, bool sign) {
        if (coefficients.size() < other.coefficients.size())
            coefficients.resize(other.coefficients.size());
        for (size_t i = 0; i < other.coefficients.size(); ++i)
            if (sign)
                coefficients[i] -= other.coefficients[i];
            else
                coefficients[i] += other.coefficients[i];
        removeLeadingZeroes();
    }

public:
    //
    // Constructors
    //
    Polynomial<T>() : coefficients() {}
    Polynomial<T>(const vector<T>& coef) : coefficients(coef) {
        removeLeadingZeroes();
    }
    Polynomial<T>(const T& coef) : coefficients(1) {
        coefficients[0] = coef;
        removeLeadingZeroes();
    }
    template<typename It>
    Polynomial<T>(It first, It last) {
        while (first != last) {
            coefficients.push_back(*first);
            ++first;
        }
        removeLeadingZeroes();
    }
    //
    // Size
    //
    Polynomial<T>(const Polynomial<T>& other) : coefficients(other.coefficients) {}
    int Degree() const {
        return coefficients.size() - 1;
    }
    size_t size() const {
        return coefficients.size();
    }
    //
    // Data access
    //

    T operator[] (size_t i) const {
        T res = T(0);
        if (i < coefficients.size())
            res = coefficients[i];
        return res;
    }
    vector<T>& GetCoef() {
        return coefficients;
    }
    vector<T> GetCoef() const {
        return coefficients;
    }
    //  Iterators
    auto begin() {
        return coefficients.begin();
    }
    auto end() {
        return coefficients.end();
    }
    auto begin() const {
        return coefficients.cbegin();
    }
    auto end() const {
        return coefficients.cend();
    }
    //
    // Math operations
    //
    //      Add / Substract
    Polynomial<T>& operator += (const Polynomial<T>& other) {
        Add(other, false);
        removeLeadingZeroes();
        return *this;
    }
    Polynomial<T>& operator -= (const Polynomial<T>& other) {
        Add(other, true);
        removeLeadingZeroes();
        return *this;
    }
    Polynomial<T> operator + (const Polynomial<T>& other) const {
        Polynomial<T> temp(*this);
        temp.Add(other, false);
        return temp;
    }
    Polynomial<T> operator - (const Polynomial<T>& other) const {
        Polynomial<T> temp(*this);
        temp.Add(other, true);
        return temp;
    }
    // Scalar version
    Polynomial<T>& operator += (const T& other) {
        if (coefficients.size())
            coefficients[0] += other;
        else
            coefficients.push_back(other);
        removeLeadingZeroes();
        return *this;
    }
    Polynomial<T>& operator -= (const T& other) {
        if (!coefficients.size())
            coefficients.push_back(T(0));
        coefficients[0] -= other;
        removeLeadingZeroes();
        return *this;
    }
    Polynomial<T> operator + (const T& other) const {
        Polynomial<T> temp(*this);
        return (temp += other);
    }
    friend Polynomial<T> operator + (const T& other, const Polynomial<T>& pol) {
        Polynomial<T> temp(other);
        temp += pol;
        return temp;
    }
    Polynomial<T> operator - (const T& other) const {
        Polynomial<T> temp(*this);
        return temp -= other;
    }
    friend Polynomial<T> operator - (const T& other, const Polynomial<T>& pol) {
        Polynomial<T> temp(other);
        temp -= pol;
        return temp;
    }
    //      Multiply
    Polynomial<T>& operator *= (const Polynomial<T>& other) {
        Polynomial<T> temp(*this);
        std::fill(coefficients.begin(), coefficients.end(), T(0));
        if (coefficients.size() + other.coefficients.size() >= 1)
            coefficients.resize(coefficients.size() + other.coefficients.size() - 1);
        for (size_t i = 0; i < temp.coefficients.size(); ++i) {
            for (size_t j = 0; j < other.coefficients.size(); ++j) {
                coefficients[i + j] += temp.coefficients[i] * other.coefficients[j];
            }
        }
        removeLeadingZeroes();
        return *this;
    }
    Polynomial<T> operator * (const Polynomial<T>& other) const {
        Polynomial<T> temp;
        if (coefficients.size() + other.coefficients.size() >= 1)
            temp.coefficients.resize(coefficients.size() + other.coefficients.size() - 1);
        for (size_t i = 0; i < coefficients.size(); ++i) {
            for (size_t j = 0; j < other.coefficients.size(); ++j) {
                temp.coefficients[i + j] += coefficients[i] * other.coefficients[j];
            }
        }
        temp.removeLeadingZeroes();
        return temp;
    }
    // Scalar version
    Polynomial<T>& operator *= (const T& other) {
        if (other == T(0)) {
            coefficients.resize(0);
            coefficients.shrink_to_fit();
        } else {
            for (size_t i = 0; i < coefficients.size(); ++i)
                coefficients[i] *= other;
        }
        removeLeadingZeroes();
        return *this;
    }
    Polynomial<T> operator * (const T& other) const {
        Polynomial<T> temp(*this);
        return (temp *= other);
    }
    friend Polynomial<T> operator * (const T& other, const Polynomial<T>& pol) {
        Polynomial<T> temp(other);
        return (temp *= pol);
    }
    // Calculate value at given point
    T operator() (const T& val) const {
        T res = T(0);
        for (size_t i = 0; i < coefficients.size(); ++i)
            res += pow(val, i) * coefficients[i];
        return res;
    }
    //
    // Composition
    //
    template<typename Z>
    friend Polynomial<Z> operator & (const Polynomial<Z>& lhs, const Polynomial<Z>& rhs);
    // Equality
    template<typename Z>
    friend bool operator == (const Polynomial<Z>& lhs, const Polynomial<Z>& rhs);
    template<typename Z>
    friend bool operator != (const Polynomial<Z>& lhs, const Polynomial<Z>& rhs);

    // Writing to ostream
};
template<typename T>
std::ostream& operator << (std::ostream& out, const Polynomial<T>& pol) {
    if (pol.size()) {
        for (size_t i = 1; i <= pol.size(); ++i) {
            if (pol[pol.size() - i] != T(0)) {
                if (pol[pol.size() - i] > T(0) && i != 1)
                    out << '+';
                if ((pol[pol.size() - i] != T(1)
                     || pol.size() - i == 0) && pol[pol.size() - i] != T(-1)) {
                    out << pol[pol.size() - i];
                    if (pol.size() - i != 0)
                        out << '*';
                } else {
                    if (pol[pol.size() - i] < T(0))
                        out << '-';
                    if (pol.size() - i == 0)
                        out << T(1);
                }
                if (pol.size() - i > 0) {
                    out << "x";
                    if (pol.size() - i > 1)
                        out << "^" << (pol.size() - i);
                }
            }
        }
    } else {
        out << T(0);
    }
    return out;
}
template<typename T>
Polynomial<T> operator & (const Polynomial<T>& lhs, const Polynomial<T>& rhs) {
    T initializer = T(0);
    if (lhs.coefficients.size() > 0)
        initializer = lhs.coefficients[0];
    Polynomial<T> res(initializer);
    Polynomial<T> temp(T(1));
    res.coefficients.resize((lhs.coefficients.size() - 1) * (rhs.coefficients.size() - 1) + 1);
    for (size_t i = 1; i != lhs.coefficients.size(); ++i) {
        temp *= rhs.coefficients;  // raise to the i-th power
        res += lhs.coefficients[i] * temp;
    }
    res.removeLeadingZeroes();
    return res;
}
template<typename T>
bool operator == (const Polynomial<T>& lhs, const Polynomial<T>& rhs) {
    return (lhs.coefficients == rhs.coefficients);
}
template<typename T>
bool operator != (const Polynomial<T>& lhs, const Polynomial<T>& rhs) {
    return (lhs.coefficients != rhs.coefficients);
}
template<typename T>
bool operator == (const Polynomial<T>& lhs, const T& rhs) {
    return ((lhs.size() == 1 && lhs[0] == rhs) || (lhs.size() == 0 && rhs == T(0)));
}
template<typename T>
bool operator != (const Polynomial<T>& lhs, const T& rhs) {
    return (lhs == rhs);
}


// (C) TechnoVirus, February 2019.
