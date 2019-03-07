#pragma once
//

//          ______     _                             _       _ 
//          | ___ \   | |                           (_)     | |
//          | |_/ ___ | |_   _ _ __   ___  _ __ ___  _  __ _| |
//          |  __/ _ \| | | | | '_ \ / _ \| '_ ` _ \| |/ _` | |
//          | | | (_) | | |_| | | | | (_) | | | | | | | (_| | |
//          \_|  \___/|_|\__, |_| |_|\___/|_| |_| |_|_|\__,_|_|
//                        __/ |                                
//                       |___/                                 

// Simple Polynomial Template classes mainly for integer coefficient Polynomials.
// Coefficients are stored in a dense format in Polynomial<T> as a vector of T and in a sparse
// format in the PolynomialSparse<T> as a map<size_t, T>.
//
// Functionality:
// - constructors: from single coeficient T, vector of coefficients or from first-last iterators.
// - internal data size: size() and Degree() (0 Polynomial has Degree of -1)
// - data access: [] operator (only onst version), begin(), end() iterators for coeficient vector
// - math operators: +; -; *, +=; -=; *=, \=, \;
// - GCF of two polynomials: operator , (recommended to use (a, b))
// - polynomial composition: & operator;
// - equality checks: operators == and != (also works for comparing Polynomials of Degree 0 and T)
// - writing to an ofstream: << operator (writes the Polynomial in a nice format)
//

#include <iostream>
#include <vector>
#include <map>

using namespace std;

namespace Polynomial {
    template<typename T>
    class PolynomialSparse;

    template<typename T>
    class Polynomial {
    private:
        //////////////////////////
        //         DATA         //
        //////////////////////////

        std::vector<T> coeficients;

        //////////////////////////
        //  Internal functions  //
        //////////////////////////

        // Removes not needed leading zero coeficients
        void removeLeadingZeroes() {
            size_t i = 1;
            for (; i <= coeficients.size(); ++i) {
                if (coeficients[coeficients.size() - i] != T(0))
                    break;
            }
            --i;
            if (coeficients.size() - i < coeficients.size() && coeficients.size() - i >= 0) {
                coeficients.resize(coeficients.size() - i);
                coeficients.shrink_to_fit();
            }
        }

        // Main Add/Substract function
        // Sign = 0 is "+="; Sign = 1 is "-="
        void Add(const Polynomial<T>& other, bool sign) {

            if (coeficients.size() < other.coeficients.size())
                coeficients.resize(other.coeficients.size());
            for (size_t i = 0; i < other.coeficients.size(); ++i)
                if (sign)
                    coeficients[i] -= other.coeficients[i];
                else
                    coeficients[i] += other.coeficients[i];
            removeLeadingZeroes();
        }

        // Main Division function
        Polynomial<T> Div(const Polynomial<T>& other) const {
            if (other.size() == 1)
                return (*this / other[0]);
            if (size() == 0)
                return (*this);
            Polynomial<T> res;
            Polynomial<T> temp(*this);
            res.coeficients.resize(size() + other.size() + 1);
            T coef;
            while (temp.size() >= other.size() && temp.size() > 0) {
                coef = temp[temp.size() - 1] / other[other.size() - 1];
                if (coef == T(0))
                    break;
                res.coeficients[temp.size() - other.size()] = coef;
                temp -= coef * other.IncreaseVarPower(temp.size() - other.size());
            }
            res.removeLeadingZeroes();
            return res;
        }

        // Main gcf function
        Polynomial<T> gcd(const Polynomial<T>& lhs, const Polynomial<T>& rhs) const {
            if (rhs == T(0))
                return lhs;
            return gcd(rhs, lhs % rhs);
        }

    public:
        //
        // Constructors
        //

        Polynomial<T>() : coeficients() {}
        Polynomial<T>(const vector<T>& coef) : coeficients(coef) {
            removeLeadingZeroes();
        }
        Polynomial<T>(const T& coef) : coeficients(1) {
            coeficients[0] = coef;
            removeLeadingZeroes();
        }
        template<class It>
        Polynomial<T>(It first, It last) {
            while (first != last) {
                coeficients.push_back(T(*first));
                ++first;
            }
            removeLeadingZeroes();
        }

        //
        // Convert to Sparse
        //
        PolynomialSparse<T> toSparse() const {
            return PolynomialSparse<T>(coeficients);
        }

        //
        // Size
        //

        Polynomial<T>(const Polynomial<T>& other) : coeficients(other.coeficients) {}
        int Degree() const {
            return coeficients.size() - 1;
        }

        size_t size() const {
            return coeficients.size();
        }

        //
        // Data access
        //

        T operator[] (size_t i) const {
            T res = T(0);
            if (i < coeficients.size())
                res = coeficients[i];
            return res;
        }

        vector<T>& GetCoef() {
            return coeficients;
        }

        vector<T> GetCoef() const {
            return coeficients;
        }

        //
        //  Iterators
        //

        auto begin() {
            return coeficients.begin();
        }

        auto end() {
            return coeficients.end();
        }

        auto begin() const {
            return coeficients.cbegin();
        }

        auto end() const {
            return coeficients.cend();
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
            temp.removeLeadingZeroes();
            return temp;
        }

        // Scalar version
        Polynomial<T>& operator += (const T& other) {
            if (coeficients.size())
                coeficients[0] += other;
            else
                coeficients.push_back(other);
            removeLeadingZeroes();
            return *this;
        }

        Polynomial<T>& operator -= (const T& other) {
            if (!coeficients.size())
                coeficients.push_back(T(0));
            coeficients[0] -= other;
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
            std::fill(coeficients.begin(), coeficients.end(), T(0));
            if (coeficients.size() + other.coeficients.size() >= 1)
                coeficients.resize(coeficients.size() + other.coeficients.size() - 1);
            for (size_t i = 0; i < temp.coeficients.size(); ++i) {
                for (size_t j = 0; j < other.coeficients.size(); ++j) {
                    coeficients[i + j] += temp.coeficients[i] * other.coeficients[j];
                }
            }
            removeLeadingZeroes();
            return *this;
        }

        Polynomial<T> operator * (const Polynomial<T>& other) const {
            Polynomial<T> temp;
            if (coeficients.size() + other.coeficients.size() >= 1)
                temp.coeficients.resize(coeficients.size() + other.coeficients.size() - 1);
            for (size_t i = 0; i < coeficients.size(); ++i) {
                for (size_t j = 0; j < other.coeficients.size(); ++j) {
                    temp.coeficients[i + j] += coeficients[i] * other.coeficients[j];
                }
            }
            temp.removeLeadingZeroes();
            return temp;
        }

        // Scalar version
        Polynomial<T>& operator *= (const T& other) {
            if (other == T(0)) {
                coeficients.resize(0);
                coeficients.shrink_to_fit();
            } else {
                for (size_t i = 0; i < coeficients.size(); ++i)
                    coeficients[i] *= other;
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
        //
        // Division
        //
        Polynomial<T> IncreaseVarPower(size_t power) const {
            Polynomial<T> temp(T(0));
            temp.coeficients.resize(size() + power + 1);
            for (size_t i = 0; i < size(); ++i) {
                if (coeficients[i] != T(0))
                    temp.coeficients[i + power] = coeficients[i];
            }
            temp.removeLeadingZeroes();
            return temp;
        }

        Polynomial<T>& operator /= (const Polynomial<T>& other) {
            *this = Div(other);
            removeLeadingZeroes();
            return *this;
        }
        Polynomial<T> operator / (const Polynomial<T>& other) const {
            return Div(other);
        }
        Polynomial<T>& operator /= (const T& other) {
            for (T& el : coeficients)
                el /= other;
            removeLeadingZeroes();
            return *this;
        }
        Polynomial<T> operator / (const T& other) const {
            Polynomial<T> temp(*this);
            return (temp /= other);
        }

        //
        // Mod
        //
        Polynomial<T> operator % (const Polynomial<T>& other) const {
            return (*this - Div(other) * other);
        }

        //
        // Calculate value at given point
        //

        T operator() (const T& val) const {
            T res = T(0);
            T temp = T(1);
            for (size_t i = 0; i < coeficients.size(); ++i) {
                res += temp * coeficients[i];
                temp *= val;
            }
            return res;
        }

        //
        // gcd
        //

        Polynomial<T> operator , (const Polynomial<T>& other) const {
            Polynomial<T> res = gcd(*this, other);
            res.removeLeadingZeroes();
            if (res.size() > 0) {
                for (T& coef : res.coeficients)
                    coef /= res.coeficients[res.size() - 1];
            }
            return res;
        }

        //
        // Composition
        //

        template<typename Z>
        friend Polynomial<Z> operator & (const Polynomial<Z>& lhs, const Polynomial<Z>& rhs);

        //
        // Equality
        //

        template<typename Z>
        friend bool operator == (const Polynomial<Z>& lhs, const Polynomial<Z>& rhs);
        template<typename Z>
        friend bool operator != (const Polynomial<Z>& lhs, const Polynomial<Z>& rhs);
    };

    //
    // Writing to ostream
    //

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
        if (lhs.coeficients.size() > 0)
            initializer = lhs.coeficients[0];
        Polynomial<T> res(initializer);
        Polynomial<T> temp(T(1));
        res.coeficients.resize((lhs.coeficients.size() - 1) * (rhs.coeficients.size() - 1) + 1);
        for (size_t i = 1; i != lhs.coeficients.size(); ++i) {
            temp *= rhs.coeficients;  // raise to the i-th power
            res += lhs.coeficients[i] * temp;
        }
        res.removeLeadingZeroes();
        return res;
    }
    template<typename T>
    bool operator == (const Polynomial<T>& lhs, const Polynomial<T>& rhs) {
        return (lhs.coeficients == rhs.coeficients);
    }
    template<typename T>
    bool operator != (const Polynomial<T>& lhs, const Polynomial<T>& rhs) {
        return (lhs.coeficients != rhs.coeficients);
    }
    template<typename T>
    bool operator == (const Polynomial<T>& lhs, const T& rhs) {
        return ((lhs.size() == 1 && lhs[0] == rhs) || (lhs.size() == 0 && rhs == T(0)));
    }
    template<typename T>
    bool operator != (const Polynomial<T>& lhs, const T& rhs) {
        return (lhs == rhs);
    }

    template<typename T>
    T pow(T val, size_t power, T one = T(1)) {
        T res = one;
        while (power) {
            if (!(power & 1)) {
                val *= T(val);
                power >>= 1;
            } else {
                res *= val;
                --power;
            }
        }
        return res;
    }

    template <typename T>
    class PolynomialSparse {
    private:
        std::map<size_t, T> coefficients;
        // Removes not needed leading zero coefficients
        void removeZeroes() {
            auto it = coefficients.begin();
            while (it != coefficients.end()) {
                if (it->second == T(0)) {
                    auto delIt = it;
                    ++it;
                    coefficients.erase(delIt);
                } else {
                    ++it;
                }
            }
        }
        // Main Add/Substract function
        // Sign = 0 is "+="; Sign = 1 is "-="
        void Add(const PolynomialSparse<T>& other, bool sign) {
            for (auto el = other.coefficients.begin(); el != other.coefficients.end(); ++el) {
                if (sign)
                    coefficients[el->first] -= el->second;
                else
                    coefficients[el->first] += el->second;
            }
            removeZeroes();
        }
        // Main Division function
        PolynomialSparse<T> Div(const PolynomialSparse<T>& other) const {
            if (other.Degree() == 0)
                return (*this / other[0]);
            if (Degree() == -1)
                return (*this);
            vector<T> res(Degree() + other.Degree() + 2);
            PolynomialSparse<T> temp(*this);
            T coef;
            while (temp.Degree() >= other.Degree() && temp.size() > 0) {
                coef = temp.coefficients.rbegin()->second / other.coefficients.rbegin()->second;
                if (coef == T(0))
                    break;
                res[temp.Degree() - other.Degree()] = coef;
                temp -= coef * other.IncreaseVarPower(temp.Degree() - other.Degree());
            }
            return res;
        }
        // Main gcf function
        PolynomialSparse<T> gcd(const PolynomialSparse<T>& lhs, const PolynomialSparse<T>& rhs) const {
            if (rhs == T(0))
                return lhs;
            return gcd(rhs, lhs % rhs);
        }

    public:
        //
        // Constructors
        //
        PolynomialSparse<T>() : coefficients() {}
        PolynomialSparse<T>(const T& coef) : coefficients() {
            if (coef != T(0))
                coefficients[0] = coef;
        }
        PolynomialSparse<T>(const vector<T>& coef) : coefficients() {
            for (size_t i = 0; i < coef.size(); ++i)
                if (coef[i] != T(0))
                    coefficients[i] = coef[i];
            removeZeroes();
        }
        PolynomialSparse<T>(const map<size_t, T>& coef) : coefficients(coef) {}
        template <typename Iter>
        PolynomialSparse<T>(Iter first, Iter last) {
            size_t i = 0;
            while (first != last) {
                coefficients[i] = T(*first);
                ++first;
                ++i;
            }
            removeZeroes();
        }
        PolynomialSparse<T>(const PolynomialSparse<T>& other) : coefficients(other.coefficients) {}
        //
        // Polynomial Conversion
        //
        Polynomial<T> toDense() const {
            vector<T> vec(Degree()+1);
            for (auto it = coefficients.begin(); it != coefficients.end(); ++it) {
                vec[it->first] = it->second;
            }
            return Polynomial<T>(vec);
        }
        //
        // Size
        //
        int Degree() const {
            if (coefficients.size())
                return static_cast<int>(coefficients.rbegin()->first);
            return -1;
        }
        size_t size() const {
            return coefficients.size();
        }
        //
        // Data access
        //
        T operator[] (size_t i) const {
            if (coefficients.count(i))
                return coefficients.at(i);
            return T(0);
        }
        auto& GetCoef() {
            return coefficients;
        }
        auto GetCoef() const {
            return coefficients;
        }
        //
        //  Iterators
        //
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
        // Add / Substract
        PolynomialSparse<T>& operator += (const PolynomialSparse<T>& other) {
            Add(other, false);
            removeZeroes();
            return *this;
        }
        PolynomialSparse<T>& operator -= (const PolynomialSparse<T>& other) {
            Add(other, true);
            removeZeroes();
            return *this;
        }
        PolynomialSparse<T> operator + (const PolynomialSparse<T>& other) const {
            PolynomialSparse<T> temp(*this);
            temp.Add(other, false);
            return temp;
        }
        PolynomialSparse<T> operator - (const PolynomialSparse<T>& other) const {
            PolynomialSparse<T> temp(*this);
            temp.Add(other, true);
            temp.removeZeroes();
            return temp;
        }
        // Scalar version
        PolynomialSparse<T>& operator += (const T& other) {
            if (Degree() < 0)
                coefficients[0] = T(0);
            coefficients[0] += other;
            removeZeroes();
            return *this;
        }
        PolynomialSparse<T>& operator -= (const T& other) {
            if (Degree() < 0)
                coefficients[0] = T(0);
            coefficients[0] -= other;
            removeZeroes();
            return *this;
        }
        PolynomialSparse<T> operator + (const T& other) const {
            PolynomialSparse<T> temp(*this);
            return (temp += other);
        }
        friend PolynomialSparse<T> operator + (const T& other, const PolynomialSparse<T>& pol) {
            PolynomialSparse<T> temp(other);
            temp += pol;
            return temp;
        }
        PolynomialSparse<T> operator - (const T& other) const {
            PolynomialSparse<T> temp(*this);
            return temp -= other;
        }
        friend PolynomialSparse<T> operator - (const T& other, const PolynomialSparse<T>& pol) {
            PolynomialSparse<T> temp(other);
            temp -= pol;
            return temp;
        }
        // Multiplication
        PolynomialSparse<T>& operator *= (const PolynomialSparse<T>& other) {
            PolynomialSparse<T> temp(*this);
            coefficients.clear();
            for (auto i = temp.coefficients.begin(); i != temp.coefficients.end(); ++i) {
                for (auto j = other.coefficients.begin(); j != other.coefficients.end(); ++j) {
                    coefficients[i->first + j->first] += (i->second) * (j->second);
                }
            }
            removeZeroes();
            return *this;
        }
        PolynomialSparse<T> operator * (const PolynomialSparse<T>& other) const {
            PolynomialSparse<T> temp;
            for (auto i = coefficients.begin(); i != coefficients.end(); ++i) {
                for (auto j = other.coefficients.begin(); j != other.coefficients.end(); ++j) {
                    temp.coefficients[i->first + j->first] += (i->second) * (j->second);
                }
            }
            temp.removeZeroes();
            return temp;
        }
        // Scalar version
        PolynomialSparse<T>& operator *= (const T& other) {
            if (other == T(0)) {
                coefficients.clear();
            } else {
                for (auto it = coefficients.begin(); it != coefficients.end(); ++it)
                    coefficients[it->first] *= other;
            }
            removeZeroes();
            return *this;
        }
        PolynomialSparse<T> operator * (const T& other) const {
            PolynomialSparse<T> temp(*this);
            return (temp *= other);
        }
        friend PolynomialSparse<T> operator * (const T& other, const PolynomialSparse<T>& pol) {
            PolynomialSparse<T> temp(other);
            return (temp *= pol);
        }
        //
        // Division
        //
        // Equivalent to multiplying Polynomial by x^power;
        PolynomialSparse<T> IncreaseVarPower(size_t power) const {
            PolynomialSparse<T> temp;
            for (auto it = coefficients.begin(); it != coefficients.end(); ++it)
                temp.coefficients[it->first + power] = it->second;
            return temp;
        }
        PolynomialSparse<T>& operator /= (const PolynomialSparse<T>& other) {
            *this = Div(other);
            removeZeroes();
            return *this;
        }
        PolynomialSparse<T> operator / (const PolynomialSparse<T>& other) const {
            return Div(other);
        }
        PolynomialSparse<T>& operator /= (const T& other) {
            for (auto it = coefficients.begin(); it != coefficients.end(); ++it)
                it->second /= other;
            removeZeroes();
            return *this;
        }
        PolynomialSparse<T> operator / (const T& other) const {
            PolynomialSparse<T> temp(*this);
            return (temp /= other);
        }
        //
        // Mod
        //
        PolynomialSparse<T> operator % (const PolynomialSparse<T>& other) const {
            return (*this - Div(other) * other);
        }
        //
        // Calculate value at given point
        //
        // Returns value of the Polynomial with variable x = val;
        T operator() (const T& val) const {
            T res = T(0);
            for (auto it = coefficients.begin(); it != coefficients.end(); ++it) {
                res += pow(val, it->first) * (it->second);
            }
            return res;
        }

        //
        // gcf
        //
        // GCF of two Polynomials [recommended synthax : (a, b)]
        PolynomialSparse<T> operator , (const PolynomialSparse<T>& other) const {
            PolynomialSparse<T> res = gcd(*this, other);
            res.removeZeroes();
            if (res.Degree() > -1) {
                for (auto it = res.coefficients.begin(); it != res.coefficients.end(); ++it)
                    res.coefficients[it->first] /= res.coefficients.rbegin()->second;
            }
            return res;
        }
        //
        // Composition
        //
        template <typename Z>
        friend PolynomialSparse<Z> operator & (const PolynomialSparse<Z>& lhs, const PolynomialSparse<Z>& rhs);
        //
        // Equality
        //
        template <typename Z>
        friend bool operator == (const PolynomialSparse<Z>& lhs, const PolynomialSparse<Z>& rhs);
        template <typename Z>
        friend bool operator != (const PolynomialSparse<Z>& lhs, const PolynomialSparse<Z>& rhs);
        template <typename Z>
        friend std::ostream& operator << (std::ostream& out, const PolynomialSparse<Z>& pol);
    };
    //
    // Writing to ostream
    //
    template <typename T>
    std::ostream& operator << (std::ostream& out, const PolynomialSparse<T>& pol) {
        if (pol.size()) {
            for (auto it = pol.coefficients.rbegin(); it != pol.coefficients.rend(); ++it) {
                if (it != pol.coefficients.rbegin() && it->second > 0)
                    out << "+";
                if ((it->second != T(1) && it->second != T(-1) && it->first > 0) || it->first == 0) {
                    out << it->second;
                    if (it->first > 0)
                        out << "*";
                }
                if (it->first > 0) {
                    if (it->second != T(-1))
                        out << "x";
                    else
                        out << "-x";
                    if (it->first > 1)
                        out << "^" << it->first;
                }
            }
        } else {
            out << T(0);
        }
        return out;
    }
    template <typename T>
    // Returns the composition of lhs and rhs Polynomials
    PolynomialSparse<T> operator & (const PolynomialSparse<T>& lhs, const PolynomialSparse<T>& rhs) {
        PolynomialSparse<T> res(T(0));
        for (auto it = lhs.coefficients.begin(); it != lhs.coefficients.end(); ++it) {
            res += it->second * pow(rhs, it->first, PolynomialSparse<T>(T(1)));
        }
        res.removeZeroes();
        return res;
    }
    template <typename T>
    bool operator == (const PolynomialSparse<T>& lhs, const PolynomialSparse<T>& rhs) {
        return (lhs.coefficients == rhs.coefficients);
    }
    template <typename T>
    bool operator != (const PolynomialSparse<T>& lhs, const PolynomialSparse<T>& rhs) {
        return (lhs.coefficients != rhs.coefficients);
    }
    template <typename T>
    bool operator == (const PolynomialSparse<T>& lhs, const T& rhs) {
        return ((lhs.Degree() == 0 && lhs[0] == rhs) || (lhs.Degree() == -1 && rhs == T(0)));
    }
    template <typename T>
    bool operator != (const PolynomialSparse<T>& lhs, const T& rhs) {
        return !(lhs == rhs);
    }
}






//  Version 0.2.2
//  (C) TechnoVirus / temporIllustribus, February 2019

//                                             
//                  ...'',;::clo,              
//         .:lodxkO0KKXNXK00Oxdl.     ..       
//         ;XWNK0OxdlxNO,..         .lK0c.     
//         .,,..     cXc          .c0NO:.      
//                   dO.         ;OXk;         
//                  'Oo        'xKx;           
//                  :O,      .o0x,             
//                  od    .'ckd'               
//                 .d;   :dxl.                 
//                 ;o. ,lc.                    
//                 ::'::.                      
//                .lc'.                        
//                .'.                          
//      