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

// Simple Polynomial Template class mainly for integer coefficient Polynomials.
// Coefficients are stored in a dense format as a vector of T.
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

using namespace std;

namespace Polynomial {
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
}






//  Version 0.2.1
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