# Polynomial

Description
-----------
Simple Polynomial Template class mainly for integer coefficient Polynomials.
Coefficients are stored in a dense format as a vector of T.

Functionality:
--------------
      • Constructors: from single coeficient T, vector of coefficients or from first-last iterators.
      • Internal data info: size() and Degree() (0 Polynomial has Degree of -1)
      • Data access: [] operator (only onst version), begin(), end() iterators for coeficient vector
      • Math operators: +; -; *, +=; -=; *=;
      • Polynomial composition: & operator;
      • Equality checks: operators == and != (also works for comparing Polynomials of Degree 0 and T)
      • Writing to an ofstream: << operator (writes the Polynomial in a nice format)

