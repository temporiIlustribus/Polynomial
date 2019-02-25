# Polynomial

Description
-----------
Simple Polynomial Template class mainly for integer coefficient Polynomials.
Coefficients are stored in a dense format as a vector of T.

Appologies for the class being written not as library but rather as a normal VS console app.
If you wish to use it as is, it is recommended you either copy the code to your project or make 
a library and paste the class code there. 

Functionality:
--------------
      • Constructors: from single coeficient T, vector of coefficients or from first-last iterators.
      • Internal data info: size() and Degree() (0 Polynomial has Degree of -1)
      • Data access: [] operator (only onst version), begin(), end() iterators for coeficient vector
      • Math operators: +; -; *, +=; -=; *=;
      • Polynomial composition: & operator;
      • Equality checks: operators == and != (also works for comparing Polynomials of Degree 0 and T)
      • Writing to an ofstream: << operator (writes the Polynomial in a nice format)

