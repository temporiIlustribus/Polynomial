# Polynomial
          ______     _                             _       _ 
          | ___ \   | |                           (_)     | |
          | |_/ ___ | |_   _ _ __   ___  _ __ ___  _  __ _| |
          |  __/ _ \| | | | | '_ \ / _ \| '_ ` _ \| |/ _` | |
          | | | (_) | | |_| | | | | (_) | | | | | | | (_| | |
          \_|  \___/|_|\__, |_| |_|\___/|_| |_| |_|_|\__,_|_|
                        __/ |                                
                       |___/                                 

Simple Polynomial Template class mainly for integer coefficient Polynomials.
Coefficients are stored in a dense format as a vector of T.

Functionality:
   • constructors: from single coeficient T, vector of coefficients or from first-last iterators.
   • internal data size: size() and Degree() (0 Polynomial has Degree of -1)
   • data access: [] operator (only onst version), begin(), end() iterators for coeficient vector
   • math operators: +; -; *, +=; -=; *=;
   • polynomial composition: & operator;
   • equality checks: operators == and != (also works for comparing Polynomials of Degree 0 and T)
   • writing to an ofstream: << operator (writes the Polynomial in a nice format)

