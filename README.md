# Calcium

Calcium (pronounced “kalkium”) is a C library for exact
computation with real and complex numbers, presently in early development.

![calcium logo](http://fredrikj.net/calcium/_images/ca2.svg)

Documentation: http://fredrikj.net/calcium/

Try Online: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fredrik-johansson/calcium/HEAD?filepath=doc%2Fintroduction.ipynb)

Author: Fredrik Johansson <fredrik.johansson@gmail.com>

Features:

* Exact real and complex numbers represented as elements of automatically extended multivariate fields
* Support for algebraic, transcendental and mixed fields
* Automatic, rigorous numerical embeddings and arbitrary-precision numerical evaluation (on top of Arb)
* Efficient field arithmetic (on top of Flint and Antic)
* Automatic, rigorous simplification (using integer relations, ideal reduction, and other methods)
* Complete decision procedures for algebraic numbers
* Partial decision procedures for transcendental numbers
* Polynomials and matrices with exact coefficients
* Exact real and complex algebraic numbers (absolute minpoly representation)
* Multivariate rational functions (on top of Flint)
* Gröbner basis computation (on top of Flint)
* Symbolic expressions (conversions, evaluation, LaTeX output)
