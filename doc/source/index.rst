Calcium
===================================

.. only:: html

    .. image:: _static/ca2.svg
        :align: center

**Calcium** (pronounced "kalkium") is a C library for exact computation
with real and complex numbers (both algebraic and transcendental),
presently in early development.

Calcium is free software (LGPL). It depends on
`GMP <https://gmplib.org/>`_, `MPFR <https://mpfr.org/>`_,
`Flint <http://flintlib.org/>`_, `Antic <https://github.com/wbhart/antic/>`_
and `Arb <http://arblib.org/>`_.

Source code: https://github.com/fredrik-johansson/calcium

Module documentation
---------------------------

.. toctree::
   :maxdepth: 2

   calcium.rst
   ca.rst
   fmpz_mpoly_q.rst
   qqbar.rst


FAQ
---------------------------

**I thought exact calculation with real numbers isn't possible. Isn't x = 0 undecidable?**

In general, yes. In practice, much of calculus
and number theory only depends on numbers that are
simple combinations of well-known elementary
and special functions, and there are heuristics that work quite well
for deciding predicates about such numbers.
Calcium will be able to give a definitive answer at least in 
simple cases (for example, proving 
`16 \operatorname{atan}(\tfrac{1}{5}) - 4 \operatorname{atan}(\tfrac{1}{239}) = \pi`
or `\sqrt{5+2\sqrt{6}} = \sqrt{2}+\sqrt{3}`),
and will simply answer "Unknown" when its heuristics are not powerful enough.

**Isn't it going to be horribly slow?**

It will definitely be too slow to replace floating-point numbers
for 99.9% of scientific computing. The target is symbolic and algebraic computation.
A big factor in making this kind of library practical is that
Flint now has extremely fast multivariate polynomial arithmetic
which can be used to implement multivariate transcendental number fields
with adequate performance.

Calcium will generally be much slower than arbitrary-precision ball arithmetic.
It will often make sense for users to first try a numerical evaluation with Arb,
and fall back on an exact calculation with Calcium only if that
fails (typically because an exact comparison is needed).

**Why this C library? Why not use a computer algebra system?**

Calculating with constant values is only a small part of what
a computer algebra system has to do, but it is actually one
of the most complex parts.
Calcium is intended to take some of the complexity out of the task
of building a computer algebra system by offering a black-box
solution to the "constant problem".
This solution will not be perfect, but it will be "good enough"
for many applications.
C is a pain, but it is the best choice for building a library that will
not be tied to a particular computer algebra system
or programming language ecosystem.
Initially, the goal of Calcium will be to implement high-performance
basic data structures in C; this will enable future experimentation with algorithms
in both C and in high-level languages using wrappers.

Indices and tables
==================

* :ref:`genindex`



