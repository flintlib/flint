.. _pyflint:

**flint_ctypes** - Python interface
===============================================================================

There is a Python wrapper (``flint_ctypes``) included with FLINT
available in the ``src/python`` directory. This wrapper is not currently
officially supported and should not be used in production, but it can be
useful for experimenting with FLINT.

.. highlight:: python3

Introduction
-------------------------------------------------------------------------------

Examples::

    >>> from flint_ctypes import *
    >>> QQ.bernoulli(50)
    495057205241079648212477525/66
    >>> sign, primes, exponents = _.factor()
    >>> sign
    1
    >>> primes
    [5, 417202699, 47464429777438199, 2, 3, 11]
    >>> exponents
    [2, 1, 1, -1, -1, -1]
    >>> sign * (primes ** exponents).product()
    495057205241079648212477525/66

Types, parents and coercions
...............................................................................

    >>> ZZ(5)
    5
    >>> _.parent()
    Integer ring (fmpz)
    >>> QQ(5)
    5
    >>> _.parent()
    Rational field (fmpq)
    >>> ZZ(10) / ZZ(6)
    Traceback (most recent call last):
      ...
    FlintDomainError: x / y is not an element of {Integer ring (fmpz)} for {x = 10}, {y = 6}
    >>> x = QQ(1) / 2; x ** x
    Traceback (most recent call last):
      ...
    FlintDomainError: x ** y is not an element of {Rational field (fmpq)} for {x = 1/2}, {y = 1/2}

    >>> ZZ(10) / QQ(6)
    5/3
    >>> x = QQbar(1) / 2; x ** x
    Root a = 0.707107 of 2*a^2-1

Real and complex numbers
...............................................................................

    >>> RR.zeta(2)
    [1.644934066848226 +/- 4.57e-16]
    >>> RR.prec = 128
    >>> RR.zeta(2)
    [1.64493406684822643647241516664602518922 +/- 2.88e-39]
    >>> RR.prec = 53       # restore default

API documentation
-------------------------------------------------------------------------------

..
   automodsumm:: flint_ctypes

..
   autoclass:: flint_ctypes.FlintDomainError

..
   autoclass:: flint_ctypes.FlintUnableError

..
   automodule:: flint_ctypes
..
   :members:
..
   :undoc-members:
..
   :special-members: __init__ , __bool__
..
   :member-order: bysource


.. raw:: latex

    \newpage
