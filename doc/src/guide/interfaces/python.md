---
title: Python interfaces
summary: Using FLINT from Python, with the separate python-flint package or the experimental flint_ctypes wrapper included with FLINT.
order: 1
---

There are two ways to use FLINT from Python.

- *python-flint* is a separate project (https://github.com/flintlib/python-flint)
  maintained alongside FLINT, with Cython bindings to the main types:
  integers, rationals, integers modulo $n$, finite fields, polynomials,
  matrices, power series, and the ball arithmetic types `arb` and `acb`
  with many special functions. It is installed with
  `pip install python-flint` and is the recommended interface for
  production use.
- *flint_ctypes* is a pure Python wrapper included with FLINT in the
  directory `src/python`. It accesses the library through `ctypes` and
  covers a large part of FLINT through the generic rings interface. It is
  not officially supported and should not be used in production, but it
  is convenient for experimenting with FLINT, including functionality
  that python-flint does not wrap yet. It is described in the rest of
  this section.

Other languages have their own interfaces: for example Nemo and
Oscar.jl use FLINT from Julia, and SageMath uses it from Python.

## Setting up flint_ctypes {#setup}

`flint_ctypes.py` loads the shared library `libflint.so`
(`libflint.dylib` on macOS) with the default search path of the dynamic
linker, so FLINT must be built as a shared library and be installed in a
standard location or be on `LD_LIBRARY_PATH`. No compilation is
needed:

```
cd src/python
LD_LIBRARY_PATH=/path/to/flint/lib python3
>>> from flint_ctypes import *
```

`make check PYTHON=1` runs the doctests of `flint_ctypes.py` against the
library that was just built. The directory also contains `flint2.py`
and the modules `RRmath.py` and `CCmath.py`, which export the functions
of the real and complex fields as module-level functions in the style of
Python's `math` and `cmath`.

## Examples {#examples}

Elements belong to parent rings such as `ZZ` (the integers, `fmpz`),
`QQ` (`fmpq`), `RR` and `CC` (real and complex balls, `arb` and `acb`),
`QQbar` (algebraic numbers, `qqbar`) and `ZZx` (polynomials). Methods on a
ring compute functions with values in that ring:

```
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
```

Operations are performed in the parent of the operands, with coercion
between parents where it makes sense. An operation whose result is not in
the parent raises `FlintDomainError`, and one that cannot be decided or
computed raises `FlintUnableError`:

```
>>> ZZ(5).parent()
Integer ring (fmpz)
>>> ZZ(10) / ZZ(6)
Traceback (most recent call last):
  ...
FlintDomainError: x / y is not an element of {Integer ring (fmpz)} for {x = 10}, {y = 6}
>>> ZZ(10) / QQ(6)
5/3
>>> x = QQ(1) / 2; x ** x
Traceback (most recent call last):
  ...
FlintDomainError: x ** y is not an element of {Rational field (fmpq)} for {x = 1/2}, {y = 1/2}
>>> x = QQbar(1) / 2; x ** x
Root a = 0.707107 of 2*a^2-1
```

The precision of the real and complex fields is an attribute of the
parent:

```
>>> RR.zeta(2)
[1.644934066848226 +/- 4.57e-16]
>>> RR.prec = 128
>>> RR.zeta(2)
[1.64493406684822643647241516664602518922 +/- 2.88e-39]
>>> RR.prec = 53       # restore the default
```

The docstrings of `flint_ctypes.py` contain many more examples, for
polynomials, matrices, power series, finite fields, number fields,
$p$-adic numbers, and other types.
