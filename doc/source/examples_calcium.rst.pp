.. _examples-calcium:

Calcium example programs
===============================================================================

.. highlight:: text

See :ref:`examples` for general information about example programs.
Running::

    make examples

will compile the programs and place the binaries in
``build/examples``. The examples related to the Calcium module are
documented below.

elementary.c
-------------------------------------------------------------------------------

This program evaluates several elementary expressions.
For some inputs,
Calcium's arithmetic should produce
a simplified result automatically.
Some inputs do not yet automatically simplify as much
as one might hope.
Calcium may still able to prove that such a number is zero or nonzero;
the output of :func:`ca_check_is_zero` is then ``T_TRUE`` or ``T_FALSE``.

Sample output::

    > build/examples/elementary
    >>> Exp(Pi*I) + 1
    0

    >>> Log(-1) / (Pi*I)
    1

    >>> Log(-I) / (Pi*I)
    -0.500000 {-1/2}

    >>> Log(1 / 10^123) / Log(100)
    -61.5000 {-123/2}

    >>> Log(1 + Sqrt(2)) / Log(3 + 2*Sqrt(2))
    0.500000 {1/2}

    >>> Sqrt(2)*Sqrt(3) - Sqrt(6)
    0

    >>> Exp(1+Sqrt(2)) * Exp(1-Sqrt(2)) / (Exp(1)^2)
    1

    >>> I^I - Exp(-Pi/2)
    0

    >>> Exp(Sqrt(3))^2 - Exp(Sqrt(12))
    0

    >>> 2*Log(Pi*I) - 4*Log(Sqrt(Pi)) - Pi*I
    0

    >>> -I*Pi/8*Log(2/3-2*I/3)^2 + I*Pi/8*Log(2/3+2*I/3)^2 + Pi^2/12*Log(-1-I) + Pi^2/12*Log(-1+I) + Pi^2/12*Log(1/3-I/3) + Pi^2/12*Log(1/3+I/3) - Pi^2/48*Log(18)
    0

    >>> Sqrt(5 + 2*Sqrt(6)) - Sqrt(2) - Sqrt(3)
    0e-1126 {a-c-d where a = 3.14626 [Sqrt(9.89898 {2*b+5})], b = 2.44949 [b^2-6=0], c = 1.73205 [c^2-3=0], d = 1.41421 [d^2-2=0]}
    >>> Is zero?
    T_TRUE

    >>> Sqrt(I) - (1+I)/Sqrt(2)
    0e-1126 + 0e-1126*I {(2*a-b*c-b)/2 where a = 0.707107 + 0.707107*I [Sqrt(1.00000*I {c})], b = 1.41421 [b^2-2=0], c = I [c^2+1=0]}
    >>> Is zero?
    T_TRUE

    >>> Exp(Pi*Sqrt(163)) - (640320^3 + 744)
    -7.49927e-13 {a-262537412640768744 where a = 2.62537e+17 [Exp(40.1092 {b*c})], b = 3.14159 [Pi], c = 12.7671 [c^2-163=0]}

    >>> Erf(2*Log(Sqrt(1/2-Sqrt(2)/4))+Log(4)) - Erf(Log(2-Sqrt(2)))
    0


    cpu/wall(s): 0.022 0.022
    virt/peak/res/peak(MB): 36.45 36.47 9.37 9.37


binet.c
-------------------------------------------------------------------------------

This program computes the *n*-th Fibonacci number using Binet's formula
`F_n = (\varphi^n - (1-\varphi)^n)/\sqrt{5}` where
`\varphi = \tfrac{1}{2} (1+\sqrt{5})`. The program takes *n* as input.

Sample output::

    > build/examples/binet 250
    7.89633e+51 {7896325826131730509282738943634332893686268675876375}

    cpu/wall(s): 0.002 0.001
    virt/peak/res/peak(MB): 36.14 36.14 5.81 5.81

This illustrates exact arithmetic in algebraic number fields.
The program also illustrates another aspect of Calcium arithmetic:
evaluation limits. For example, trying
to compute the index `n = 10^6`
Fibonacci number hits an evaluation limit, so the value is
not expanded to an explicit integer::

    > build/examples/binet 1000000
    1.95328e+208987 {(a*c-b*c)/5 where a = 4.36767e+208987 [Pow(1.61803 {(c+1)/2}, 1.00000e+6 {1000000})], b = 2.28955e-208988 [Pow(-0.618034 {(-c+1)/2}, 1.00000e+6 {1000000})], c = 2.23607 [c^2-5=0]}

    cpu/wall(s): 0.006 0.005
    virt/peak/res/peak(MB): 36.14 36.14 9.05 9.05

Calling the program with ``-limit B n`` raises the bit evaluation
limit to *B*. Setting this large enough allows `F_{10^6}` to expand
to an integer (the following output has been truncated to avoid
reproducing all 208988 digits)::

    > build/examples/binet -limit 10000000 1000000
    1.95328e+208987 {1953282128...8242546875}

    cpu/wall(s): 0.229 0.242
    virt/peak/res/peak(MB): 36.79 37.29 7.13 7.13

The exact mechanisms and interfaces for evaluation limits are still a
work in progress.

machin.c
-------------------------------------------------------------------------------

This program checks several variations of Machin's formula

.. math::

    \frac{\pi}{4} = 4 \operatorname{atan}\left(\frac{1}{5}\right) - \operatorname{atan}\left(\frac{1}{239}\right)

expressing `\pi` or logarithms of small integers in terms of
arctangents or hyperbolic arctangents of rational numbers.
The program actually evaluates
`4 \operatorname{atan}\left(\tfrac{1}{5}\right) - \operatorname{atan}\left(\tfrac{1}{239}\right) - \tfrac{\pi}{4}`
(etc.) and prints the result, which should be precisely 0, proving the identity.
Inverse trigonometric functions are not yet implemented in Calcium,
so the example program evaluates them using logarithms.

Sample output::

    > build/examples/machin
    [(1)*Atan(1/1) - Pi/4]   =   0
    [(1)*Atan(1/2) + (1)*Atan(1/3) - Pi/4]   =   0
    [(2)*Atan(1/2) + (-1)*Atan(1/7) - Pi/4]   =   0
    [(2)*Atan(1/3) + (1)*Atan(1/7) - Pi/4]   =   0
    [(4)*Atan(1/5) + (-1)*Atan(1/239) - Pi/4]   =   0
    [(1)*Atan(1/2) + (1)*Atan(1/5) + (1)*Atan(1/8) - Pi/4]   =   0
    [(1)*Atan(1/3) + (1)*Atan(1/4) + (1)*Atan(1/7) + (1)*Atan(1/13) - Pi/4]   =   0
    [(12)*Atan(1/49) + (32)*Atan(1/57) + (-5)*Atan(1/239) + (12)*Atan(1/110443) - Pi/4]   =   0

    [(14)*Atanh(1/31) + (10)*Atanh(1/49) + (6)*Atanh(1/161) - Log(2)]   =   0
    [(22)*Atanh(1/31) + (16)*Atanh(1/49) + (10)*Atanh(1/161) - Log(3)]   =   0
    [(32)*Atanh(1/31) + (24)*Atanh(1/49) + (14)*Atanh(1/161) - Log(5)]   =   0
    [(144)*Atanh(1/251) + (54)*Atanh(1/449) + (-38)*Atanh(1/4801) + (62)*Atanh(1/8749) - Log(2)]   =   0
    [(228)*Atanh(1/251) + (86)*Atanh(1/449) + (-60)*Atanh(1/4801) + (98)*Atanh(1/8749) - Log(3)]   =   0
    [(334)*Atanh(1/251) + (126)*Atanh(1/449) + (-88)*Atanh(1/4801) + (144)*Atanh(1/8749) - Log(5)]   =   0
    [(404)*Atanh(1/251) + (152)*Atanh(1/449) + (-106)*Atanh(1/4801) + (174)*Atanh(1/8749) - Log(7)]   =   0

    cpu/wall(s): 0.016 0.016
    virt/peak/res/peak(MB): 35.57 35.57 8.80 8.80

swinnerton_dyer_poly.c
-------------------------------------------------------------------------------

This program computes the coefficients of the Swinnerton-Dyer polynomial

.. math::

    S_n = \prod (x \pm \sqrt{2} \pm \sqrt{3} \pm \sqrt{5} \pm \ldots \pm \sqrt{p_n})

where `p_n` denotes the `n`-th prime number and all combinations
of signs are taken. This polynomial has degree `2^n`.
The polynomial is expanded from its roots
using naive polynomial multiplication over :type:`ca_t` coefficients.
There are far more efficient ways to construct this polynomial;
this program simply illustrates that arithmetic in
multivariate number fields works smoothly.

The program prints the coefficients of `S_n`, from the constant
term to the coefficient of `x^{2^n}`.

Sample output::

    > build/examples/swinnerton_dyer_poly 3
    576
    0
    -960
    0
    352
    0
    -40
    0
    1

    cpu/wall(s): 0.002 0.002
    virt/peak/res/peak(MB): 35.07 35.11 5.40 5.40

A big benchmark problem (output truncated)::

    > build/examples/swinnerton_dyer_poly 10
    4.35675e+809 {43567450015...212890625}
    0
    ...
    0
    1

    cpu/wall(s): 9.296 9.307
    virt/peak/res/peak(MB): 38.95 38.95 10.01 10.01

huge_expr.c
-------------------------------------------------------------------------------

This program proves equality of two complicated algebraic numbers.
More precisely, the program verifies
that `N = -(1 - |M|^2)^2` where *N* and *M* are given by huge symbolic
expressions involving nested square roots (about 7000
operations in total).

By default, the program runs the computation using :type:`qqbar_t` arithmetic::

    > build/examples/huge_expr
    Evaluating N...
    cpu/wall(s): 7.205 7.206
    Evaluating M...
    cpu/wall(s): 0.933 0.934
    Evaluating E = -(1-|M|^2)^2...
    cpu/wall(s): 0.391 0.391
    N ~ -0.16190853053311203695842869991458578203473645660641
    E ~ -0.16190853053311203695842869991458578203473645660641
    Testing E = N...
    cpu/wall(s): 0.001 0

    Equal = T_TRUE

    Total: cpu/wall(s): 8.53 8.531
    virt/peak/res/peak(MB): 54.50 64.56 24.64 34.61

To run the computation using :type:`ca_t` arithmetic instead,
pass the -ca flag::

    > build/examples/huge_expr -ca
    Evaluating N...
    cpu/wall(s): 0.193 0.193
    Evaluating M...
    cpu/wall(s): 0.024 0.024
    Evaluating E = -(1-|M|^2)^2...
    cpu/wall(s): 0.008 0.009
    N ~ -0.16190853053311203695842869991458578203473645660641
    E ~ -0.16190853053311203695842869991458578203473645660641
    Testing E = N...
    cpu/wall(s): 8.017 8.019

    Equal = T_TRUE

    Total: cpu/wall(s): 8.243 8.246
    virt/peak/res/peak(MB): 61.67 65.29 33.97 37.54

This simplification problem was posted in a help request for Sage
(https://ask.sagemath.org/question/52653).
The C code has been generated from the symbolic expressions
using a Python script.


hilbert_matrix.c
-------------------------------------------------------------------------------

This program constructs the Hilbert matrix
`H_n = (1/(i+j-1))_{i=1,j=1}^n`, computes its
eigenvalues `\lambda_1, \ldots, \lambda_n`,
as exact algebraic numbers, and verifies
the exact trace and determinant formulas

.. math::

    \lambda_1 + \lambda_2 + \ldots + \lambda_n = \operatorname{tr}(H_n), \quad
    \lambda_1 \lambda_2 \cdots \lambda_n = \operatorname{det}(H_n).

Sample output::

    > build/examples/hilbert_matrix 6
    Trace:
    1.87821 {6508/3465}
    1.87821 {6508/3465}
    Equal: T_TRUE

    Det:
    5.36730e-18 {1/186313420339200000}
    5.36730e-18 {1/186313420339200000}
    Equal: T_TRUE


    cpu/wall(s): 0.07 0.069
    virt/peak/res/peak(MB): 36.56 36.66 9.69 9.69

The program accepts the following optional arguments:

* With ``-vieta``, force use of Vieta's formula internally (by default, Calcium
  uses Vieta's formulas when working with algebraic conjugates,
  but only up to some bound on the degree).
* With ``-novieta``, force Calcium not to use Vieta's formulas internally.
* With ``-qqbar``, do a similar computation using :type:`qqbar_t`
  arithmetic.

dft.c
-------------------------------------------------------------------------------

This program demonstrates the
discrete Fourier transform (DFT) in exact arithmetic.
For the input vector `\textbf{x} = (x_n)_{n=0}^{N-1}`, it verifies
the identity

.. math::

    \textbf{x} - \operatorname{DFT}^{-1}(\operatorname{DFT}(\textbf{x})) = 0

where

.. math::

    \operatorname{DFT}(\textbf{x})_n = \sum_{k=0}^{N-1} \omega^{-kn} x_k, \quad
    \operatorname{DFT}^{-1}(\textbf{x})_n = \frac{1}{N} \sum_{k=0}^{N-1} \omega^{kn} x_k,
    \quad \omega = e^{2 \pi i / N}.

The program computes the DFT by naive `O(N^2)` summation (not using FFT).
It uses repeated multiplication of `\omega`
to precompute an array of roots of unity
`1,\omega,\omega^2,\ldots,\omega^{2N-1}`
for use in both the DFT and the inverse DFT.

Usage::

    build/examples/dft [-verbose] [-input i] [-limit B] [-timing T] N

The required parameter ``N`` selects the length of the vector.

The optional flag ``-verbose`` chooses whether to print the arrays.

The optional parameter ``-timing T`` selects a timing method (default = 0).

* 0: run the computation once and time it
* 1: run the computation repeatedly if needed to get an accurate timing, creating a new context object for each iteration so that fields are not cached
* 2: run the computation once, then run the computation at least one more time (repeatedly if needed to get an accurate timing), recycling the same context object to measure the performance with cached fields

The optional parameter ``-input i`` selects an input sequence (default = 0).

* 0: `x_n = n+2`
* 1: `x_n = \sqrt{n+2}`
* 2: `x_n = \log(n+2)`
* 3: `x_n = e^{2 \pi i / (n+2)}`

The optional parameter ``-limit B`` sets the internal degree limit for algebraic numbers.

Sample output::

    > build/examples/dft 4 -input 1 -verbose
    DFT benchmark, length N = 4

    [x] =
    1.41421 {a where a = 1.41421 [a^2-2=0]}
    1.73205 {a where a = 1.73205 [a^2-3=0]}
    2
    2.23607 {a where a = 2.23607 [a^2-5=0]}

    DFT([x]) =
    7.38233 {a+b+c+2 where a = 2.23607 [a^2-5=0], b = 1.73205 [b^2-3=0], c = 1.41421 [c^2-2=0]}
    -0.585786 + 0.504017*I {a*d-b*d+c-2 where a = 2.23607 [a^2-5=0], b = 1.73205 [b^2-3=0], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}
    -0.553905 {-a-b+c+2 where a = 2.23607 [a^2-5=0], b = 1.73205 [b^2-3=0], c = 1.41421 [c^2-2=0]}
    -0.585786 - 0.504017*I {-a*d+b*d+c-2 where a = 2.23607 [a^2-5=0], b = 1.73205 [b^2-3=0], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}

    IDFT(DFT([x])) =
    1.41421 {c where a = 2.23607 [a^2-5=0], b = 1.73205 [b^2-3=0], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}
    1.73205 {b where a = 2.23607 [a^2-5=0], b = 1.73205 [b^2-3=0], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}
    2
    2.23607 {a where a = 2.23607 [a^2-5=0], b = 1.73205 [b^2-3=0], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}

    [x] - IDFT(DFT([x])) =
    0       (= 0   T_TRUE)
    0       (= 0   T_TRUE)
    0       (= 0   T_TRUE)
    0       (= 0   T_TRUE)

    cpu/wall(s): 0.009 0.009
    virt/peak/res/peak(MB): 36.28 36.28 9.14 9.14




.. raw:: latex

    \newpage

