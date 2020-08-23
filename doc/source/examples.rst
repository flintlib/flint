.. _examples:

Example programs
===============================================================================

.. highlight:: text

The *examples* directory
(https://github.com/fredrik-johansson/calcium/tree/master/examples)
contains complete C programs illustrating use of Calcium.
Running::

    make examples

will compile the programs and place the binaries in ``build/examples``.

elementary.c
-------------------------------------------------------------------------------

This program evaluates several elementary expressions.
For some inputs,
Calcium's arithmetic should produce
a simplified result automatically (e.g. ``p/q  in  QQ`` if the
value is a rational number `p/q`).
Some inputs do not yet automatically simplify as much
as one might hope. 
Calcium may still able to prove that such a number is zero or nonzero;
the output of :func:`ca_check_is_zero` is then ``T_TRUE`` or ``T_FALSE``.

Sample output::

    > build/examples/elementary 
    [Exp(Pi*I) + 1]   =   0  in  QQ

    [Log(-1) / (Pi*I)]   =   1  in  QQ

    [Log(-I) / (Pi*I)]   =   -1/2  in  QQ

    [Log(1 / 10^123) / Log(100)]   =   -123/2  in  QQ

    [Log(1 + Sqrt(2)) / Log(3 + 2*Sqrt(2))]   =   1/2  in  QQ

    [Sqrt(2)*Sqrt(3) - Sqrt(6)]   =   0  in  QQ

    [Exp(1+Sqrt(2)) * Exp(1-Sqrt(2)) / (Exp(1)^2)]   =   1  in  QQ

    [I^I - Exp(-Pi/2)]   =   0  in  QQ

    [Exp(Sqrt(3))^2 - Exp(Sqrt(12))]   =   0  in  QQ

    [Sqrt(5 + 2*Sqrt(6)) - Sqrt(2) - Sqrt(3)]   =
    (x1-x3-x4)/(1)  in  QQ(x1, x2, x3, x4) where {x1 = Sqrt(2*x1+5  in  QQ(x1) where {x1 = Algebraic 2.4494897}), x2 = Algebraic 2.4494897, x3 = Algebraic 1.7320508, x4 = Algebraic 1.4142136} with ideal {x1^2-2*x2-5, x2^2-6, x3^2-3, x4^2-2, x2-x3*x4}
              Is zero?   T_TRUE

    [Sqrt(I) - (1+I)/Sqrt(2)]   =
    (2*x1-x2*x3-x2)/(2)  in  QQ(x1, x2, x3) where {x1 = Sqrt(1*x1+0  in  QQ(x1) where {x1 = Algebraic I}), x2 = Algebraic 1.4142136, x3 = Algebraic I} with ideal {x1^2-x3, x2^2-2, x3^2+1}
              Is zero?   T_TRUE

    [Exp(Pi*Sqrt(163)) - (640320^3 + 744)]   =
    (x1-262537412640768744)/(1)  in  QQ(x1) where {x1 = Exp((x1*x2)/(1)  in  QQ(x1, x2) where {x1 = Pi, x2 = Algebraic 12.767145} with ideal {x2^2-163})}
              Is zero?   T_FALSE

    cpu/wall(s): 0.014 0.013
    virt/peak/res/peak(MB): 36.16 36.16 8.93 8.93


machin.c
-------------------------------------------------------------------------------

This program checks several variations of Machin's formula

.. math ::

    \frac{\pi}{4} = 4 \operatorname{atan}\left(\frac{1}{5}\right) - \operatorname{atan}\left(\frac{1}{239}\right)

expressing `\pi` or logarithms of small integers in terms of
arctangents or hyperbolic arctangents of rational numbers.
The program actually evaluates 
`4 \operatorname{atan}\left(\tfrac{1}{5}\right) - \operatorname{atan}\left(\tfrac{1}{239}\right) - \tfrac{\pi}{4}`
(etc.) and prints the result, which should be precisely 0
(``0  in  QQ``), proving the identity.
Inverse trigonometric functions are not yet implemented in Calcium,
so the example program evaluates them using logarithms.

Sample output::

    > build/examples/machin 
    [(1)*Atan(1/1) - Pi/4]   =   0  in  QQ
    [(1)*Atan(1/2) + (1)*Atan(1/3) - Pi/4]   =   0  in  QQ
    [(2)*Atan(1/2) + (-1)*Atan(1/7) - Pi/4]   =   0  in  QQ
    [(2)*Atan(1/3) + (1)*Atan(1/7) - Pi/4]   =   0  in  QQ
    [(4)*Atan(1/5) + (-1)*Atan(1/239) - Pi/4]   =   0  in  QQ
    [(1)*Atan(1/2) + (1)*Atan(1/5) + (1)*Atan(1/8) - Pi/4]   =   0  in  QQ
    [(1)*Atan(1/3) + (1)*Atan(1/4) + (1)*Atan(1/7) + (1)*Atan(1/13) - Pi/4]   =   0  in  QQ
    [(12)*Atan(1/49) + (32)*Atan(1/57) + (-5)*Atan(1/239) + (12)*Atan(1/110443) - Pi/4]   =   0  in  QQ

    [(14)*Atanh(1/31) + (10)*Atanh(1/49) + (6)*Atanh(1/161) - Log(2)]   =   0  in  QQ
    [(22)*Atanh(1/31) + (16)*Atanh(1/49) + (10)*Atanh(1/161) - Log(3)]   =   0  in  QQ
    [(32)*Atanh(1/31) + (24)*Atanh(1/49) + (14)*Atanh(1/161) - Log(5)]   =   0  in  QQ
    [(144)*Atanh(1/251) + (54)*Atanh(1/449) + (-38)*Atanh(1/4801) + (62)*Atanh(1/8749) - Log(2)]   =   0  in  QQ
    [(228)*Atanh(1/251) + (86)*Atanh(1/449) + (-60)*Atanh(1/4801) + (98)*Atanh(1/8749) - Log(3)]   =   0  in  QQ
    [(334)*Atanh(1/251) + (126)*Atanh(1/449) + (-88)*Atanh(1/4801) + (144)*Atanh(1/8749) - Log(5)]   =   0  in  QQ
    [(404)*Atanh(1/251) + (152)*Atanh(1/449) + (-106)*Atanh(1/4801) + (174)*Atanh(1/8749) - Log(7)]   =   0  in  QQ

    cpu/wall(s): 0.03 0.029
    virt/peak/res/peak(MB): 35.57 35.57 8.80 8.80

sdpoly.c
-------------------------------------------------------------------------------

This program computes the coefficients of the Swinnerton-Dyer polynomial

.. math ::

    S_n = \prod (x \pm \sqrt{2} \pm \sqrt{3} \pm \sqrt{5} \pm \ldots \pm \sqrt{p_n})

where `p_n` denotes the `n`-th prime number and all combinations
of signs are taken. This polynomial has degree `2^n`.
The polynomial is expanded from its roots
using naive polynomial multiplication over :type:`ca_t` coefficients.
There are far more efficient ways to construct this polynomial;
this program simply illustrates that arithmetic in
multivariate number fields works smoothly.

The program prints the coefficients from `S_n`, from the constant
term to the coefficient of `x^{2^n}`.

Sample output::

    > build/examples/sdpoly 3
    576  in  QQ
    0  in  QQ
    -960  in  QQ
    0  in  QQ
    352  in  QQ
    0  in  QQ
    -40  in  QQ
    0  in  QQ
    1  in  QQ

    cpu/wall(s): 0.002 0.002
    virt/peak/res/peak(MB): 35.07 35.11 5.40 5.40

A big benchmark problem (output truncated)::

    > build/examples/sdpoly 10
    43567450015...212890625  in  QQ
    ...
    0  in  QQ
    1  in  QQ

    cpu/wall(s): 9.296 9.307
    virt/peak/res/peak(MB): 38.95 38.95 10.01 10.01

huge_expr.c
-------------------------------------------------------------------------------

This program proves equality of two complicated algebraic numbers.
More precisely, the program verifies
that `N = -(1 - |M|^2)^2` where *N* and *M* are given by huge symbolic
expressions involving nested square roots (about 7000
operations in total).

By default, the program runs the computation using :type:`qqbar_t` arithmetic.
This takes half a minute::

    > build/examples/huge_expr 
    Evaluating N...
    cpu/wall(s): 18.279 18.279
    Evaluating M...
    cpu/wall(s): 6.049 6.051
    Evaluating E = -(1-|M|^2)^2...
    cpu/wall(s): 0.595 0.595
    N ~ -0.16190853053311203695842869991458578203473645660641
    E ~ -0.16190853053311203695842869991458578203473645660641
    Testing E = N...
    cpu/wall(s): 0 0

    Equal = T_TRUE

    Total: cpu/wall(s): 24.927 24.93
    virt/peak/res/peak(MB): 56.61 68.64 28.73 40.70

To run the computation using :type:`ca_t` arithmetic instead, one
may pass the ``-ca`` flag. This currently takes much longer::

    > build/examples/huge_expr -ca
    Evaluating N...
    cpu/wall(s): 2.116 2.116
    Evaluating M...
    cpu/wall(s): 0.068 0.068
    Evaluating E = -(1-|M|^2)^2...
    cpu/wall(s): 0.043 0.043
    N ~ -0.16190853053311203695842869991458578203473645660641
    E ~ -0.16190853053311203695842869991458578203473645660641
    Testing E = N...
    cpu/wall(s): 176.235 176.242

    Equal = T_TRUE

    Total: cpu/wall(s): 178.465 178.472
    virt/peak/res/peak(MB): 55.92 67.88 29.80 41.76

This should be possible to improve significantly;
we keep this program as a benchmark for future optimizations
to the :type:`ca_t` type.

This simplification problem was posted in a help request for Sage
(https://ask.sagemath.org/question/52653/equality-of-algebraic-numbers-given-by-huge-symbolic-expressions/).
The C code has been generated from the symbolic expressions
using a Python script.



.. raw:: latex

    \newpage

