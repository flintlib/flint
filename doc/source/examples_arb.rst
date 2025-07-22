.. _examples-arb:

Arb example programs
===============================================================================

.. highlight:: text

See :ref:`examples` for general information about example programs.
Running::

    make examples

will compile the programs and place the binaries in
``build/examples``. The examples related to the Arb module are
documented below.

pi.c
-------------------------------------------------------------------------------

This program computes `\pi` to an accuracy of roughly *n* decimal digits
by calling the :func:`arb_const_pi` function with a
working precision of roughly `n \log_2(10)` bits.

Sample output, computing `\pi` to one million digits::

    > build/examples/pi 1000000
    precision = 3321933 bits... cpu/wall(s): 0.243 0.244
    virt/peak/res/peak(MB): 24.46 30.44 8.73 14.42
    [3.14159265358979323846{...999959 digits...}42209010610577945815 +/- 1.38e-1000000]

The program prints an interval guaranteed to contain `\pi`, and where
all displayed digits are correct up to an error of plus or minus
one unit in the last place (see :func:`arb_printn`).
By default, only the first and last few digits are printed.
Pass 0 as a second argument to print all digits (or pass *m* to
print *m* + 1 leading and *m* trailing digits, as above with
the default *m* = 20).

The program can optionally compute various other constants, and can
use multiple threads::

    > build/examples/pi 1000000 -threads 4
    precision = 3321933 bits... cpu/wall(s): 0.265 0.147
    virt/peak/res/peak(MB): 241.95 422.15 13.33 17.54
    [3.14159265358979323846{...999959 digits...}42209010610577945815 +/- 1.38e-1000000]
    > build/examples/pi 1000000 -constant e
    precision = 3321933 bits... cpu/wall(s): 0.09 0.09
    virt/peak/res/peak(MB): 25.56 29.19 9.58 13.11
    [2.71828182845904523536{...999959 digits...}01379817644769422819 +/- 1.39e-1000000]

pi_agm.c
-------------------------------------------------------------------------------

This program implements the Brent-Salamin AGM iteration to compute `\pi`.
This algorithm is not used by :func:`arb_const_pi` since the Chudnovsky
algorithm is faster in practice.

    > build/examples/pi_agm 100000
    precision = 332197 bits...
    Iteration 1 / 16: error bound 0.00391
    Iteration 2 / 16: error bound 2.98e-8
    Iteration 3 / 16: error bound 8.67e-19
    Iteration 4 / 16: error bound 1.84e-40
    Iteration 5 / 16: error bound 8.24e-84
    Iteration 6 / 16: error bound 8.28e-171
    Iteration 7 / 16: error bound 4.18e-345
    Iteration 8 / 16: error bound 5.34e-694
    Iteration 9 / 16: error bound 2.18e-1392
    Iteration 10 / 16: error bound 3.62e-2789
    Iteration 11 / 16: error bound 5.01e-5583
    Iteration 12 / 16: error bound 2.39e-11171
    Iteration 13 / 16: error bound 5.46e-22348
    Iteration 14 / 16: error bound 1.42e-44701
    Iteration 15 / 16: error bound 4.82e-89409
    Iteration 16 / 16: error bound 1.38e-178824
    cpu/wall(s): 0.028 0.028
    virt/peak/res/peak(MB): 19.19 19.19 7.17 7.20
    [3.14159265358979323846{...99959 digits...}76742080565549362465 +/- 3.91e-100000]

zeta_zeros.c
-------------------------------------------------------------------------------

This program computes one or several consecutive zeros of the
Riemann zeta function on the critical line::

    > build/examples/zeta_zeros -n 1 -count 10 -digits 30
    1	14.1347251417346937904572519836
    2	21.0220396387715549926284795939
    3	25.0108575801456887632137909926
    4	30.4248761258595132103118975306
    5	32.9350615877391896906623689641
    6	37.5861781588256712572177634807
    7	40.9187190121474951873981269146
    8	43.3270732809149995194961221654
    9	48.0051508811671597279424727494
    10	49.7738324776723021819167846786
    cpu/wall(s): 0.01 0.01
    virt/peak/res/peak(MB): 21.28 21.28 7.29 7.29

Five zeros starting with the millionth::

    > build/examples/zeta_zeros -n 1000000 -count 5 -digits 20
    1000000	600269.67701244495552
    1000001	600270.30109071169866
    1000002	600270.74787059436613
    1000003	600271.48637367364820
    1000004	600271.76148042593778
    cpu/wall(s): 0.03 0.03
    virt/peak/res/peak(MB): 21.41 21.41 7.41 7.41

The program supports the following options::

    zeta_zeros [-n n] [-count n] [-prec n] [-digits n] [-threads n] [-platt] [-noplatt] [-v] [-verbose] [-h] [-help]

With ``-platt``, Platt's algorithm is used, which may be faster when
computing many zeros of large index simultaneously.

bernoulli.c
-------------------------------------------------------------------------------

This program benchmarks computing the nth Bernoulli number exactly::

    > build/examples/bernoulli 1000000 -threads 8
    cpu/wall(s): 27.227 5.836
    virt/peak/res/peak(MB): 573.47 731.39 73.23 165.13

class_poly.c
-------------------------------------------------------------------------------

This program benchmarks computing Hilbert class polynomials::

    > build/examples/class_poly -1000004 -threads 8
    cpu/wall(s): 6.932 1.478
    virt/peak/res/peak(MB): 535.27 653.18 71.02 100.65
    degree = 624, bits = -37823

hilbert_matrix.c
-------------------------------------------------------------------------------

Given an input integer *n*, this program accurately computes the
determinant of the *n* by *n* Hilbert matrix.
Hilbert matrices are notoriously ill-conditioned: although the
entries are close to unit magnitude, the determinant `h_n`
decreases superexponentially (nearly as `1/4^{n^2}`) as
a function of *n*.
This program automatically doubles the working precision
until the ball computed for `h_n` by :func:`arb_mat_det`
does not contain zero.

Sample output::

    $ build/examples/hilbert_matrix 200
    prec=20: [+/- 1.32e-335]
    prec=40: [+/- 1.63e-545]
    prec=80: [+/- 1.30e-933]
    prec=160: [+/- 3.62e-1926]
    prec=320: [+/- 1.81e-4129]
    prec=640: [+/- 3.84e-8838]
    prec=1280: [2.955454297e-23924 +/- 8.29e-23935]
    success!
    cpu/wall(s): 8.494 8.513
    virt/peak/res/peak(MB): 134.98 134.98 111.57 111.57

Called with ``-eig n``, instead of computing the determinant,
the program computes the smallest eigenvalue of the Hilbert matrix
(in fact, it isolates all eigenvalues and prints the smallest eigenvalue)::

    $ build/examples/hilbert_matrix -eig 50
    prec=20: nan
    prec=40: nan
    prec=80: nan
    prec=160: nan
    prec=320: nan
    prec=640: [1.459157797e-74 +/- 2.49e-84]
    success!
    cpu/wall(s): 1.84 1.841
    virt/peak/res/peak(MB): 33.97 33.97 10.51 10.51

keiper_li.c
-------------------------------------------------------------------------------

Given an input integer *n*, this program rigorously computes numerical
values of the Keiper-Li coefficients
`\lambda_0, \ldots, \lambda_n`. The Keiper-Li coefficients
have the property that `\lambda_n > 0` for all `n > 0` if and only if the
Riemann hypothesis is true. This program was used for the record
computations described in [Joh2013]_ (the paper describes
the algorithm in some more detail).

The program takes the following parameters::

    keiper_li n [-prec prec] [-threads num_threads] [-out out_file]

The program prints the first and last few coefficients. It can optionally
write all the computed data to a file. The working precision defaults
to a value that should give all the coefficients to a few digits of
accuracy, but can optionally be set higher (or lower).
On a multicore system, using several threads results in faster
execution.

Sample output::

    > build/examples/keiper_li 1000 -threads 2
    zeta: cpu/wall(s): 0.4 0.244
    virt/peak/res/peak(MB): 167.98 294.69 5.09 7.43
    log: cpu/wall(s): 0.03 0.038
    gamma: cpu/wall(s): 0.02 0.016
    binomial transform: cpu/wall(s): 0.01 0.018
    0: -0.69314718055994530941723212145817656807550013436026 +/- 6.5389e-347
    1: 0.023095708966121033814310247906495291621932127152051 +/- 2.0924e-345
    2: 0.046172867614023335192864243096033943387066108314123 +/- 1.674e-344
    3: 0.0692129735181082679304973488726010689942120263932 +/- 5.0219e-344
    4: 0.092197619873060409647627872409439018065541673490213 +/- 2.0089e-343
    5: 0.11510854289223549048622128109857276671349132303596 +/- 1.0044e-342
    6: 0.13792766871372988290416713700341666356138966078654 +/- 6.0264e-342
    7: 0.16063715965299421294040287257385366292282442046163 +/- 2.1092e-341
    8: 0.18321945964338257908193931774721859848998098273432 +/- 8.4368e-341
    9: 0.20565733870917046170289387421343304741236553410044 +/- 7.5931e-340
    10: 0.22793393631931577436930340573684453380748385942738 +/- 7.5931e-339
    991: 2.3196617961613367928373899656994682562101430813341 +/- 2.461e-11
    992: 2.3203766239254884035349896518332550233162909717288 +/- 9.5363e-11
    993: 2.321092061239733282811659116333262802034375592414 +/- 1.8495e-10
    994: 2.3218073540188462110258826121503870112747188888893 +/- 3.5907e-10
    995: 2.3225217392815185726928702951225314023773358152533 +/- 6.978e-10
    996: 2.3232344485814623873333223609413703912358283071281 +/- 1.3574e-09
    997: 2.3239447114886014522889542667580382034526509232475 +/- 2.6433e-09
    998: 2.3246517591032700808344143240352605148856869322209 +/- 5.1524e-09
    999: 2.3253548275861382119812576052060526988544993162101 +/- 1.0053e-08
    1000: 2.3260531616864664574065046940832238158044982041872 +/- 3.927e-08
    virt/peak/res/peak(MB): 170.18 294.69 7.51 7.51

logistic.c
-------------------------------------------------------------------------------

This program computes the *n*-th iterate of the logistic map defined
by `x_{n+1} = r x_n (1 - x_n)` where `r` and `x_0` are given.
It takes the following parameters::

    logistic n [x_0] [r] [digits]

The inputs `x_0`, *r* and *digits* default to 0.5, 3.75 and 10 respectively.
The computation is automatically restarted with doubled precision
until the result is accurate to *digits* decimal digits.

Sample output::

    > build/examples/logistic 10
    Trying prec=64 bits...success!
    cpu/wall(s): 0 0.001
    x_10 = [0.6453672908 +/- 3.10e-11]

    > build/examples/logistic 100
    Trying prec=64 bits...ran out of accuracy at step 18
    Trying prec=128 bits...ran out of accuracy at step 53
    Trying prec=256 bits...success!
    cpu/wall(s): 0 0
    x_100 = [0.8882939923 +/- 1.60e-11]

    > build/examples/logistic 10000
    Trying prec=64 bits...ran out of accuracy at step 18
    Trying prec=128 bits...ran out of accuracy at step 53
    Trying prec=256 bits...ran out of accuracy at step 121
    Trying prec=512 bits...ran out of accuracy at step 256
    Trying prec=1024 bits...ran out of accuracy at step 525
    Trying prec=2048 bits...ran out of accuracy at step 1063
    Trying prec=4096 bits...ran out of accuracy at step 2139
    Trying prec=8192 bits...ran out of accuracy at step 4288
    Trying prec=16384 bits...ran out of accuracy at step 8584
    Trying prec=32768 bits...success!
    cpu/wall(s): 0.859 0.858
    x_10000 = [0.8242048008 +/- 4.35e-11]

    > build/examples/logistic 1234 0.1 3.99 30
    Trying prec=64 bits...ran out of accuracy at step 0
    Trying prec=128 bits...ran out of accuracy at step 10
    Trying prec=256 bits...ran out of accuracy at step 76
    Trying prec=512 bits...ran out of accuracy at step 205
    Trying prec=1024 bits...ran out of accuracy at step 461
    Trying prec=2048 bits...ran out of accuracy at step 974
    Trying prec=4096 bits...success!
    cpu/wall(s): 0.009 0.009
    x_1234 = [0.256445391958651410579677945635 +/- 3.92e-31]

real_roots.c
-------------------------------------------------------------------------------

This program isolates the roots of a function on the interval `(a,b)`
(where *a* and *b* are input as double-precision literals)
using the routines in the :ref:`arb_calc <arb-calc>` module.
The program takes the following arguments::

    real_roots function a b [-refine d] [-verbose] [-maxdepth n] [-maxeval n] [-maxfound n] [-prec n]

The following functions (specified by an integer code) are implemented:

  * 0 - `Z(x)` (Riemann-Siegel Z-function)
  * 1 - `\sin(x)`
  * 2 - `\sin(x^2)`
  * 3 - `\sin(1/x)`
  * 4 - `\operatorname{Ai}(x)` (Airy function)
  * 5 - `\operatorname{Ai}'(x)` (Airy function)
  * 6 - `\operatorname{Bi}(x)` (Airy function)
  * 7 - `\operatorname{Bi}'(x)` (Airy function)

The following options are available:

  * ``-refine d``: If provided, after isolating the roots, attempt to refine
    the roots to *d* digits of accuracy using a few bisection steps followed
    by Newton's method with adaptive precision, and then print them.

  * ``-verbose``: Print more information.

  * ``-maxdepth n``: Stop searching after *n* recursive subdivisions.

  * ``-maxeval n``: Stop searching after approximately *n* function evaluations
    (the actual number evaluations will be a small multiple of this).

  * ``-maxfound n``: Stop searching after having found *n* isolated roots.

  * ``-prec n``: Working precision to use for the root isolation.

With *function* 0, the program isolates roots of the Riemann zeta function
on the critical line, and guarantees that no roots are missed
(see `zeta_zeros.c` for a far more efficient way to do this)::

    > build/examples/real_roots 0 0.0 50.0 -verbose
    interval: [0, 50]
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    found isolated root in: [14.111328125, 14.16015625]
    found isolated root in: [20.99609375, 21.044921875]
    found isolated root in: [25, 25.048828125]
    found isolated root in: [30.419921875, 30.4443359375]
    found isolated root in: [32.91015625, 32.958984375]
    found isolated root in: [37.548828125, 37.59765625]
    found isolated root in: [40.91796875, 40.966796875]
    found isolated root in: [43.310546875, 43.3349609375]
    found isolated root in: [47.998046875, 48.0224609375]
    found isolated root in: [49.755859375, 49.7802734375]
    ---------------------------------------------------------------
    Found roots: 10
    Subintervals possibly containing undetected roots: 0
    Function evaluations: 3058
    cpu/wall(s): 0.202 0.202
    virt/peak/res/peak(MB): 26.12 26.14 2.76 2.76

Find just one root and refine it to approximately 75 digits::

    > build/examples/real_roots 0 0.0 50.0 -maxfound 1 -refine 75
    interval: [0, 50]
    maxdepth = 30, maxeval = 100000, maxfound = 1, low_prec = 30
    refined root (0/8):
    [14.134725141734693790457251983562470270784257115699243175685567460149963429809 +/- 2.57e-76]

    ---------------------------------------------------------------
    Found roots: 1
    Subintervals possibly containing undetected roots: 7
    Function evaluations: 761
    cpu/wall(s): 0.055 0.056
    virt/peak/res/peak(MB): 26.12 26.14 2.75 2.75

Find the first few roots of an Airy function and refine them to 50 digits each::

    > build/examples/real_roots 4 -10 0 -refine 50
    interval: [-10, 0]
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    refined root (0/6):
    [-9.022650853340980380158190839880089256524677535156083 +/- 4.85e-52]

    refined root (1/6):
    [-7.944133587120853123138280555798268532140674396972215 +/- 1.92e-52]

    refined root (2/6):
    [-6.786708090071758998780246384496176966053882477393494 +/- 3.84e-52]

    refined root (3/6):
    [-5.520559828095551059129855512931293573797214280617525 +/- 1.05e-52]

    refined root (4/6):
    [-4.087949444130970616636988701457391060224764699108530 +/- 2.46e-52]

    refined root (5/6):
    [-2.338107410459767038489197252446735440638540145672388 +/- 1.48e-52]

    ---------------------------------------------------------------
    Found roots: 6
    Subintervals possibly containing undetected roots: 0
    Function evaluations: 200
    cpu/wall(s): 0.003 0.003
    virt/peak/res/peak(MB): 26.12 26.14 2.24 2.24

Find roots of `\sin(x^2)` on `(0,100)`. The algorithm cannot isolate
the root at `x = 0` (it is at the endpoint of the interval, and in any
case a root of multiplicity higher than one). The failure is reported::

    > build/examples/real_roots 2 0 100
    interval: [0, 100]
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    ---------------------------------------------------------------
    Found roots: 3183
    Subintervals possibly containing undetected roots: 1
    Function evaluations: 34058
    cpu/wall(s): 0.032 0.032
    virt/peak/res/peak(MB): 26.32 26.37 2.04 2.04

This does not miss any roots::

    > build/examples/real_roots 2 1 100
    interval: [1, 100]
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    ---------------------------------------------------------------
    Found roots: 3183
    Subintervals possibly containing undetected roots: 0
    Function evaluations: 34039
    cpu/wall(s): 0.023 0.023
    virt/peak/res/peak(MB): 26.32 26.37 2.01 2.01

Looking for roots of `\sin(1/x)` on `(0,1)`, the algorithm finds many roots,
but will never find all of them since there are infinitely many::

    > build/examples/real_roots 3 0.0 1.0
    interval: [0, 1]
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    ---------------------------------------------------------------
    Found roots: 10198
    Subintervals possibly containing undetected roots: 24695
    Function evaluations: 202587
    cpu/wall(s): 0.171 0.171
    virt/peak/res/peak(MB): 28.39 30.38 4.05 4.05

Remark: the program always computes rigorous containing intervals
for the roots, but the accuracy after refinement could be less than *d* digits.

poly_roots.c
-------------------------------------------------------------------------------

This program finds the complex roots of an integer polynomial
by calling :func:`arb_fmpz_poly_complex_roots`, which in turn calls
:func:`acb_poly_find_roots` with increasing
precision until the roots certainly have been isolated.
The program takes the following arguments::

    poly_roots [-refine d] [-print d] <poly>

    Isolates all the complex roots of a polynomial with integer coefficients.

    If -refine d is passed, the roots are refined to a relative tolerance
    better than 10^(-d). By default, the roots are only computed to sufficient
    accuracy to isolate them. The refinement is not currently done efficiently.

    If -print d is passed, the computed roots are printed to d decimals.
    By default, the roots are not printed.

    The polynomial can be specified by passing the following as <poly>:

    a <n>          Easy polynomial 1 + 2x + ... + (n+1)x^n
    t <n>          Chebyshev polynomial T_n
    u <n>          Chebyshev polynomial U_n
    p <n>          Legendre polynomial P_n
    c <n>          Cyclotomic polynomial Phi_n
    s <n>          Swinnerton-Dyer polynomial S_n
    b <n>          Bernoulli polynomial B_n
    w <n>          Wilkinson polynomial W_n
    e <n>          Taylor series of exp(x) truncated to degree n
    m <n> <m>      The Mignotte-like polynomial x^n + (100x+1)^m, n > m
    coeffs <c0 c1 ... cn>        c0 + c1 x + ... + cn x^n

    Concatenate to multiply polynomials, e.g.: p 5 t 6 coeffs 1 2 3
    for P_5(x)*T_6(x)*(1+2x+3x^2)

This finds the roots of the Wilkinson polynomial with roots at the
positive integers 1, 2, ..., 100::

    > build/examples/poly_roots -print 15 w 100
    computing squarefree factorization...
    cpu/wall(s): 0.001 0.001
    roots with multiplicity 1
    searching for 100 roots, 100 deflated
    prec=32: 0 isolated roots | cpu/wall(s): 0.098 0.098
    prec=64: 0 isolated roots | cpu/wall(s): 0.247 0.247
    prec=128: 0 isolated roots | cpu/wall(s): 0.498 0.497
    prec=256: 0 isolated roots | cpu/wall(s): 0.713 0.713
    prec=512: 100 isolated roots | cpu/wall(s): 0.104 0.105
    done!
    [1.00000000000000 +/- 3e-20]
    [2.00000000000000 +/- 3e-19]
    [3.00000000000000 +/- 1e-19]
    [4.00000000000000 +/- 1e-19]
    [5.00000000000000 +/- 1e-19]
    ...
    [96.0000000000000 +/- 1e-17]
    [97.0000000000000 +/- 1e-17]
    [98.0000000000000 +/- 3e-17]
    [99.0000000000000 +/- 3e-17]
    [100.000000000000 +/- 3e-17]
    cpu/wall(s): 1.664 1.664

This finds the roots of a Bernoulli polynomial which has both real
and complex roots::

    > build/examples/poly_roots -refine 100 -print 20 b 16
    computing squarefree factorization...
    cpu/wall(s): 0.001 0
    roots with multiplicity 1
    searching for 16 roots, 16 deflated
    prec=32: 16 isolated roots | cpu/wall(s): 0.006 0.006
    prec=64: 16 isolated roots | cpu/wall(s): 0.001 0.001
    prec=128: 16 isolated roots | cpu/wall(s): 0.001 0.001
    prec=256: 16 isolated roots | cpu/wall(s): 0.001 0.002
    prec=512: 16 isolated roots | cpu/wall(s): 0.002 0.001
    done!
    [-0.94308706466055783383 +/- 2.02e-21]
    [-0.75534059252067985752 +/- 2.70e-21]
    [-0.24999757119077421009 +/- 4.27e-21]
    [0.24999757152512726002 +/- 4.43e-21]
    [0.75000242847487273998 +/- 4.43e-21]
    [1.2499975711907742101 +/- 1.43e-20]
    [1.7553405925206798575 +/- 1.74e-20]
    [1.9430870646605578338 +/- 3.21e-20]
    [-0.99509334829256233279 +/- 9.42e-22] + [0.44547958157103608805 +/- 3.59e-21]*I
    [-0.99509334829256233279 +/- 9.42e-22] + [-0.44547958157103608805 +/- 3.59e-21]*I
    [1.9950933482925623328 +/- 1.10e-20] + [0.44547958157103608805 +/- 3.59e-21]*I
    [1.9950933482925623328 +/- 1.10e-20] + [-0.44547958157103608805 +/- 3.59e-21]*I
    [-0.92177327714429290564 +/- 4.68e-21] + [-1.0954360955079385542 +/- 1.71e-21]*I
    [-0.92177327714429290564 +/- 4.68e-21] + [1.0954360955079385542 +/- 1.71e-21]*I
    [1.9217732771442929056 +/- 3.54e-20] + [1.0954360955079385542 +/- 1.71e-21]*I
    [1.9217732771442929056 +/- 3.54e-20] + [-1.0954360955079385542 +/- 1.71e-21]*I
    cpu/wall(s): 0.011 0.012

Roots are automatically separated by multiplicity by performing an initial
squarefree factorization::

    > build/examples/poly_roots -print 5 p 5 p 5 t 7 coeffs 1 5 10 10 5 1
    computing squarefree factorization...
    cpu/wall(s): 0 0
    roots with multiplicity 1
    searching for 6 roots, 3 deflated
    prec=32: 3 isolated roots | cpu/wall(s): 0 0.001
    done!
    [-0.97493 +/- 2.10e-6]
    [-0.78183 +/- 1.49e-6]
    [-0.43388 +/- 3.75e-6]
    [0.43388 +/- 3.75e-6]
    [0.78183 +/- 1.49e-6]
    [0.97493 +/- 2.10e-6]
    roots with multiplicity 2
    searching for 4 roots, 2 deflated
    prec=32: 2 isolated roots | cpu/wall(s): 0 0
    done!
    [-0.90618 +/- 1.56e-7]
    [-0.53847 +/- 6.91e-7]
    [0.53847 +/- 6.91e-7]
    [0.90618 +/- 1.56e-7]
    roots with multiplicity 3
    searching for 1 roots, 0 deflated
    prec=32: 0 isolated roots | cpu/wall(s): 0 0
    done!
    0
    roots with multiplicity 5
    searching for 1 roots, 1 deflated
    prec=32: 1 isolated roots | cpu/wall(s): 0 0
    done!
    -1.0000
    cpu/wall(s): 0 0.001

zeta_zeros.c
-------------------------------------------------------------------------------

This program finds the imaginary parts of consecutive nontrivial zeros
of the Riemann zeta function by calling either
:func:`acb_dirichlet_hardy_z_zeros` or
:func:`acb_dirichlet_platt_local_hardy_z_zeros` depending on the height
of the zeros and the number of zeros requested.
The program takes the following arguments::

    zeta_zeros [-n n] [-count n] [-prec n] [-threads n] [-platt] [-noplatt] [-v] [-verbose] [-h] [-help]

    > build/examples/zeta_zeros -n 1048449114 -count 2
    1048449114      [388858886.0022851217767970582 +/- 7.46e-20]
    1048449115      [388858886.0023936897027167201 +/- 7.59e-20]
    cpu/wall(s): 0.255 0.255
    virt/peak/res/peak(MB): 26.77 26.77 7.88 7.88

complex_plot.c
-------------------------------------------------------------------------------

This program plots one of the predefined functions over a complex
interval `[x_a, x_b] + [y_a, y_b]i` using domain coloring, at
a resolution of *xn* times *yn* pixels.

The program takes the parameters::

    complex_plot [-range xa xb ya yb] [-size xn yn] [-color n] [-threads n] <func>

Defaults parameters are `[-10,10] + [-10,10]i` and *xn* = *yn* = 512.

A color function can be selected with -color. Valid options
are 0 (phase=hue, magnitude=brightness) and 1 (phase only,
white-gold-black-blue-white counterclockwise).

The output is written to ``arbplot.ppm``. If you have ImageMagick,
run ``convert arbplot.ppm arbplot.png`` to get a PNG.

Function codes ``<func>`` are:

  * ``gamma``   - Gamma function
  * ``digamma`` - Digamma function
  * ``lgamma``  - Logarithmic gamma function
  * ``zeta``    - Riemann zeta function
  * ``erf``     - Error function
  * ``ai``      - Airy function Ai
  * ``bi``      - Airy function Bi
  * ``besselj`` - Bessel function `J_0`
  * ``bessely`` - Bessel function `Y_0`
  * ``besseli`` - Bessel function `I_0`
  * ``besselk`` - Bessel function `K_0`
  * ``modj``    - Modular j-function
  * ``modeta``  - Dedekind eta function
  * ``barnesg`` - Barnes G-function
  * ``agm``     - Arithmetic geometric mean

The function is just sampled at point values; no attempt is made to resolve
small features by adaptive subsampling.

For example, the following plots the Riemann zeta function around
a portion of the critical strip with imaginary part between 100 and 140::

    > build/examples/complex_plot zeta -range -10 10 100 140 -size 256 512

For parallel computation on a multicore system, use ``-threads n``.

lvalue.c
-------------------------------------------------------------------------------

This program evaluates Dirichlet L-functions. It takes the following input::

    > build/examples/lvalue
    lvalue [-character q n] [-re a] [-im b] [-prec p] [-z] [-deflate] [-len l]

    Print value of Dirichlet L-function at s = a+bi.
    Default a = 0.5, b = 0, p = 53, (q, n) = (1, 0) (Riemann zeta)
    [-z]       - compute Z(s) instead of L(s)
    [-deflate] - remove singular term at s = 1
    [-len l]   - compute l terms in Taylor series at s

Evaluating the Riemann zeta function and
the Dirichlet beta function at `s = 2`::

    > build/examples/lvalue -re 2 -prec 128
    L(s) = [1.64493406684822643647241516664602518922 +/- 4.37e-39]
    cpu/wall(s): 0.001 0.001
    virt/peak/res/peak(MB): 26.86 26.88 2.05 2.05

    > build/examples/lvalue -character 4 3 -re 2 -prec 128
    L(s) = [0.91596559417721901505460351493238411077 +/- 7.86e-39]
    cpu/wall(s): 0.002 0.003
    virt/peak/res/peak(MB): 26.86 26.88 2.31 2.31

Evaluating the L-function for character number 101 modulo 1009
at `s = 1/2` and `s = 1`::

    > build/examples/lvalue -character 1009 101
    L(s) = [-0.459256562383872 +/- 5.24e-16] + [1.346937111206009 +/- 3.03e-16]*I
    cpu/wall(s): 0.012 0.012
    virt/peak/res/peak(MB): 26.86 26.88 2.30 2.30

    > build/examples/lvalue -character 1009 101 -re 1
    L(s) = [0.657952586112728 +/- 6.02e-16] + [1.004145273214022 +/- 3.10e-16]*I
    cpu/wall(s): 0.017 0.018
    virt/peak/res/peak(MB): 26.86 26.88 2.30 2.30

Computing the first few coefficients in the Laurent series of the
Riemann zeta function at `s = 1`::

    > build/examples/lvalue -re 1 -deflate -len 8
    L(s) = [0.577215664901532861 +/- 5.29e-19]
    L'(s) = [0.072815845483676725 +/- 2.68e-19]
    [x^2] L(s+x) = [-0.004845181596436159 +/- 3.87e-19]
    [x^3] L(s+x) = [-0.000342305736717224 +/- 4.20e-19]
    [x^4] L(s+x) = [9.6890419394471e-5 +/- 2.40e-19]
    [x^5] L(s+x) = [-6.6110318108422e-6 +/- 4.51e-20]
    [x^6] L(s+x) = [-3.316240908753e-7 +/- 3.85e-20]
    [x^7] L(s+x) = [1.0462094584479e-7 +/- 7.78e-21]
    cpu/wall(s): 0.003 0.004
    virt/peak/res/peak(MB): 26.86 26.88 2.30 2.30

Evaluating the Riemann zeta function near the first nontrivial root::

    > build/examples/lvalue -re 0.5 -im 14.134725
    L(s) = [1.76743e-8 +/- 1.93e-14] + [-1.110203e-7 +/- 2.84e-14]*I
    cpu/wall(s): 0.001 0.001
    virt/peak/res/peak(MB): 26.86 26.88 2.31 2.31

    > build/examples/lvalue -z -re 14.134725 -prec 200
    Z(s) = [-1.12418349839417533300111494358128257497862927935658e-7 +/- 4.62e-58]
    cpu/wall(s): 0.001 0.001
    virt/peak/res/peak(MB): 26.86 26.88 2.57 2.57

    > build/examples/lvalue -z -re 14.134725 -len 4
    Z(s) = [-1.124184e-7 +/- 7.00e-14]
    Z'(s) = [0.793160414884 +/- 4.09e-13]
    [x^2] Z(s+x) = [0.065164586492 +/- 5.39e-13]
    [x^3] Z(s+x) = [-0.020707762705 +/- 5.37e-13]
    cpu/wall(s): 0.002 0.003
    virt/peak/res/peak(MB): 26.86 26.88 2.57 2.57

lcentral.c
-------------------------------------------------------------------------------

This program computes the central value `L(1/2)` for each Dirichlet L-function
character modulo *q* for each *q* in the range *qmin* to *qmax*. Usage::

    > build/examples/lcentral
    Computes central values (s = 0.5) of Dirichlet L-functions.

    usage: build/examples/lcentral [--quiet] [--check] [--prec <bits>] qmin qmax

The first few values::

    > build/examples/lcentral 1 8
    3,2: [0.48086755769682862618122006324 +/- 7.35e-30]
    4,3: [0.66769145718960917665869092930 +/- 1.62e-30]
    5,2: [0.76374788011728687822451215264 +/- 2.32e-30] + [0.21696476751886069363858659310 +/- 3.06e-30]*I
    5,4: [0.23175094750401575588338366176 +/- 2.21e-30]
    5,3: [0.76374788011728687822451215264 +/- 2.32e-30] + [-0.21696476751886069363858659310 +/- 3.06e-30]*I
    7,3: [0.71394334376831949285993820742 +/- 1.21e-30] + [0.47490218277139938263745243935 +/- 4.52e-30]*I
    7,2: [0.31008936259836766059195052534 +/- 5.29e-30] + [-0.07264193137017790524562171245 +/- 5.48e-30]*I
    7,6: [1.14658566690370833367712697646 +/- 1.95e-30]
    7,4: [0.31008936259836766059195052534 +/- 5.29e-30] + [0.07264193137017790524562171245 +/- 5.48e-30]*I
    7,5: [0.71394334376831949285993820742 +/- 1.21e-30] + [-0.47490218277139938263745243935 +/- 4.52e-30]*I
    8,5: [0.37369171291254730738158695002 +/- 4.01e-30]
    8,3: [1.10042140952554837756713576997 +/- 3.37e-30]
    cpu/wall(s): 0.002 0.003
    virt/peak/res/peak(MB): 26.32 26.34 2.35 2.35

Testing a large *q*::

    > build/examples/lcentral --quiet --check --prec 256 100000 100000
    cpu/wall(s): 1.668 1.667
    virt/peak/res/peak(MB): 35.67 46.66 11.67 22.61

It is conjectured that the central value never vanishes. Running with ``--check``
verifies that the interval certainly is nonzero. This can fail with
insufficient precision::

    > build/examples/lcentral --check --prec 15 100000 100000
    100000,71877: [0.1 +/- 0.0772] + [+/- 0.136]*I
    100000,90629: [2e+0 +/- 0.106] + [+/- 0.920]*I
    100000,28133: [+/- 0.811] + [-2e+0 +/- 0.501]*I
    100000,3141: [0.8 +/- 0.0407] + [-0.1 +/- 0.0243]*I
    100000,53189: [4.0 +/- 0.0826] + [+/- 0.107]*I
    100000,53253: [1.9 +/- 0.0855] + [-3.9 +/- 0.0681]*I
    Value could be zero!
    100000,53381: [+/- 0.0329] + [+/- 0.0413]*I
    Aborted

integrals.c
-------------------------------------------------------------------------------

This program computes integrals using :func:`acb_calc_integrate`.
Invoking the program without parameters shows usage::

    > build/examples/integrals
    Compute integrals using acb_calc_integrate.
    Usage: integrals -i n [-prec p] [-tol eps] [-twice] [...]

    -i n       - compute integral n (0 <= n <= 23), or "-i all"
    -prec p    - precision in bits (default p = 64)
    -goal p    - approximate relative accuracy goal (default p)
    -tol eps   - approximate absolute error goal (default 2^-p)
    -twice     - run twice (to see overhead of computing nodes)
    -heap      - use heap for subinterval queue
    -verbose   - show information
    -verbose2  - show more information
    -deg n     - use quadrature degree up to n
    -eval n    - limit number of function evaluations to n
    -depth n   - limit subinterval queue size to n
    -threads n - use parallel computation with n threads

    Implemented integrals:
    I0 = int_0^100 sin(x) dx
    I1 = 4 int_0^1 1/(1+x^2) dx
    I2 = 2 int_0^{inf} 1/(1+x^2) dx   (using domain truncation)
    I3 = 4 int_0^1 sqrt(1-x^2) dx
    I4 = int_0^8 sin(x+exp(x)) dx
    I5 = int_1^101 floor(x) dx
    I6 = int_0^1 |x^4+10x^3+19x^2-6x-6| exp(x) dx
    I7 = 1/(2 pi i) int zeta(s) ds  (closed path around s = 1)
    I8 = int_0^1 sin(1/x) dx  (slow convergence, use -heap and/or -tol)
    I9 = int_0^1 x sin(1/x) dx  (slow convergence, use -heap and/or -tol)
    I10 = int_0^10000 x^1000 exp(-x) dx
    I11 = int_1^{1+1000i} gamma(x) dx
    I12 = int_{-10}^{10} sin(x) + exp(-200-x^2) dx
    I13 = int_{-1020}^{-1010} exp(x) dx  (use -tol 0 for relative error)
    I14 = int_0^{inf} exp(-x^2) dx   (using domain truncation)
    I15 = int_0^1 sech(10(x-0.2))^2 + sech(100(x-0.4))^4 + sech(1000(x-0.6))^6 dx
    I16 = int_0^8 (exp(x)-floor(exp(x))) sin(x+exp(x)) dx  (use higher -eval)
    I17 = int_0^{inf} sech(x) dx   (using domain truncation)
    I18 = int_0^{inf} sech^3(x) dx   (using domain truncation)
    I19 = int_0^1 -log(x)/(1+x) dx   (using domain truncation)
    I20 = int_0^{inf} x exp(-x)/(1+exp(-x)) dx   (using domain truncation)
    I21 = int_C wp(x)/x^(11) dx   (contour for 10th Laurent coefficient of Weierstrass p-function)
    I22 = N(1000) = count zeros with 0 < t <= 1000 of zeta(s) using argument principle
    I23 = int_0^{1000} W_0(x) dx
    I24 = int_0^pi max(sin(x), cos(x)) dx
    I25 = int_{-1}^1 erf(x/sqrt(0.0002)*0.5+1.5)*exp(-x) dx
    I26 = int_{-10}^10 Ai(x) dx
    I27 = int_0^10 (x-floor(x)-1/2) max(sin(x),cos(x)) dx
    I28 = int_{-1-i}^{-1+i} sqrt(x) dx
    I29 = int_0^{inf} exp(-x^2+ix) dx   (using domain truncation)
    I30 = int_0^{inf} exp(-x) Ai(-x) dx   (using domain truncation)
    I31 = int_0^pi x sin(x) / (1 + cos(x)^2) dx

A few examples::

    build/examples/integrals -i 4
    I4 = int_0^8 sin(x+exp(x)) dx ...
    cpu/wall(s): 0.02 0.02
    I4 = [0.34740017265725 +/- 3.95e-15]

    > build/examples/integrals -i 3 -prec 333 -tol 1e-80
    I3 = 4 int_0^1 sqrt(1-x^2) dx ...
    cpu/wall(s): 0.024 0.024
    I3 = [3.141592653589793238462643383279502884197169399375105820974944592307816406286209 +/- 4.24e-79]

    > build/examples/integrals -i 9 -heap
    I9 = int_0^1 x sin(1/x) dx  (slow convergence, use -heap and/or -tol) ...
    cpu/wall(s): 0.019 0.018
    I9 = [0.3785300 +/- 3.17e-8]

fpwrap.c
-------------------------------------------------------------------------------

This program demonstrates calling the floating-point wrapper::

    > build/examples/fpwrap
    zeta(2) = 1.644934066848226
    zeta(0.5 + 123i) = 0.006252861175594465 + 0.08206030514520983i

functions_benchmark.c
-------------------------------------------------------------------------------

This program benchmarks performance of some standard functions.


.. highlight:: c

