.. _using:

Using ball arithmetic
===============================================================================

This section gives an introduction to working with
real numbers in Arb (see :ref:`arb` for the API and technical documentation).
The general principles carry over to complex numbers, polynomials and
matrices.

Ball semantics
-------------------------------------------------------------------------------

Let `f : A \to B` be a function.
A ball implementation of `f` is a function `F` that maps sets `X \subseteq A`
to sets `F(X) \subseteq B` subject to the following rule:

    For all `x \in X`,
    we have `f(x) \in F(X)`.

In other words, `F(X)` is an *enclosure* for the set `\{f(x) : x \in X\}`.
This rule is sometimes called the *inclusion principle*.

Throughout the documentation (except where otherwise noted),
we will simply write `f(x)` instead of `F(X)`
when describing ball implementations of pointwise-defined mathematical
functions, understanding that the input is a set of point values and that
the output is an enclosure.

General subsets of `\mathbb{R}` are not possible to
represent on a computer. Instead, we work with subsets of the form
`[m \pm r] = [m-r, m+r]` where the midpoint *m* and radius *r* are binary
floating-point numbers, i.e. numbers of the form `u 2^v` with `u, v \in \mathbb{Z}`
(to make this scheme complete,
we also need to adjoin the special floating-point
values `-\infty`, `+\infty` and `\operatorname{NaN}`).

Given a ball `[m \pm r]` with `m \in \mathbb{R}` (not necessarily a
floating-point number),
we can always round *m* to a nearby floating-point number that has at most
most *prec* bits in the component *u*,
and add an upper bound for the rounding error to *r*.
In Arb, ball functions that take a *prec* argument as input
(e.g. :func:`arb_add`) always round their output to *prec* bits.
Some functions are always exact (e.g. :func:`arb_neg`), and thus do not take a *prec* argument.

The programming interface resembles that of GMP.
Each :type:`arb_t` variable must be initialized with :func:`arb_init` before use
(this also sets its value to zero), and deallocated with :func:`arb_clear`
after use. Variables have pass-by-reference semantics.
In the list of arguments to a function, output variables come first,
followed by input variables, and finally the precision:

.. code-block:: c

    #include "arb.h"

    int main()
    {
        arb_t x, y;
        arb_init(x); arb_init(y);
        arb_set_ui(x, 3);       /* x = 3 */
        arb_const_pi(y, 128);   /* y = pi, to 128 bits */
        arb_sub(y, y, x, 53);   /* y = y - x, to 53 bits */
        arb_clear(x); arb_clear(y);
    }

Binary and decimal
-------------------------------------------------------------------------------

While the internal representation uses binary floating-point numbers,
it is usually preferable to print numbers in decimal. The binary-to-decimal
conversion generally requires rounding. Three different methods
are available for printing a number to standard output:

* :func:`arb_print` shows the exact internal representation of a ball, with binary exponents.
* :func:`arb_printd` shows an inexact view of the internal representation, approximated by decimal floating-point numbers.
* :func:`arb_printn` shows a *decimal ball* that is guaranteed to be an enclosure of the binary
  floating-point ball. By default, it only prints digits in the midpoint that are certain to
  be correct, up to an error of at most one unit in the last place.
  Converting from binary to decimal is generally inexact, and the output of this
  method takes this rounding into account when printing the radius.

This snippet computes a 53-bit enclosure of `\pi` and prints it
in three ways:

.. code-block:: c

    arb_const_pi(x, 53);
    arb_print(x); printf("\n");
    arb_printd(x, 20); printf("\n");
    arb_printn(x, 20, 0); printf("\n");

The output is:

.. code-block:: text

    (884279719003555 * 2^-48) +/- (536870913 * 2^-80)
    3.141592653589793116 +/- 4.4409e-16
    [3.141592653589793 +/- 5.61e-16]

The :func:`arb_get_str` and :func:`arb_set_str` methods are useful for
converting rigorously between decimal strings and binary balls
(:func:`arb_get_str` produces the same string as :func:`arb_printn`,
and :func:`arb_set_str` can parse such strings back).

A potential mistake is to create a ball from a ``double`` constant
such as ``2.3``, when this actually represents
``2.29999999999999982236431605997495353221893310546875``.
To produce a ball containing the rational number
`23/10`, one of the following can be used:

.. code-block:: c

    arb_set_str(x, "2.3", prec)

    arb_set_ui(x, 23);
    arb_div_ui(x, x, 10, prec)

    fmpq_set_si(q, 23, 10);   /* q is a FLINT fmpq_t */
    arb_set_fmpq(x, q, prec);

Quality of enclosures
-------------------------------------------------------------------------------

The main problem when working with ball arithmetic (or interval arithmetic)
is *overestimation*. In general, the enclosure of a value or set
of values as computed with ball arithmetic will be larger than the smallest
possible enclosure.

Overestimation results naturally from rounding errors and cancellations
in the individual steps of a calculation.
As a general principle, formula rewriting techniques that make
floating-point code more numerically stable also make ball arithmetic code
more numerically stable, in the sense of producing tighter enclosures.

As a result of the *dependency problem*, ball or interval
arithmetic can produce error
bounds that are much larger than the actual numerical errors
resulting from doing floating-point arithmetic.
Consider the expression `(x + 1) - x` as an example.
When evaluated in floating-point
arithmetic, `x` may have a large initial error. However, that error will
cancel itself out in the subtraction, so that the result equals 1
(except perhaps for a small rounding error left from the operation `x + 1`).
In ball arithmetic, dependent errors add up instead of cancelling out.
If `x = [3 \pm 0.1]`, the result will be `[1 \pm 0.2]`, where
the error bound has doubled.
In unfavorable circumstances, error bounds can grow exponentially
with the number of steps.

If all inputs to a calculation are "point values", i.e.
exact numbers and known mathematical constants that can
be approximated arbitrarily closely (such as `\pi`), then an error
of order `2^n` can typically be overcome by working with *n* extra bits of
precision, increasing the computation time by an amount
that is polynomial in *n*.
In certain situations, however, overestimation leads to exponential
slowdown or even failure of an algorithm to converge.
For example, root-finding algorithms that refine the result iteratively
may fail to converge in ball arithmetic, even if they do converge in plain
floating-point arithmetic.

Therefore, ball arithmetic is not a silver bullet: there will always
be situations where some amount of numerical or mathematical analysis
is required. Some experimentation may be required to find whether
(and how) it can be used effectively for a given problem.

Predicates
-------------------------------------------------------------------------------

A ball implementation of a predicate 
`f : \mathbb{R} \to \{\operatorname{True}, \operatorname{False}\}`
would need to be able to return a third logical value indicating
that the result could be either True or False.
In most cases, predicates in Arb are implemented as 
functions that return the *int* value 1 to indicate that the
result certainly is True, and the *int* value 0 to indicate
that the result could be either True or False.
To test whether a predicate certainly is False, the user must
test whether the negated predicate certainly is True.

For example, the following code would *not* be correct in general:

.. code-block:: c

    if (arb_is_positive(x))
    {
        ...  /* do things assuming that x > 0 */
    }
    else
    {
        ...  /* do things assuming that x <= 0 */
    }

Instead, the following can be used:

.. code-block:: c

    if (arb_is_positive(x))
    {
        ...  /* do things assuming that x > 0 */
    }
    else if (arb_is_nonpositive(x))
    {
        ...  /* do things assuming that x <= 0 */
    }
    else
    {
        ... /* do things assuming that the sign of x is unknown */
    }

Likewise, we will write `x \le y` in mathematical notation with the meaning
that `x \le y` holds for all `x \in X, y \in Y` where `X` and `Y` are balls.

Note that some predicates such as :func:`arb_overlaps` and :func:`arb_contains`
actually are predicates on balls viewed as sets, and not ball implementations
of pointwise predicates.

Some predicates are also complementary.
For example :func:`arb_contains_zero` tests whether the input ball
contains the point zero.
Negated, it is equivalent to :func:`arb_is_nonzero`,
and complementary to :func:`arb_is_zero` as a pointwise predicate:

.. code-block:: c

    if (arb_is_zero(x))
    {
        ...  /* do things assuming that x = 0 */
    }
    #if 1
    else if (arb_is_nonzero(x))
    #else
    else if (!arb_contains_zero(x))      /* equivalent */
    #endif
    {
        ...  /* do things assuming that x != 0 */
    }
    else
    {
        ... /* do things assuming that the sign of x is unknown */
    }

A worked example: the sine function
-------------------------------------------------------------------------------

We implement the function `\sin(x)` naively using
the Taylor series `\sum_{k=0}^{\infty} (-1)^k x^{2k+1} / (2k+1)!`
and :type:`arb_t` arithmetic.
Since there are infinitely many terms, we need to split the series
in two parts: a finite sum that can be evaluated directly, and
a tail that has to be bounded.

We stop as soon as we reach a term `t` bounded by `|t| \le 2^{-prec} < 1`.
The terms are alternating and must have decreasing magnitude
from that point, so the tail of the series
is bounded by `|t|`. We add this magnitude to the radius
of the output. Since ball arithmetic automatically bounds the numerical errors
resulting from all arithmetic operations, the output *res* is a
ball guaranteed to contain `\sin(x)`.

.. code-block:: c

    #include "arb.h"

    void arb_sin_naive(arb_t res, const arb_t x, slong prec)
    {
        arb_t s, t, u, tol;
        slong k;
        arb_init(s); arb_init(t); arb_init(u); arb_init(tol);

        arb_one(tol);
        arb_mul_2exp_si(tol, tol, -prec);  /* tol = 2^-prec */

        for (k = 0; ; k++)
        {
            arb_pow_ui(t, x, 2 * k + 1, prec);
            arb_fac_ui(u, 2 * k + 1, prec);
            arb_div(t, t, u, prec);  /* t = x^(2k+1) / (2k+1)! */

            arb_abs(u, t);
            if (arb_le(u, tol))   /* if |t| <= 2^-prec */
            {
                arb_add_error(s, u);    /* add |t| to the radius and stop */
                break;
            }

            if (k % 2 == 0)
                arb_add(s, s, t, prec);
            else
                arb_sub(s, s, t, prec);

        }

        arb_set(res, s);
        arb_clear(s); arb_clear(t); arb_clear(u); arb_clear(tol);
    }

This algorithm is naive, because the Taylor series is slow to converge
and suffers from catastrophic cancellation when `|x|` is large
(we could also improve the efficiency of the code slightly by
computing the terms using recurrence relations instead of
computing `x^k` and `k!` from scratch each iteration).

As a test, we compute `\sin(2016.1)`.
The largest term in the Taylor series for `\sin(x)` reaches
a magnitude of about `x^x / x!`, or about `10^{873}` in this case.
Therefore, we need over 873 digits (about 3000 bits) of precision
to overcome the catastrophic cancellation and determine
the result with sufficient accuracy to tell whether it is positive
or negative.

.. code-block:: c

    int main()
    {
        arb_t x, y;
        slong prec;
        arb_init(x); arb_init(y);

        for (prec = 64; ; prec *= 2)
        {
            arb_set_str(x, "2016.1", prec);
            arb_sin_naive(y, x, prec);
            printf("Using %5ld bits, sin(x) = ", prec);
            arb_printn(y, 10, 0); printf("\n");
            if (!arb_contains_zero(y))  /* stopping condition */
                break;
        }

        arb_clear(x); arb_clear(y);
    }

The program produces the following output:

.. code-block:: text

    Using    64 bits, sin(x) = [+/- 2.67e+859]
    Using   128 bits, sin(x) = [+/- 1.30e+840]
    Using   256 bits, sin(x) = [+/- 3.60e+801]
    Using   512 bits, sin(x) = [+/- 3.01e+724]
    Using  1024 bits, sin(x) = [+/- 2.18e+570]
    Using  2048 bits, sin(x) = [+/- 1.22e+262]
    Using  4096 bits, sin(x) = [-0.7190842207 +/- 1.20e-11]

As an exercise, the reader may improve the naive algorithm by making it
subtract a well-chosen multiple of `2 \pi` from `x` before invoking
the Taylor series (hint: use :func:`arb_const_pi`, :func:`arb_div`
and :func:`arf_get_fmpz`).
If done correctly, 64 bits of precision should be more than enough to
compute `\sin(2016.1)`, and with minor adjustments
to the code, the user should be able to compute
`\sin(\exp(2016.1))` quite easily as well.

This example illustrates how ball arithmetic can be used to perform
nontrivial calculations. To evaluate an infinite series, the user
needs to know how to bound the tail of the series, but everything
else is automatic.
When evaluating a finite formula that can be expressed
completely using built-in functions, all error bounding is automatic
from the point of view of the user.
In particular, the :func:`arb_sin` method should be used to compute the sine
of a real number; it uses a much more efficient algorithm
than the naive code above.

This example also illustrates the "guess-and-verify" paradigm:
instead of determining *a priori* the floating-point precision necessary
to get a correct result, we *guess* some initial precision, use ball arithmetic
to *verify* that the result is accurate enough, and restart with
higher precision (or signal failure) if it is not.

If we think of rounding errors as essentially random processes,
then a floating-point computation is analogous to a
*Monte Carlo algorithm*. Using ball arithmetic to get a verified result
effectively turns it into the analog of a *Las Vegas algorithm*,
which is a randomized algorithm that always gives a correct result if it terminates, but
may fail to terminate (alternatively, instead of actually looping forever,
it might signal failure after a certain number of iterations).

The loop will fail to terminate if we attempt to determine the sign of
`\sin(\pi)`:

.. code-block:: text

    Using    64 bits, sin(x) = [+/- 3.96e-18]
    Using   128 bits, sin(x) = [+/- 2.17e-37]
    Using   256 bits, sin(x) = [+/- 6.10e-76]
    Using   512 bits, sin(x) = [+/- 5.13e-153]
    Using  1024 bits, sin(x) = [+/- 4.01e-307]
    Using  2048 bits, sin(x) = [+/- 2.13e-615]
    Using  4096 bits, sin(x) = [+/- 6.85e-1232]
    Using  8192 bits, sin(x) = [+/- 6.46e-2465]
    Using 16384 bits, sin(x) = [+/- 5.09e-4931]
    Using 32768 bits, sin(x) = [+/- 5.41e-9863]
    ...

The sign of a nonzero real number can be
decided by computing it to sufficiently high accuracy, but the sign
of an expression that is exactly equal to zero cannot be decided
by a numerical computation unless the entire computation
happens to be exact (in this example, we could use the :func:`arb_sin_pi` 
function which computes `\sin(\pi x)` in one step, with the input `x = 1`).

It is up to the user to implement a stopping criterion appropriate for
the circumstances of a given application. For example, breaking
when it is clear that `|\sin(x)| < 10^{-10000}` would allow the program
to terminate and convey some meaningful information about the input `x = \pi`,
though this would not constitute a mathematical proof that
`\sin(\pi) = 0`.

More on precision and accuracy
-------------------------------------------------------------------------------

The relation between the working precision and the accuracy of the output
is not always easy predict. The following remarks might help
to choose *prec* optimally.

For a ball `[m \pm r]` it is convenient to define the following notions:

* Absolute error: `e_{abs} = |r|`
* Relative error: `e_{rel} = |r| / \max(0, |m| - |r|)` (or `e_{rel} = 0` if `r = m = 0`)
* Absolute accuracy: `a_{abs} = 1 / e_{abs}`
* Relative accuracy: `a_{rel} = 1 / e_{rel}`

Expressed in bits, one takes the corresponding `\log_2` values.

Of course, if `x` is the exact value being approximated, then
the "absolute error" so defined is an upper bound for the
actual absolute error `|x-m|` and "absolute accuracy"
a lower bound for `1/|x-m|`, etc.

The *prec* argument in Arb should be thought of as controlling
the working precision.
Generically, when evaluating a fixed expression (that is, when the
sequence of operations does not depend on the precision), the
absolute or relative error will be bounded by

.. math::

    2^{O(1) - prec}

where the `O(1)` term depends on the expression and implementation
details of the ball functions used to evaluate it.
Accordingly, for an accuracy of *p* bits, we need to use a working precision
`O(1) + p`.
If the expression is numerically well-behaved, then the `O(1)` term
will be small, which leads to the heuristic of "adding a few guard bits"
(for most basic calculations, 10 or 20 guard bits is enough).
If the `O(1)` term is unknown, then increasing the number of guard
bits in exponential steps until the result is accurate enough
is generally a good heuristic.

Sometimes, a partially accurate result can be used to estimate the `O(1)`
term. For example, if the goal is to achieve 100 bits of accuracy
and a precision of 120 bits yields 80 bits of accuracy, then
it is plausible that a precision of just over
140 bits yields 100 bits of accuracy.

Built-in functions in Arb can roughly be characterized as
belonging to one of two extremes (though there is actually a spectrum):

* Simple operations, including basic arithmetic operations and many
  elementary functions. In most cases, for an input `x = [m \pm r]`,
  `f(x)` is evaluated by computing `f(m)` and then separately bounding the
  *propagated error* `|f(m) - f(m + \varepsilon)|, |\varepsilon| \le r`.
  The working precision is automatically increased internally
  so that `f(m)` is computed to *prec* bits of relative accuracy
  with an error of at most a few units in the last place (perhaps with
  rare exceptions).
  The propagated error can generally be bounded quite tightly as well (see :ref:`general_formulas`).
  As a result, the enclosure will be close to the best possible
  at the given precision, and the user can estimate the precision to use
  accordingly.

* Complex operations, such as certain higher
  transcendental functions (for example, the Riemann zeta function).
  The function is evaluated by performing a sequence of simpler operations,
  each using ball arithmetic with a working precision of roughly *prec*
  bits. The sequence of operations might depend on *prec*;
  for example, an infinite series might be truncated
  so that the remainder is smaller than `2^{-prec}`.
  The final result can be far from tight, and it is not guaranteed
  that the error converges to zero as `prec \to \infty`, though
  in practice, it should do so in most cases.

In short, the *inclusion principle* is the fundamental contract in Arb.
Enclosures computed by built-in functions may or may not be tight
enough to be useful, but the hope is that they will be sufficient
for most purposes.
Tightening the error bounds for more complex operations is a long
term optimization goal, which in many cases will require a
fair amount of research.
A tradeoff also has to be made for efficiency: tighter error bounds
allow the user to work with a lower precision, but they may
also be much more expensive to compute.

Polynomial time guarantee
-------------------------------------------------------------------------------

Arb provides a soft guarantee that the time used to evaluate a ball
function will depend polynomially on *prec* and the bit size
of the input, uniformly regardless of the numerical value of the input.

The idea behind this soft guarantee is to allow Arb to be used as a
black box to evaluate expressions numerically without potentially
slowing down, hanging indefinitely or crashing
because of "bad" input such as nested exponentials.
By controlling the precision, the user can cancel
a computation before it uses up
an unreasonable amount of resources,
without having to rely on other timeout or exception mechanisms.
A result that is feasible but very expensive to compute
can still be forced by setting the precision high enough.

As motivation, consider evaluating `\sin(x)` or `\exp(x)` with
the exact floating-point number
`x = 2^{2^n}` as input.
The time and space required to compute an accurate floating-point
approximation of `\sin(x)` or `\exp(x)` increases as `2^n`,
in the first case because because of the need to subtract an accurate
multiple of `2\pi` and in the second case due to the size of the
output exponent and the internal subtraction of an accurate multiple of `\log(2)`.
This is despite the fact that the size of `x` as an object in memory only
increases linearly with `n`.
Already `n = 33` would require at least 1 GB of memory, and
`n = 100` would be physically impossible to process.
For functions that are computed by direct use of power series expansions,
e.g. `f(x) = \sum_{k=0}^{\infty} c_k x^k`,
without having fast argument-reduction techniques
like those for elementary functions,
the time would be exponential in `n` already when `x = 2^n`.

Therefore, Arb caps internal work parameters
(the internal working precision,
the number terms of an infinite series to add, etc.) by polynomial,
usually linear, functions of *prec*.
When the limit is exceeded, the output is set to a crude bound.
For example, if `x` is too large, :func:`arb_sin` will
simply return `[\pm 1]`, and :func:`arb_exp`
will simply return `[\pm \infty]` if `x` is positive
or `[\pm 2^{-m}]` if `x` is negative.

This is not just a failsafe, but occasionally a useful optimization.
It is not entirely uncommon to have formulas where one term
is modest and another term decreases exponentially, such as:

.. math::

    \log(x) + \sin(x) \exp(-x).

For example, the reflection formula of the digamma function has
a similar structure.
When `x` is large, the right term would be expensive to compute
to high relative accuracy. Doing so is unnecessary, however,
since a crude bound of `[\pm 1] \cdot [\pm 2^{-m}]` is enough to evaluate
the expression as a whole accurately.

The polynomial time guarantee is "soft" in that there are a few exceptions.
For example, the complexity of computing the Riemann zeta function
`\zeta(\sigma+it)` increases linearly with the imaginary height `|t|`
in the current implementation, and all known algorithms
have a complexity of `|t|^{\alpha}` where the best known value for `\alpha`
is about `0.3`.
Input with large `|t|` is most likely to be given deliberately
by users with the explicit intent of evaluating the zeta
function itself, so the evaluation is not cut off automatically.

