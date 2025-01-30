.. _issues:

Technical conventions and potential issues
===============================================================================

Integer overflow
-------------------------------------------------------------------------------

When machine-size integers are used for precisions, sizes of integers in
bits, lengths of polynomials, and similar quantities that relate
to sizes in memory, very few internal checks are performed to verify that
such quantities do not overflow.

Precisions and lengths exceeding a small fraction
of *LONG_MAX*, say `2^{24} \approx 10^7` on 32-bit systems,
should be regarded as resulting in undefined behavior.
On 64-bit systems this should generally not be an issue,
since most calculations will exhaust the available memory
(or the user's patience waiting for the computation to complete)
long before running into integer overflows.
However, the user needs to be wary of unintentionally passing input
parameters of order *LONG_MAX* or negative parameters where
positive parameters are expected, for example due to a runaway loop
that repeatedly increases the precision.

Currently, no hard upper limit on the precision is defined, but
`2^{24} \approx 10^7` bits on 32-bit system
and `2^{36} \approx 10^{11}` bits on a 64-bit system
can be considered safe for most purposes.
The relatively low limit on 64-bit systems is due to the fact that GMP
integers are used internally in some algorithms, and GMP integers
are limited to `2^{37}` bits.
The minimum allowed precision is 2 bits.

This caveat does not apply to exponents of floating-point numbers,
which are represented as arbitrary-precision integers, nor to
integers used as numerical scalars (e.g. :func:`arb_mul_si`).
However, it still applies to conversions and operations where
the result is requested exactly and sizes become an issue.
For example, trying to convert
the floating-point number `2^{2^{100}}` to an integer could
result in anything from a silent wrong value to thrashing followed
by a crash, and it is the user's responsibility not
to attempt such a thing.

Aliasing
-------------------------------------------------------------------------------

As a rule, Arb allows aliasing of operands. For example, in the function call
``arb_add(z, x, y, prec)``,
which performs `z \gets x + y`, any two (or all three) of the variables *x*,
*y* and *z* are allowed to be the same. Exceptions to this rule are
documented explicitly.

The general rule that input and output variables can be aliased with each
other only applies to variables *of the same type*
(ignoring *const* qualifiers on input variables -- a special case is that
:type:`arb_srcptr` is considered the *const* version of :type:`arb_ptr`).
This is a natural extension of the so-called *strict aliasing rule* in C.

For example, in :func:`arb_poly_evaluate` which evaluates
`y = f(x)` for a polynomial *f*, the output variable *y* is
not allowed to be a pointer to one of the coefficients of *f* (but
aliasing between *x* and *y* or between *x* and the coefficients
of *f* is allowed).
This also applies to :func:`_arb_poly_evaluate`:
for the purposes of aliasing,
:type:`arb_srcptr` (the type of the coefficient array within *f*) and :type:`arb_t`
(the type of *x*) are *not* considered
to be the same type, and therefore must not be aliased
with each other,
even though an :type:`arb_ptr`/:type:`arb_srcptr` variable pointing
to a length 1 array would otherwise be interchangeable with an :type:`arb_t`/*const* :type:`arb_t`.

Moreover, in functions that allow aliasing between an input
array and an output array, the arrays must either be identical or
completely disjoint, never partially overlapping.

There are natural exceptions to these aliasing restrictions, which may
used internally without being documented explicitly.
However, third party code should avoid relying on such exceptions.

An important caveat applies to **aliasing of input variables**.
Identical pointers are understood to
give permission for **algebraic simplification**.
This assumption is made to improve performance.
For example, the call ``arb_mul(z, x, x, prec)``
sets *z* to a ball enclosing the set

.. math::

    \{ t^2 \,:\, t \in x \}

and not the (generally larger) set

.. math::

    \{ t u \,:\, t \in x, u \in x \}.

If the user knows that two values *x* and *y*
both lie in the interval `[-1,1]` and wants to compute an
enclosure for `f(x,y)`, then it would be a mistake to 
create an :type:`arb_t` variable *x* enclosing `[-1,1]`
and reusing the same variable for *y*, calling `f(x,x)`.
Instead, the user has to create a
distinct variable *y* also enclosing `[-1,1]`.

Algebraic simplification is not guaranteed to occur.
For example, ``arb_add(z, x, x, prec)`` and ``arb_sub(z, x, x, prec)``
currently do not implement this optimization.
It is better to use ``arb_mul_2exp_si(z, x, 1)`` and
``arb_zero(z)``, respectively.

Thread safety and caches
-------------------------------------------------------------------------------

Arb should be fully threadsafe, provided that both MPFR and FLINT have
been built in threadsafe mode.
Use ``flint_set_num_threads()`` to set the number of threads that
Arb is allowed to use internally for single computations
(this is currently only exploited by a handful of operations).
Please note that thread safety is
only tested minimally, and extra caution when developing
multithreaded code is therefore recommended.

Arb may cache some data (such as the value of `\pi` and
Bernoulli numbers) to speed up various computations. In threadsafe mode,
caches use thread-local storage. There is currently no way to save memory
and avoid recomputation by having several threads share the same cache.
Caches can be freed by calling the ``flint_cleanup()`` function. To avoid
memory leaks, the user should call ``flint_cleanup()`` when exiting a thread.
It is also recommended to call ``flint_cleanup()`` when exiting the main
program (this should result in a clean output when running
`Valgrind <http://valgrind.org/>`_, and can help catching memory issues).

There does not seem to be an obvious way to make sure that ``flint_cleanup()``
is called when exiting a thread using OpenMP.
A possible solution to this problem is to use OpenMP sections,
or to use C++ and create a thread-local object whose destructor
invokes ``flint_cleanup()``.

Use of hardware floating-point arithmetic
-------------------------------------------------------------------------------

Arb uses hardware floating-point arithmetic (the ``double`` type in C) in two
different ways.

First, ``double`` arithmetic as well as transcendental ``libm`` functions
(such as ``exp``, ``log``) are used to select parameters heuristically
in various algorithms. Such heuristic use of approximate arithmetic does not
affect correctness: when any error bounds depend on the parameters, the error
bounds are evaluated separately using rigorous methods. At worst, flaws
in the floating-point arithmetic on a particular machine could cause an
algorithm to become inefficient due to inefficient parameters being
selected.

Second, ``double`` arithmetic is used internally for some rigorous error bound
calculations. To guarantee correctness, we make the following assumptions.
With the stated exceptions, these should hold on all commonly used platforms.

* A ``double`` uses the standard IEEE 754 format (with a 53-bit significand,
  11-bit exponent, encoding of infinities and NaNs, etc.)
* We assume that the compiler does not perform "unsafe" floating-point
  optimizations, such as reordering of operations. Unsafe optimizations are
  disabled by default in most modern C compilers, including GCC and Clang.
  The exception appears to be the Intel C++ compiler, which does some
  unsafe optimizations by default. These must be disabled by the user.
* We do not assume that floating-point operations are correctly rounded
  (a counterexample is the x87 FPU), or that rounding is done in any
  particular direction (the rounding mode may have been changed by the user).
  We assume that any floating-point operation is done with at most 1.1 ulp
  error.
* We do not assume that underflow or overflow behaves in a particular way (we
  only use doubles that fit in the regular exponent range, or explicit
  infinities).
* We do not use transcendental ``libm`` functions, since these can have errors
  of several ulps, and there is unfortunately no way to get guaranteed
  bounds. However, we do use functions such as ``ldexp`` and ``sqrt``, which we
  assume to be correctly implemented.

Interface changes
-------------------------------------------------------------------------------

Most of the core API should be stable at this point,
and significant compatibility-breaking changes will be specified in the
release notes.

In general, Arb does not distinguish between "private" and "public"
parts of the API. The implementation is meant to be transparent by design.
All methods are intended to be fully documented and tested
(exceptions to this are mainly due to lack of time on part of the
author).
The user should use common sense to determine whether a function is
concerned with implementation details, making it likely
to change as the implementation changes in the future.
The interface of :func:`arb_add` is probably not going to change in
the next version, but :func:`_arb_get_mpn_fixed_mod_pi4` just might.

General note on correctness
-------------------------------------------------------------------------------

Except where otherwise specified, Arb is designed to produce
provably correct error bounds. The code has been written carefully,
and the library is extensively tested.
However, like any complex mathematical software, Arb is virtually certain to
contain bugs, so the usual precautions are advised:

* Do sanity checks. For example, check that the result satisfies an expected
  mathematical relation, or compute the same result in two different ways,
  with different settings, and with different levels of precision.
  Arb's unit tests already do such checks, but they are not guaranteed to
  catch every possible bug, and they provide no protection against
  the user accidentally using the interface incorrectly.
* Compare results with other mathematical software.
* Read the source code to verify that it really does what it is supposed to do.

All bug reports are highly appreciated.

