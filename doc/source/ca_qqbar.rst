.. _ca_qqbar:

**ca_qqbar.h** -- algebraic numbers
===============================================================================

A :type:`ca_qqbar_t` represents a real or complex algebraic number
(an element of `\overline{\mathbb{Q}}`) by its unique reduced
minimal polynomial in `\mathbb{Z}[x]` and an isolating complex interval.
The precision of isolating intervals is maintained automatically to
ensure that all operations on :type:`ca_qqbar_t` instances are exact.

This type is useful for representing individual algebraic numbers
of moderate degree (up to 100, say). An arithmetic operation
on numbers of degrees *m* and *n* involves computing
and then factoring an annihilating polynomial of degree *mn*;
this will generally be the bottleneck.
In rare cases, numerical root-finding can be the more expensive step.
For doing lots of arithmetic with algebraic numbers, it is
often a good idea to work in a fixed number field instead of
doing repeated operations on :type:`ca_qqbar_t` instances.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: ca_qqbar_struct

.. type:: ca_qqbar_t

    A *ca_qqbar_struct* consists of an *fmpz_poly_struct* and an *acb_struct*.
    A *ca_qqbar_t* is defined as an array of length one of type
    *ca_qqbar_struct*, permitting a *ca_qqbar_t* to be passed by
    reference.

.. macro:: CA_QQBAR_POLY(x)

    Macro returning a pointer to the minimal polynomial of *x* as an *fmpz_poly_t*.

.. macro:: CA_QQBAR_COEFFS(x)

    Macro returning a pointer to the array of *fmpz* coefficients of the
    minimal polynomial of *x*.

.. macro:: CA_QQBAR_ENCLOSURE(x)

    Macro returning a pointer to the enclosure of *x* as an *acb_t*.

Memory management
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_init(ca_qqbar_t res)

    Initializes the variable *res* for use, and sets its value to zero.

.. function:: void ca_qqbar_clear(ca_qqbar_t res)

    Clears the variable *res*, freeing or recycling its allocated memory.

Assignment
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_swap(ca_qqbar_t x, ca_qqbar_t y)

    Swaps the values of *x* and *y* efficiently.

.. function:: void ca_qqbar_set(ca_qqbar_t res, const ca_qqbar_t x)

.. function:: void ca_qqbar_set_si(ca_qqbar_t res, slong x)

.. function:: void ca_qqbar_set_ui(ca_qqbar_t res, ulong x)

.. function:: void ca_qqbar_set_fmpz(ca_qqbar_t res, const fmpz_t x)

.. function:: void ca_qqbar_set_fmpq(ca_qqbar_t res, const fmpq_t x)

    Sets *res* to the value *x*.

Properties
-------------------------------------------------------------------------------

.. function:: slong ca_qqbar_degree(const ca_qqbar_t x)

    Returns the degree of *x*, i.e. the degree of the minimal polynomial.

.. function:: int ca_qqbar_is_rational(const ca_qqbar_t x)

    Returns whether *x* is a rational number.

.. function:: int ca_qqbar_is_integer(const ca_qqbar_t x)

    Returns whether *x* is an integer (an element of `\mathbb{Z}`).

.. function:: int ca_qqbar_is_algebraic_integer(const ca_qqbar_t x)

    Returns whether *x* is an algebraic integer, i.e. whether its minimal
    polynomial has leading coefficient 1.

.. function:: int ca_qqbar_is_zero(const ca_qqbar_t x)

    Returns whether *x* is the number 0.

.. function:: int ca_qqbar_is_one(const ca_qqbar_t x)

    Returns whether *x* is the number 1.

.. function:: int ca_qqbar_is_real(const ca_qqbar_t x)

    Returns whether *x* is a real number.

.. function:: int ca_qqbar_real_sgn(const ca_qqbar_t x)

    Returns the sign of the real part of *x* (-1, 0 or +1).

.. function:: int ca_qqbar_imag_sgn(const ca_qqbar_t x)

    Returns the sign of the imaginary part of *x* (-1, 0 or +1).

Special values
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_zero(ca_qqbar_t res)

    Sets *res* to the number 0.

.. function:: void ca_qqbar_one(ca_qqbar_t res)

    Sets *res* to the number 1.

.. function:: void ca_qqbar_i(ca_qqbar_t res)

    Sets *res* to the imaginary unit `i`.

.. function:: void ca_qqbar_phi(ca_qqbar_t res)

    Sets *res* to the golden ratio `\varphi = \tfrac{1}{2}(\sqrt{5} + 1)`.

Input and output
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_print(const ca_qqbar_t x)

    Prints *res* to standard output. The output shows the list of coefficients
    of the minimal polynomial followed by a decimal representation of
    the enclosing interval. This function is mainly intended for debugging.

Random generation
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_randtest(ca_qqbar_t res, flint_rand_t state, slong deg, slong bits)

    Sets *res* to a random algebraic number with degree up to *deg* and
    with height (measured in bits) up to *bits*.

.. function:: void ca_qqbar_randtest_real(ca_qqbar_t res, flint_rand_t state, slong deg, slong bits)

    Sets *res* to a random real algebraic number with degree up to *deg* and
    with height (measured in bits) up to *bits*.

.. function:: void ca_qqbar_randtest_nonreal(ca_qqbar_t res, flint_rand_t state, slong deg, slong bits)

    Sets *res* to a random nonreal algebraic number with degree up to *deg* and
    with height (measured in bits) up to *bits*. Since all algebraic numbers
    of degree 1 are real, *deg* must be at least 2.

Comparisons
-------------------------------------------------------------------------------

.. function:: int ca_qqbar_equal(const ca_qqbar_t x, const ca_qqbar_t y)

    Returns whether *x* and *y* are equal.

.. function:: int ca_qqbar_cmp_re(const ca_qqbar_t x, const ca_qqbar_t y)

    Compares the real parts of *x* and *y*, returning -1, 0 or +1.

.. function:: int ca_qqbar_cmp_im(const ca_qqbar_t x, const ca_qqbar_t y)

    Compares the imaginary parts of *x* and *y*, returning -1, 0 or +1.

.. function:: int ca_qqbar_cmpabs_re(const ca_qqbar_t x, const ca_qqbar_t y)

    Compares the absolute values of the real parts of *x* and *y*, returning -1, 0 or +1.

.. function:: int ca_qqbar_cmpabs_im(const ca_qqbar_t x, const ca_qqbar_t y)

    Compares the absolute values of the imaginary parts of *x* and *y*, returning -1, 0 or +1.

.. function:: int ca_qqbar_cmpabs(const ca_qqbar_t x, const ca_qqbar_t y)

    Compares the absolute values of *x* and *y*, returning -1, 0 or +1.

Complex parts
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_conj(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the complex conjugate of *x*.

.. function:: void ca_qqbar_re(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the real part of *x*.

.. function:: void ca_qqbar_im(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the imaginary part of *x*.

.. function:: void ca_qqbar_re_im(ca_qqbar_t res1, ca_qqbar_t res2, const ca_qqbar_t x)

    Sets *res1* to the real part of *x* and *res2* to the imaginary part of *x*.

.. function:: void ca_qqbar_abs(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the absolute value of *x*:

.. function:: void ca_qqbar_abs2(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the square of the absolute value of *x*.

.. function:: void ca_qqbar_sgn(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the complex sign of *x*, defined as 0 if *x* is zero
    and as `x / |x|` otherwise.

.. function:: int ca_qqbar_csgn(const ca_qqbar_t x)

    Returns the extension of the real sign function taking the
    value 1 for *x* strictly in the right half plane, -1 for *x* strictly
    in the left half plane, and the sign of the imaginary part when *x* is on
    the imaginary axis. Equivalently, `\operatorname{csgn}(x) = x / \sqrt{x^2}`
    except that the value is 0 when *x* is zero.

Integer parts
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_floor(fmpz_t res, const ca_qqbar_t x)

    Sets *res* to the floor function of *x*. If *x* is not real, the
    value is defined as the floor function of the real part of *x*.

.. function:: void ca_qqbar_ceil(fmpz_t res, const ca_qqbar_t x)

    Sets *res* to the ceiling function of *x*. If *x* is not real, the
    value is defined as the ceiling function of the real part of *x*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_neg(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the negation of *x*.

.. function:: void ca_qqbar_add(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_add_fmpq(ca_qqbar_t res, const ca_qqbar_t x, const fmpq_t y)

.. function:: void ca_qqbar_add_fmpz(ca_qqbar_t res, const ca_qqbar_t x, const fmpz_t y)

.. function:: void ca_qqbar_add_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong y)

.. function:: void ca_qqbar_add_si(ca_qqbar_t res, const ca_qqbar_t x, slong y)

    Sets *res* to the sum of *x* and *y*.

.. function:: void ca_qqbar_sub(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_sub_fmpq(ca_qqbar_t res, const ca_qqbar_t x, const fmpq_t y)

.. function:: void ca_qqbar_sub_fmpz(ca_qqbar_t res, const ca_qqbar_t x, const fmpz_t y)

.. function:: void ca_qqbar_sub_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong y)

.. function:: void ca_qqbar_sub_si(ca_qqbar_t res, const ca_qqbar_t x, slong y)

    Sets *res* to the difference of *x* and *y*.

.. function:: void ca_qqbar_mul(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_mul_fmpq(ca_qqbar_t res, const ca_qqbar_t x, const fmpq_t y)

.. function:: void ca_qqbar_mul_fmpz(ca_qqbar_t res, const ca_qqbar_t x, const fmpz_t y)

.. function:: void ca_qqbar_mul_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong y)

.. function:: void ca_qqbar_mul_si(ca_qqbar_t res, const ca_qqbar_t x, slong y)

    Sets *res* to the product of *x* and *y*.

.. function:: void ca_qqbar_mul_2exp_si(ca_qqbar_t res, const ca_qqbar_t x, slong e)

    Sets *res* to *x* multiplied by `2^e`.

.. function:: void ca_qqbar_sqr(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the square of *x*.

.. function:: void ca_qqbar_inv(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

    Sets *res* to the multiplicative inverse of *y*.
    Division by zero calls *flint_abort*.

.. function:: void ca_qqbar_div(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_div_fmpq(ca_qqbar_t res, const ca_qqbar_t x, const fmpq_t y)

.. function:: void ca_qqbar_div_fmpz(ca_qqbar_t res, const ca_qqbar_t x, const fmpz_t y)

.. function:: void ca_qqbar_div_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong y)

.. function:: void ca_qqbar_div_si(ca_qqbar_t res, const ca_qqbar_t x, slong y)

.. function:: void ca_qqbar_fmpq_div(ca_qqbar_t res, const fmpq_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_fmpz_div(ca_qqbar_t res, const fmpz_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_ui_div(ca_qqbar_t res, ulong x, const ca_qqbar_t y)

.. function:: void ca_qqbar_si_div(ca_qqbar_t res, slong x, const ca_qqbar_t y)

    Sets *res* to the quotient of *x* and *y*.
    Division by zero calls *flint_abort*.

.. function:: void ca_qqbar_affine_transform(ca_qqbar_t res, const ca_qqbar_t x, const fmpz_t a, const fmpz_t b, const fmpz_t c)

    Sets *res* to the rational affine transformation `(ax+b)/c`, performed as
    a single operation. There are no restrictions on *a*, *b* and *c*
    except that *c* must be nonzero. Division by zero calls *flint_abort*.

.. function:: void ca_qqbar_sqrt(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the principal square root of *x*.

.. function:: void ca_qqbar_rsqrt(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the reciprocal of the principal square root of *x*.
    Division by zero calls *flint_abort*.

.. function:: void ca_qqbar_pow_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong n)

    Sets *res* to *x* raised to the *n*-th power.

.. function:: void ca_qqbar_root_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong n)

    Sets *res* to the principal *n*-th root of *x*. The order *n*
    must be positive.

Numerical enclosures
-------------------------------------------------------------------------------

The following functions guarantee a polished output in which both the real
and imaginary parts are accurate to *prec* bits and exact when exactly
representable (that is, when a real or imaginary part is a sufficiently
small dyadic number).

In some cases, the computations needed to polish the output may be
expensive. When polish is unnecessary, :func:`ca_qqbar_enclosure_raw`
can be used instead.

.. function:: void ca_qqbar_get_acb(acb_t res, const ca_qqbar_t x, slong prec)

    Sets *res* to an enclosure of *x* rounded to *prec* bits.

.. function:: void ca_qqbar_get_arb(arb_t res, const ca_qqbar_t x, slong prec)

    Sets *res* to an enclosure of *x* rounded to *prec* bits, assuming that
    *x* is a real number. If *x* is not real, *res* is set to
    `[\operatorname{NaN} \pm \infty]`.

.. function:: void ca_qqbar_get_arb_re(arb_t res, const ca_qqbar_t x, slong prec)

    Sets *res* to an enclosure of the real part of *x* rounded to *prec* bits.

.. function:: void ca_qqbar_get_arb_im(arb_t res, const ca_qqbar_t x, slong prec)

    Sets *res* to an enclosure of the imaginary part of *x* rounded to *prec* bits.


Internal functions
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_fmpz_poly_composed_op(fmpz_poly_t res, const fmpz_poly_t A, const fmpz_poly_t B, int op)

    Given nonconstant polynomials *A* and *B*, sets *res* to a polynomial
    whose roots are `a+b`, `a-b`, `ab` or `a/b` for all roots *a* of *A*
    and all roots *b* of *B*. The parameter *op* selects the arithmetic
    operation: 0 for addition, 1 for subtraction, 2 for multiplication
    and 3 for division. If *op* is 3, *B* must not have zero as a root.

.. function:: void ca_qqbar_binary_op(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y, int op)

    Performs a binary operation using a generic algorithm. This does not
    check for special cases.

.. function:: int _ca_qqbar_validate_enclosure(acb_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec)

    Given *z* known to be an enclosure of at least one root of *poly*,
    certifies that the enclosure contains a unique root, and in that
    case sets *res* to a new (possibly improved) enclosure for the same
    root, returning 1. Returns 0 if uniqueness cannot be certified.

    The enclosure is validated by performing a single step with the
    interval Newton method. The working precision is determined from the
    accuracy of *z*, but limited by *max_prec* bits.

.. function:: void _ca_qqbar_enclosure_raw(acb_t res, const fmpz_poly_t poly, const acb_t z, slong prec)

.. function:: void ca_qqbar_enclosure_raw(acb_t res, const ca_qqbar_t x, slong prec)

    Sets *res* to an enclosure of *x* accurate to about *prec* bits
    (the actual accuracy can be slightly lower, or higher).

    This function uses repeated interval Newton steps to polish the initial
    enclosure *z*, doubling the working precision each time. If any step
    fails to improve the accuracy significantly, the root is recomputed
    from scratch to higher precision.

    If the initial enclosure is accurate enough, *res* is set to this value
    without rounding and without further computation.

