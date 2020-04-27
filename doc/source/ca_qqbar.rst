.. _ca_qqbar:

**ca_qqbar.h** -- algebraic numbers
===============================================================================

A :type:`ca_qqbar_t` represents a real or complex algebraic number
(an element of `\overline{\mathbb{Q}}`) by its unique reduced
minimal polynomial in `\mathbb{Z}[x]` and an isolating complex interval.

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

    Returns whether *x* is an integer.

.. function:: int ca_qqbar_is_zero(const ca_qqbar_t x)

    Returns whether *x* is the number 0.

.. function:: int ca_qqbar_is_one(const ca_qqbar_t x)

    Returns whether *x* is the number 1.

Special values
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_zero(ca_qqbar_t res)

    Sets *res* to the number 0.

.. function:: void ca_qqbar_one(ca_qqbar_t res)

    Sets *res* to the number 1.

Random generation
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_randtest(ca_qqbar_t res, flint_rand_t state, slong deg, slong bits)

    Sets *res* to a random algebraic number with degree up to *deg* and
    with height (measured in bits) up to *bits*.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void ca_qqbar_conj(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the complex conjugate of *x*.

.. function:: void ca_qqbar_neg(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the negation of *x*.

.. function:: void ca_qqbar_add(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_sub(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_mul(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

.. function:: void ca_qqbar_div(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

    Sets *res* to the sum, difference, product or quotient of *x* and *y*.
    Division by zero calls *flint_abort*.

.. function:: void ca_qqbar_inv(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)

    Sets *res* to the multiplicative inverse of *y*.
    Division by zero calls *flint_abort*.

.. function:: void ca_qqbar_sqrt(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the principal square root of *x*.

.. function:: void ca_qqbar_rsqrt(ca_qqbar_t res, const ca_qqbar_t x)

    Sets *res* to the reciprocal of the principal square root of *x*.
    Division by zero calls *flint_abort*.

.. function:: void ca_qqbar_root_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong n)

    Sets *res* to the principal *n*-th root of *x*. The order *n*
    must be positive.

Internal polynomial and enclosure functions
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

