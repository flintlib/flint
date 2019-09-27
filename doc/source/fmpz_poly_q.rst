.. _fmpz-poly-q:

**fmpz_poly_q.h** -- rational functions over the rational numbers
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_poly_q_struct

.. type:: fmpz_poly_q_t

    Description.

Memory management
--------------------------------------------------------------------------------

We represent a rational function over `\mathbf{Q}` as the quotient 
of two coprime integer polynomials of type ``fmpz_poly_t``, 
enforcing that the leading coefficient of the denominator is 
positive.  The zero function is represented as `0/1`.

.. function:: void fmpz_poly_q_init(fmpz_poly_q_t rop)

    Initialises ``rop``.

.. function:: void fmpz_poly_q_clear(fmpz_poly_q_t rop)

    Clears the object ``rop``.

.. function:: fmpz_poly_struct * fmpz_poly_q_numref(const fmpz_poly_q_t op)

    Returns a reference to the numerator of ``op``.

.. function:: fmpz_poly_struct * fmpz_poly_q_denref(const fmpz_poly_q_t op)

    Returns a reference to the denominator of ``op``.

.. function:: void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop)

    Brings ``rop`` into canonical form, only assuming that 
    the denominator is non-zero.

.. function:: int fmpz_poly_q_is_canonical(const fmpz_poly_q_t op)

    Checks whether the rational function ``op`` is in 
    canonical form.


Randomisation
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state, slong len1, flint_bitcnt_t bits1, slong len2, flint_bitcnt_t bits2)

    Sets ``poly`` to a random rational function.

.. function:: void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, slong len1, flint_bitcnt_t bits1, slong len2, flint_bitcnt_t bits2)

    Sets ``poly`` to a random non-zero rational function.


Assignment
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

    Sets the element ``rop`` to the same value as the element ``op``.

.. function:: void fmpz_poly_q_set_si(fmpz_poly_q_t rop, slong op)

    Sets the element ``rop`` to the value given by the ``slong`` 
    ``op``.

.. function:: void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2)

    Swaps the elements ``op1`` and ``op2``.

    This is done efficiently by swapping pointers.

.. function:: void fmpz_poly_q_zero(fmpz_poly_q_t rop)

    Sets ``rop`` to zero.

.. function:: void fmpz_poly_q_one(fmpz_poly_q_t rop)

    Sets ``rop`` to one.

.. function:: void fmpz_poly_q_neg(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

    Sets the element ``rop`` to the additive inverse of ``op``.

.. function:: void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

    Sets the element ``rop`` to the multiplicative inverse of ``op``.

    Assumes that the element ``op`` is non-zero.


Comparison
--------------------------------------------------------------------------------


.. function:: int fmpz_poly_q_is_zero(const fmpz_poly_q_t op)

    Returns whether the element ``op`` is zero.

.. function:: int fmpz_poly_q_is_one(const fmpz_poly_q_t op)

    Returns whether the element ``rop`` is equal to the constant 
    polynomial `1`.

.. function:: int fmpz_poly_q_equal(const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    Returns whether the two elements ``op1`` and ``op2`` are equal.


Addition and subtraction
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_q_add(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    Sets ``rop`` to the sum of ``op1`` and ``op2``.

.. function:: void fmpz_poly_q_sub(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    Sets ``rop`` to the difference of ``op1`` and ``op2``.

.. function:: void fmpz_poly_q_addmul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    Adds the product of ``op1`` and ``op2`` to ``rop``.

.. function:: void fmpz_poly_q_submul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    Subtracts the product of ``op1`` and ``op2`` from ``rop``.


Scalar multiplication and division
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x)

    Sets ``rop`` to the product of the rational function ``op`` 
    and the ``slong`` integer `x`.

.. function:: void fmpz_poly_q_scalar_mul_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x)

    Sets ``rop`` to the product of the rational function ``op`` 
    and the ``mpz_t`` integer `x`.

.. function:: void fmpz_poly_q_scalar_mul_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x)

    Sets ``rop`` to the product of the rational function ``op`` 
    and the ``mpq_t`` rational `x`.

.. function:: void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x)

    Sets ``rop`` to the quotient of the rational function ``op`` 
    and the ``slong`` integer `x`.

.. function:: void fmpz_poly_q_scalar_div_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x)

    Sets ``rop`` to the quotient of the rational function ``op`` 
    and the ``mpz_t`` integer `x`.

.. function:: void fmpz_poly_q_scalar_div_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x)

    Sets ``rop`` to the quotient of the rational function ``op`` 
    and the ``mpq_t`` rational `x`.


Multiplication and division
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_q_mul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    Sets ``rop`` to the product of ``op1`` and ``op2``.

.. function:: void fmpz_poly_q_div(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    Sets ``rop`` to the quotient of ``op1`` and ``op2``.


Powering
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_q_pow(fmpz_poly_q_t rop, const fmpz_poly_q_t op, ulong exp)

    Sets ``rop`` to the ``exp``-th power of ``op``.

    The corner case of ``exp == 0`` is handled by setting ``rop`` to 
    the constant function `1`.  Note that this includes the case `0^0 = 1`.


Derivative
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_q_derivative(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

    Sets ``rop`` to the derivative of ``op``.


Evaluation
--------------------------------------------------------------------------------


.. function:: int fmpz_poly_q_evaluate(mpq_t rop, const fmpz_poly_q_t f, const mpq_t a)

    Sets ``rop`` to `f` evaluated at the rational `a`.

    If the denominator evaluates to zero at `a`, returns non-zero and 
    does not modify any of the variables.  Otherwise, returns `0` and 
    sets ``rop`` to the rational `f(a)`.


Input and output
--------------------------------------------------------------------------------

The following three methods enable users to construct elements of type\\ 
``fmpz_poly_q_t`` from strings or to obtain string representations of 
such elements.
The format used is based on the FLINT format for integer polynomials of 
type ``fmpz_poly_t``, which we recall first: 
A non-zero polynomial `a_0 + a_1 X + \dotsb + a_n X^n` of length 
`n + 1` is represented by the string ``"n+1  a_0 a_1 ... a_n"``, 
where there are two space characters following the length and single 
space characters separating the individual coefficients.  There is no 
leading or trailing white-space.  The zero polynomial is simply 
represented by ``"0"``.
We adapt this notation for rational functions as follows.  We denote the 
zero function by ``"0"``.  Given a non-zero function with numerator 
and denominator string representations ``num`` and ``den``, 
respectively, we use the string ``num/den`` to represent the rational 
function, unless the denominator is equal to one, in which case we simply 
use ``num``.
There is also a ``_pretty`` variant available, which bases the string 
parts for the numerator and denominator on the output of the function 
``fmpz_poly_get_str_pretty`` and introduces parentheses where 
necessary.
Note that currently these functions are not optimised for performance and 
are intended to be used only for debugging purposes or one-off input and 
output, rather than as a low-level parser.

.. function:: int fmpz_poly_q_set_str(fmpz_poly_q_t rop, const char *s)

    Sets ``rop`` to the rational function given 
    by the string ``s``.

.. function:: char * fmpz_poly_q_get_str(const fmpz_poly_q_t op)

    Returns the string representation of 
    the rational function ``op``.

.. function:: char * fmpz_poly_q_get_str_pretty(const fmpz_poly_q_t op, const char *x)

    Returns the pretty string representation of 
    the rational function ``op``.

.. function:: int fmpz_poly_q_print(const fmpz_poly_q_t op)

    Prints the representation of the rational 
    function ``op`` to ``stdout``.

.. function:: int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x)

    Prints the pretty representation of the rational 
    function ``op`` to ``stdout``.
