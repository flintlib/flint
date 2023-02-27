.. _arb_fmpz_poly:

**arb_fmpz_poly.h** -- extra methods for integer polynomials
=========================================================================================

This module provides methods for FLINT polynomials
with integer and rational coefficients
(:type:`fmpz_poly_t`) and (:type:`fmpq_poly_t`)
requiring use of Arb real or complex numbers.

Some methods output real or complex numbers while others
use real and complex numbers internally to produce an exact result.
This module also contains some useful helper functions not specifically
related to real and complex numbers.

Note that methods that combine Arb *polynomials* and FLINT polynomials
are found in the respective Arb polynomial modules, such as
:func:`arb_poly_set_fmpz_poly`
and :func:`arb_poly_get_unique_fmpz_poly`.

Evaluation
-------------------------------------------------------------------------------

.. function:: void _arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)

    Evaluates *poly* (given by a polynomial object or an array with *len* coefficients)
    at the given real or complex number, respectively using Horner's rule, rectangular
    splitting, or a default algorithm choice.

Utility methods
-------------------------------------------------------------------------------

.. function:: ulong arb_fmpz_poly_deflation(const fmpz_poly_t poly)

    Finds the maximal exponent by which *poly* can be deflated.

.. function:: void arb_fmpz_poly_deflate(fmpz_poly_t res, const fmpz_poly_t poly, ulong deflation)

    Sets *res* to a copy of *poly* deflated by the exponent *deflation*.

Polynomial roots
-------------------------------------------------------------------------------

.. function:: void arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, slong prec)

    Writes to *roots* all the real and complex roots of the polynomial *poly*,
    computed to at least *prec* accurate bits.
    The root enclosures are guaranteed to be disjoint, so that
    all roots are isolated.

    The real roots are written first in ascending order (with
    the imaginary parts set exactly to zero). The following
    nonreal roots are written in arbitrary order, but with conjugate pairs
    grouped together (the root in the upper plane leading
    the root in the lower plane).

    The input polynomial *must* be squarefree. For a general polynomial,
    compute the squarefree part `f / \gcd(f,f')` or do a full squarefree
    factorization to obtain the multiplicities of the roots::

        fmpz_poly_factor_t fac;
        fmpz_poly_factor_init(fac);
        fmpz_poly_factor_squarefree(fac, poly);

        for (i = 0; i < fac->num; i++)
        {
            deg = fmpz_poly_degree(fac->p + i);
            flint_printf("%wd roots of multiplicity %wd\n", deg, fac->exp[i]);
            roots = _acb_vec_init(deg);
            arb_fmpz_poly_complex_roots(roots, fac->p + i, 0, prec);
            _acb_vec_clear(roots, deg);
        }

        fmpz_poly_factor_clear(fac);

    All roots are refined to a relative accuracy of at least *prec* bits.
    The output values will generally have higher actual precision,
    depending on the precision needed for isolation and the
    precision used internally by the algorithm.

    This implementation should be adequate for general use, but it is not
    currently competitive with state-of-the-art isolation
    methods for finding real roots alone.

    The following *flags* are supported:

    * *ARB_FMPZ_POLY_ROOTS_VERBOSE*

Special polynomials
-------------------------------------------------------------------------------

Note: see also the methods available in FLINT (e.g. for cyclotomic polynomials).

.. function:: void arb_fmpz_poly_cos_minpoly(fmpz_poly_t res, ulong n)

    Sets *res* to the monic minimal polynomial of `2 \cos(2 \pi / n)`.
    This is a wrapper of FLINT's *fmpz_poly_cos_minpoly*, provided here
    for backward compatibility.

.. function:: void arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, ulong q, ulong n)

    Sets *res* to the minimal polynomial of the Gaussian periods
    `\sum_{a \in H} \zeta^a` where `\zeta = \exp(2 \pi i / q)`
    and *H* are the cosets of the subgroups of order `d = (q - 1) / n` of
    `(\mathbb{Z}/q\mathbb{Z})^{\times}`.
    The resulting polynomial has degree *n*.
    When `d = 1`, the result is the cyclotomic polynomial `\Phi_q`.

    The implementation assumes that *q* is prime, and that *n* is a divisor of
    `q - 1` such that *n* is coprime with *d*. If any condition is not met,
    *res* is set to the zero polynomial.

    This method provides a fast (in practice) way to
    construct finite field extensions of prescribed degree.
    If *q* satisfies the conditions stated above and `(q-1)/f` additionally
    is coprime with *n*, where *f* is the multiplicative order of *p* mod *q*, then
    the Gaussian period minimal polynomial is irreducible over
    `\operatorname{GF}(p)` [CP2005]_.

