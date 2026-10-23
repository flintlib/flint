/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"
#include "impl.h"

#define OP_ADD 0
#define OP_SUB 1
#define OP_MUL 2
#define OP_DIV 3

static void
fmpq_poly_hadamard_product(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    slong i, len;

    len = FLINT_MIN(fmpq_poly_length(poly1), fmpq_poly_length(poly2));

    fmpq_poly_fit_length(res, len);

    for (i = 0; i < len; i++)
        fmpz_mul(res->coeffs + i, poly1->coeffs + i, poly2->coeffs + i);

    fmpz_mul(res->den, poly1->den, poly2->den);
    _fmpq_poly_set_length(res, len);
    _fmpq_poly_canonicalise(res->coeffs, res->den, len);
}

static void
fmpq_poly_borel_transform(fmpq_poly_t res, const fmpq_poly_t poly)
{
    slong i, len;

    len = fmpq_poly_length(poly);

    if (len <= 2)
    {
        fmpq_poly_set(res, poly);
    }
    else
    {
        fmpz_t c;
        fmpz_init(c);

        fmpz_one(c);
        fmpq_poly_fit_length(res, len);

        for (i = len - 1; i >= 0; i--)
        {
            fmpz_mul(res->coeffs + i, poly->coeffs + i, c);
            if (i > 1)
                fmpz_mul_ui(c, c, i);
        }

        fmpz_mul(fmpq_poly_denref(res), fmpq_poly_denref(poly), c);

        _fmpq_poly_set_length(res, len);
        _fmpq_poly_canonicalise(res->coeffs, res->den, len);

        fmpz_clear(c);
    }
}

static void
fmpq_poly_inv_borel_transform(fmpq_poly_t res, const fmpq_poly_t poly)
{
    slong i, len;

    len = fmpq_poly_length(poly);

    if (len <= 2)
    {
        fmpq_poly_set(res, poly);
    }
    else
    {
        fmpz_t c;
        fmpz_init(c);

        fmpz_one(c);

        fmpq_poly_fit_length(res, len);
        fmpz_set(fmpq_poly_denref(res), fmpq_poly_denref(poly));
        fmpz_set(res->coeffs, poly->coeffs);
        fmpz_set(res->coeffs + 1, poly->coeffs + 1);

        for (i = 2; i < len; i++)
        {
            fmpz_mul_ui(c, c, i);
            fmpz_mul(res->coeffs + i, poly->coeffs + i, c);
        }

        _fmpq_poly_set_length(res, len);
        _fmpq_poly_canonicalise(res->coeffs, res->den, len);

        fmpz_clear(c);
    }
}

void
qqbar_fmpz_poly_composed_op(fmpz_poly_t res, const fmpz_poly_t A, const fmpz_poly_t B, int op)
{
    slong d1, d2, n, i;
    fmpq_poly_t P1, P2, P1rev, P1drev, P2rev, P2drev;

    d1 = fmpz_poly_degree(A);
    d2 = fmpz_poly_degree(B);

    if (d1 <= 0 || d2 <= 0)
    {
        flint_throw(FLINT_ERROR, "composed_op: inputs must not be constants\n");
    }

    n = d1 * d2 + 1;

    fmpq_poly_init(P1);
    fmpq_poly_init(P2);
    fmpq_poly_init(P1rev);
    fmpq_poly_init(P1drev);
    fmpq_poly_init(P2rev);
    fmpq_poly_init(P2drev);

    fmpq_poly_set_fmpz_poly(P1, A);
    fmpq_poly_set_fmpz_poly(P2, B);

    if (op == OP_DIV)
    {
        if (fmpz_is_zero(P2->coeffs))
        {
            flint_throw(FLINT_DIVZERO, "composed_op: division by zero\n");
        }

        fmpq_poly_reverse(P2, P2, d2 + 1);
    }

    if (op == OP_SUB)
        for (i = 1; i <= d2; i += 2)
            fmpz_neg(P2->coeffs + i, P2->coeffs + i);

    fmpq_poly_reverse(P1rev, P1, d1 + 1);
    fmpq_poly_derivative(P1drev, P1);
    fmpq_poly_reverse(P1drev, P1drev, d1);

    fmpq_poly_reverse(P2rev, P2, d2 + 1);
    fmpq_poly_derivative(P2drev, P2);
    fmpq_poly_reverse(P2drev, P2drev, d2);

    fmpq_poly_div_series(P1, P1drev, P1rev, n);
    fmpq_poly_div_series(P2, P2drev, P2rev, n);

    if (op == OP_MUL || op == OP_DIV)
    {
        fmpq_poly_hadamard_product(P1, P1, P2);
        fmpq_poly_shift_right(P1, P1, 1);
        fmpq_poly_neg(P1, P1);
        fmpq_poly_integral(P1, P1);
    }
    else
    {
        fmpq_poly_borel_transform(P1, P1);
        fmpq_poly_borel_transform(P2, P2);
        fmpq_poly_mullow(P1, P1, P2, n);
        fmpq_poly_shift_right(P1, P1, 1);
        fmpq_poly_inv_borel_transform(P1, P1);
        fmpq_poly_neg(P1, P1);
        fmpq_poly_shift_left(P1, P1, 1);
    }

    fmpq_poly_exp_series(P1, P1, n);
    fmpq_poly_reverse(P1, P1, n);

    fmpq_poly_get_numerator(res, P1);

    fmpq_poly_clear(P1);
    fmpq_poly_clear(P2);
    fmpq_poly_clear(P1rev);
    fmpq_poly_clear(P1drev);
    fmpq_poly_clear(P2rev);
    fmpq_poly_clear(P2drev);
}

/*
    Symmetric composed operations. Given A with roots a_1, ..., a_d
    (d >= 2), computes a nonzero integer multiple of

        op = OP_ADD:  prod_{i<j} (z - (a_i + a_j))
        op = OP_MUL:  prod_{i<j} (z - a_i a_j)
        op = OP_SUB:  prod_{i<j} (z - (a_i - a_j)^2)

    each of degree d(d-1)/2. If A is irreducible and x != y are two roots
    of A, then x + y, x y and (x - y)^2 are roots of these polynomials.
    Compared to the general composed operations (degree d^2, with the
    square of the polynomial above as a factor), this halves the degree
    and removes the need for squarefree factorization.
    With power sums p_k = sum_i a_i^k:

        sum_{i<j} (a_i + a_j)^k = (sum_{m} binom(k,m) p_m p_{k-m} - 2^k p_k) / 2,
        sum_{i<j} (a_i a_j)^k = (p_k^2 - p_{2k}) / 2,
        sum_{i<j} (a_i - a_j)^(2k) = (1/2) sum_{m} binom(2k,m) (-1)^m p_m p_{2k-m}.

    The first and third are computed using exponential generating
    functions: E(t)^2 - E(2t) and E(t) E(-t) where E(t) = sum p_k t^k / k!.
*/
void
qqbar_fmpz_poly_symmetric_composed_op(fmpz_poly_t res, const fmpz_poly_t A, int op)
{
    slong d, m, n, len, i;
    fmpq_poly_t P, Prev, Pdrev, E, F;
    fmpq_t c;

    d = fmpz_poly_degree(A);

    if (d <= 1)
    {
        flint_throw(FLINT_ERROR, "symmetric_composed_op: input must have degree >= 2\n");
    }

    m = d * (d - 1) / 2;
    n = m + 1;

    /* number of power sums needed */
    len = (op == OP_ADD) ? n : 2 * n - 1;

    fmpq_poly_init(P);
    fmpq_poly_init(Prev);
    fmpq_poly_init(Pdrev);
    fmpq_poly_init(E);
    fmpq_poly_init(F);
    fmpq_init(c);

    fmpq_poly_set_fmpz_poly(P, A);
    fmpq_poly_reverse(Prev, P, d + 1);
    fmpq_poly_derivative(Pdrev, P);
    fmpq_poly_reverse(Pdrev, Pdrev, d);

    /* power sums p_0, ..., p_{len-1} */
    fmpq_poly_div_series(P, Pdrev, Prev, len);

    if (op == OP_MUL)
    {
        /* E = sum_k (p_k^2 - p_{2k}) / 2 t^k */
        for (i = 0; i < n; i++)
        {
            fmpq_poly_get_coeff_fmpq(c, P, i);
            fmpq_mul(c, c, c);
            fmpq_poly_set_coeff_fmpq(E, i, c);
            fmpq_poly_get_coeff_fmpq(c, P, 2 * i);
            fmpq_poly_set_coeff_fmpq(F, i, c);
        }

        fmpq_poly_sub(E, E, F);
        fmpq_poly_scalar_div_ui(E, E, 2);
    }
    else if (op == OP_ADD)
    {
        /* E = (B(t)^2 - B(2t)) / 2, B(t) = sum_k p_k t^k / k!,
           gives the exponential generating function for the power sums */
        fmpq_poly_borel_transform(P, P);
        fmpq_poly_mullow(E, P, P, n);
        fmpz_set_ui(fmpq_numref(c), 2);
        fmpz_one(fmpq_denref(c));
        fmpq_poly_rescale(F, P, c);
        fmpq_poly_sub(E, E, F);
        fmpq_poly_scalar_div_ui(E, E, 2);
        /* convert to ordinary generating function */
        fmpq_poly_inv_borel_transform(E, E);
    }
    else
    {
        /* B(t) B(-t) gives the exponential generating function for
           sum_{i,j} (a_i - a_j)^k; we want (1/2) of the even terms */
        fmpq_poly_borel_transform(P, P);
        fmpz_set_si(fmpq_numref(c), -1);
        fmpz_one(fmpq_denref(c));
        fmpq_poly_rescale(F, P, c);
        fmpq_poly_mullow(F, P, F, 2 * n - 1);
        fmpq_poly_inv_borel_transform(F, F);

        for (i = 0; i < n; i++)
        {
            fmpq_poly_get_coeff_fmpq(c, F, 2 * i);
            fmpq_poly_set_coeff_fmpq(E, i, c);
        }

        fmpq_poly_scalar_div_ui(E, E, 2);
    }

    /* E = sum_k s_k t^k (s_k = power sums of the roots). Recover the
       polynomial as exp(-sum_{k>=1} s_k t^k / k), reversed. */
    fmpq_poly_shift_right(E, E, 1);
    fmpq_poly_neg(E, E);
    fmpq_poly_integral(E, E);
    fmpq_poly_exp_series(E, E, n);
    fmpq_poly_reverse(E, E, n);

    fmpq_poly_get_numerator(res, E);

    fmpq_poly_clear(P);
    fmpq_poly_clear(Prev);
    fmpq_poly_clear(Pdrev);
    fmpq_poly_clear(E);
    fmpq_poly_clear(F);
    fmpq_clear(c);
}

/* Given a factorization containing the minimal polynomial of w as a
   factor, where w is computed from x and y as follows, identifies the
   correct factor and sets res to w.

    op = 0: w = x + y
    op = 1: w = x - y
    op = 2: w = x * y
    op = 3: w = x / y
    op = 4: w = (x - y)^2
    op = 5: w = (x - y) / (2i)
*/
void
_qqbar_binary_op_select_factor(qqbar_t res, const fmpz_poly_factor_t fac, const qqbar_t x, const qqbar_t y, int op)
{
    slong i, prec, found;
    acb_t z1, z2, w, t;

    acb_init(z1);
    acb_init(z2);
    acb_init(w);
    acb_init(t);

    acb_set(z1, QQBAR_ENCLOSURE(x));
    acb_set(z2, QQBAR_ENCLOSURE(y));

    for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        _qqbar_enclosure_raw(z1, QQBAR_POLY(x), z1, prec);
        _qqbar_enclosure_raw(z2, QQBAR_POLY(y), z2, prec);

        if (op == 0)
            acb_add(w, z1, z2, prec);
        else if (op == 1)
            acb_sub(w, z1, z2, prec);
        else if (op == 2)
            acb_mul(w, z1, z2, prec);
        else if (op == 3)
            acb_div(w, z1, z2, prec);
        else if (op == 4)
        {
            acb_sub(w, z1, z2, prec);
            acb_sqr(w, w, prec);
        }
        else
        {
            acb_sub(w, z1, z2, prec);
            acb_div_onei(w, w);
            acb_mul_2exp_si(w, w, -1);
        }

        /* Look for potential roots -- we want exactly one */
        found = -1;
        for (i = 0; i < fac->num && found != -2; i++)
        {
            arb_fmpz_poly_evaluate_acb(t, fac->p + i, w, prec);
            if (acb_contains_zero(t))
            {
                if (found == -1)
                    found = i;
                else
                    found = -2;
            }
        }

        /* Check if the enclosure is good enough */
        if (found >= 0)
        {
            if (_qqbar_validate_uniqueness(t, fac->p + found, w, 2 * prec))
            {
                fmpz_poly_set(QQBAR_POLY(res), fac->p + found);
                acb_set(QQBAR_ENCLOSURE(res), t);
                break;
            }
        }
    }

    acb_clear(z1);
    acb_clear(z2);
    acb_clear(w);
    acb_clear(t);
}

static void
qqbar_binary_op_without_guess(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)
{
    fmpz_poly_t H;
    fmpz_poly_factor_t fac;

    fmpz_poly_init(H);
    fmpz_poly_factor_init(fac);

    qqbar_fmpz_poly_composed_op(H, QQBAR_POLY(x), QQBAR_POLY(y), op);
    fmpz_poly_factor(fac, H);
    _qqbar_binary_op_select_factor(res, fac, x, y, op);

    fmpz_poly_clear(H);
    fmpz_poly_factor_clear(fac);
}

/* Given T irreducible (primitive, positive leading coefficient), factor
   T(z^2), using a Capelli certificate to avoid refactoring T when possible. */
void
_qqbar_factor_inflate2_irreducible(fmpz_poly_factor_t fac, const fmpz_poly_t T)
{
    fmpz_poly_t U;
    fmpz_poly_init(U);
    fmpz_poly_inflate(U, T, 2);

    if (fmpz_poly_degree(T) >= 1 && !fmpz_is_zero(T->coeffs) &&
        _fmpz_poly_factor_inflation_is_irreducible_capelli(T, 2))
    {
        fmpz_poly_factor_fit_length(fac, 1);
        fmpz_one(&fac->c);
        fmpz_poly_swap(fac->p, U);
        fac->exp[0] = 1;
        fac->num = 1;
    }
    else
    {
        fmpz_poly_factor(fac, U);
    }

    fmpz_poly_clear(U);
}

/* x != y are roots of the same irreducible polynomial.
   Computes x + y (op = 0), x - y (op = 1), x * y (op = 2),
   (x - y)^2 (op = 4), or (x - y) / (2i) (op = 5). */
void
_qqbar_conjugate_pair_op(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)
{
    fmpz_poly_t H;
    fmpz_poly_factor_t fac;

    fmpz_poly_init(H);
    fmpz_poly_factor_init(fac);

    if (op == 0 || op == 2)
    {
        qqbar_fmpz_poly_symmetric_composed_op(H, QQBAR_POLY(x), op);
        fmpz_poly_factor(fac, H);
        _qqbar_binary_op_select_factor(res, fac, x, y, op);
    }
    else
    {
        qqbar_t v;
        qqbar_init(v);

        /* v = (x - y)^2 */
        qqbar_fmpz_poly_symmetric_composed_op(H, QQBAR_POLY(x), OP_SUB);
        fmpz_poly_factor(fac, H);
        _qqbar_binary_op_select_factor(v, fac, x, y, 4);

        if (op == 4)
        {
            qqbar_swap(res, v);
        }
        else
        {
            /* w = x - y is a root of f(z^2), (x - y) / (2i) of f(-4 z^2),
               where f = minpoly(v) */
            if (op == 5)
            {
                slong i;
                for (i = 1; i < QQBAR_POLY(v)->length; i++)
                    fmpz_mul_2exp(QQBAR_COEFFS(v) + i, QQBAR_COEFFS(v) + i, 2 * i);
                for (i = 1; i < QQBAR_POLY(v)->length; i += 2)
                    fmpz_neg(QQBAR_COEFFS(v) + i, QQBAR_COEFFS(v) + i);
                fmpz_poly_primitive_part(QQBAR_POLY(v), QQBAR_POLY(v));
                if (fmpz_sgn(QQBAR_COEFFS(v) + qqbar_degree(v)) < 0)
                    fmpz_poly_neg(QQBAR_POLY(v), QQBAR_POLY(v));
            }

            _qqbar_factor_inflate2_irreducible(fac, QQBAR_POLY(v));
            _qqbar_binary_op_select_factor(res, fac, x, y, op);
        }

        qqbar_clear(v);
    }

    fmpz_poly_clear(H);
    fmpz_poly_factor_clear(fac);
}

void
qqbar_binary_op(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)
{
    slong dx, dy;
    int same_poly;

    dx = qqbar_degree(x);
    dy = qqbar_degree(y);

    same_poly = (dx == dy) && fmpz_poly_equal(QQBAR_POLY(x), QQBAR_POLY(y));

    if (same_poly && qqbar_equal(x, y))
    {
        if (op == 0)
            qqbar_mul_2exp_si(res, x, 1);
        else if (op == 1)
            qqbar_zero(res);
        else if (op == 2)
            qqbar_sqr(res, x);
        else
            qqbar_one(res);
        return;
    }

    /* Specialized algorithms when one operand is quadratic. */
    if (dx == 2 || dy == 2)
    {
        _qqbar_binary_op_quadratic(res, x, y, op);
        return;
    }

    /* Guess and verify rational result; this could be generalized to
       higher degree results. */
    if (dx >= 4 && dy >= 4 && dx == dy)
    {
        qqbar_t t, u;
        acb_t z;
        slong prec;
        int found;

        found = 0;
        prec = QQBAR_DEFAULT_PREC;  /* Could be set dynamically (with recomputation...) */
        qqbar_init(t);
        qqbar_init(u);
        acb_init(z);

        if (op == 0)
            acb_add(z, QQBAR_ENCLOSURE(x), QQBAR_ENCLOSURE(y), prec);
        else if (op == 1)
            acb_sub(z, QQBAR_ENCLOSURE(x), QQBAR_ENCLOSURE(y), prec);
        else if (op == 2)
            acb_mul(z, QQBAR_ENCLOSURE(x), QQBAR_ENCLOSURE(y), prec);
        else if (op == 3)
            acb_div(z, QQBAR_ENCLOSURE(x), QQBAR_ENCLOSURE(y), prec);

        if (qqbar_guess(t, z, 1, prec, 0, prec))
        {
            /* x + y = t  <=>  x = t - y */
            /* x - y = t  <=>  x = t + y */
            /* x * y = t  <=>  x = t / y */
            /* x / y = t  <=>  x = t * y */
            if (op == 0)
                qqbar_sub(u, t, y);
            else if (op == 1)
                qqbar_add(u, t, y);
            else if (op == 2)
                qqbar_div(u, t, y);
            else if (op == 3)
                qqbar_mul(u, t, y);

            if (qqbar_equal(x, u))
            {
                qqbar_swap(res, t);
                found = 1;
            }
        }

        qqbar_clear(t);
        qqbar_clear(u);
        acb_clear(z);

        if (found)
            return;
    }

    /* If one operand lies in the field generated by the other, or both lie
       in a common cyclotomic field, we can compute the result by linear
       algebra in that field instead of factoring a composed polynomial
       of degree dx * dy. */
    if (_qqbar_binary_op_structured(res, x, y, op))
        return;

    /* Distinct conjugates: use a composed polynomial of half the degree. */
    if (same_poly && op != 3)
    {
        _qqbar_conjugate_pair_op(res, x, y, op);
        return;
    }

    qqbar_binary_op_without_guess(res, x, y, op);
}
