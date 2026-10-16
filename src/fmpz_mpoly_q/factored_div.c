/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"

/* 1 / (xn / xd) = xd / xn. Only the numerator needs to be factored;
   the result is automatically canonical. */
void
fmpz_mpoly_q_factored_inv(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_factor_t F;
    fmpz_mpoly_t N;
    fmpz_t s;

    if (fmpz_mpoly_q_factored_is_zero(x, ctx))
        flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored_inv: division by zero\n");

    fmpz_mpoly_factor_init(F, ctx);
    fmpz_mpoly_init(N, ctx);
    fmpz_init(s);

    _fmpz_mpoly_q_factored_factor_num(s, F, fmpz_mpoly_q_factored_numref(x), ctx);

    if (!fmpz_mpoly_factor_expand(N, fmpz_mpoly_q_factored_denref(x), ctx))
        flint_throw(FLINT_ERROR, "fmpz_mpoly_q_factored_inv: expansion failed\n");

    if (fmpz_sgn(s) < 0)
    {
        fmpz_neg(s, s);
        fmpz_mpoly_neg(N, N, ctx);
    }

    fmpz_swap(F->constant, s);

    fmpz_mpoly_swap(fmpz_mpoly_q_factored_numref(res), N, ctx);
    fmpz_mpoly_factor_swap(fmpz_mpoly_q_factored_denref(res), F, ctx);

    fmpz_mpoly_factor_clear(F, ctx);
    fmpz_mpoly_clear(N, ctx);
    fmpz_clear(s);
}

/*
    (xn / xd) / (yn / yd) = (xn yd) / (xd yn).

    Only yn needs to be factored, as s * Q with Q a product of irreducible
    factors q. Since gcd(xn, xd) = gcd(yn, yd) = 1, the possible
    cancellations are:

      * a factor q of Q not in xd may divide xn (trial division),
      * a factor of yd may coincide with a factor of xd (subtract exponents),

    and the integer content of xn may cancel against s while the constants
    of xd and yd may cancel against each other. Note that a factor q of Q
    cannot occur in yd.
*/
void
fmpz_mpoly_q_factored_div(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)
{
    const fmpz_mpoly_factor_struct * xd = fmpz_mpoly_q_factored_denref(x);
    const fmpz_mpoly_factor_struct * yd = fmpz_mpoly_q_factored_denref(y);
    fmpz_mpoly_factor_t Q, A, D;
    fmpz_mpoly_t xn, M, t;
    int M_is_one;
    fmpz_t s, e, a, b;
    slong i, j;
    int c;

    if (fmpz_mpoly_q_factored_is_zero(y, ctx))
        flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored_div: division by zero\n");

    if (fmpz_mpoly_q_factored_is_zero(x, ctx))
    {
        fmpz_mpoly_q_factored_zero(res, ctx);
        return;
    }

    fmpz_mpoly_factor_init(Q, ctx);
    fmpz_mpoly_factor_init(A, ctx);
    fmpz_mpoly_factor_init(D, ctx);
    fmpz_mpoly_init(xn, ctx);
    fmpz_mpoly_init(M, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_init(s);
    fmpz_init(e);
    fmpz_init(a);
    fmpz_init(b);

    _fmpz_mpoly_q_factored_factor_num(s, Q, fmpz_mpoly_q_factored_numref(y), ctx);
    fmpz_mpoly_set(xn, fmpz_mpoly_q_factored_numref(x), ctx);

    /* A = xd * Q, removing factors of Q from xn */
    fmpz_mpoly_factor_fit_length(A, xd->num + Q->num, ctx);
    i = j = 0;
    while (i < xd->num || j < Q->num)
    {
        if (j >= Q->num)
            c = -1;
        else if (i >= xd->num)
            c = 1;
        else
            c = fmpz_mpoly_cmp(xd->poly + i, Q->poly + j, ctx);

        if (c < 0)
        {
            _fmpz_mpoly_q_factored_append(A, xd->poly + i, xd->exp + i, ctx);
            i++;
        }
        else if (c > 0)
        {
            fmpz_set(e, Q->exp + j);
            _fmpz_mpoly_q_factored_divide_out(xn, e, Q->poly + j, ctx);
            if (!fmpz_is_zero(e))
                _fmpz_mpoly_q_factored_append(A, Q->poly + j, e, ctx);
            j++;
        }
        else
        {
            fmpz_add(e, xd->exp + i, Q->exp + j);
            _fmpz_mpoly_q_factored_append(A, xd->poly + i, e, ctx);
            i++;
            j++;
        }
    }

    /* D = A / yd, moving factors with net negative exponent to M */
    fmpz_mpoly_factor_fit_length(D, A->num, ctx);
    M_is_one = 1;
    i = j = 0;
    while (i < A->num || j < yd->num)
    {
        if (j >= yd->num)
            c = -1;
        else if (i >= A->num)
            c = 1;
        else
            c = fmpz_mpoly_cmp(A->poly + i, yd->poly + j, ctx);

        if (c < 0)
        {
            _fmpz_mpoly_q_factored_append(D, A->poly + i, A->exp + i, ctx);
            i++;
            continue;
        }

        if (c > 0)
        {
            fmpz_set(e, yd->exp + j);
        }
        else
        {
            fmpz_sub(e, A->exp + i, yd->exp + j);

            if (fmpz_sgn(e) > 0)
                _fmpz_mpoly_q_factored_append(D, A->poly + i, e, ctx);

            fmpz_neg(e, e);
            i++;
        }

        if (fmpz_sgn(e) > 0)
        {
            if (!fmpz_mpoly_pow_fmpz(t, yd->poly + j, e, ctx))
                flint_throw(FLINT_ERROR, "fmpz_mpoly_q_factored_div: exponent too large\n");

            if (M_is_one)
            {
                fmpz_mpoly_swap(M, t, ctx);
                M_is_one = 0;
            }
            else
            {
                fmpz_mpoly_mul(M, M, t, ctx);
            }
        }

        j++;
    }

    /* Constants: numerator gets sgn(s) * yd->constant,
       denominator gets xd->constant * |s|. */
    if (fmpz_sgn(s) < 0)
    {
        fmpz_neg(s, s);
        fmpz_mpoly_neg(xn, xn, ctx);
    }

    fmpz_gcd(a, xd->constant, yd->constant);
    _fmpz_vec_content2(b, xn->coeffs, xn->length, s);

    if (!fmpz_is_one(b))
    {
        fmpz_mpoly_scalar_divexact_fmpz(xn, xn, b, ctx);
        fmpz_divexact(s, s, b);
    }

    fmpz_divexact(D->constant, xd->constant, a);
    fmpz_mul(D->constant, D->constant, s);

    if (!fmpz_equal(a, yd->constant))
    {
        fmpz_divexact(e, yd->constant, a);
        fmpz_mpoly_scalar_mul_fmpz(xn, xn, e, ctx);
    }

    if (M_is_one)
        fmpz_mpoly_swap(fmpz_mpoly_q_factored_numref(res), xn, ctx);
    else
        fmpz_mpoly_mul(fmpz_mpoly_q_factored_numref(res), xn, M, ctx);

    fmpz_mpoly_factor_swap(fmpz_mpoly_q_factored_denref(res), D, ctx);

    fmpz_mpoly_factor_clear(Q, ctx);
    fmpz_mpoly_factor_clear(A, ctx);
    fmpz_mpoly_factor_clear(D, ctx);
    fmpz_mpoly_clear(xn, ctx);
    fmpz_mpoly_clear(M, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_clear(s);
    fmpz_clear(e);
    fmpz_clear(a);
    fmpz_clear(b);
}
