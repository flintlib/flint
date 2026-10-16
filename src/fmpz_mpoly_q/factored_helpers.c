/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"

/* Append p^e to f. */
void
_fmpz_mpoly_q_factored_append(fmpz_mpoly_factor_t f, const fmpz_mpoly_t p,
    const fmpz_t e, const fmpz_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mpoly_set(f->poly + i, p, ctx);
    fmpz_set(f->exp + i, e);
    f->num = i + 1;
}

/* Cheap necessary condition for p | a in Z[x]: the leading and trailing
   coefficients of p must divide those of a. */
static int
_maybe_divisible(const fmpz_mpoly_t a, const fmpz_mpoly_t p)
{
    if (a->length == 0)
        return 1;

    if (!fmpz_is_pm1(p->coeffs) && !fmpz_divisible(a->coeffs, p->coeffs))
        return 0;

    if (!fmpz_is_pm1(p->coeffs + p->length - 1) &&
        !fmpz_divisible(a->coeffs + a->length - 1, p->coeffs + p->length - 1))
        return 0;

    return 1;
}

/* While e > 0 and p | a, set a = a / p and decrement e. Assumes that
   a is nonzero and that p is not constant. */
void
_fmpz_mpoly_q_factored_divide_out(fmpz_mpoly_t a, fmpz_t e,
    const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t q;

    if (fmpz_sgn(e) <= 0 || fmpz_mpoly_is_zero(a, ctx))
        return;

    if (!_maybe_divisible(a, p))
        return;

    fmpz_mpoly_init(q, ctx);

    while (fmpz_sgn(e) > 0 && _maybe_divisible(a, p) &&
                fmpz_mpoly_divides(q, a, p, ctx))
    {
        fmpz_mpoly_swap(a, q, ctx);
        fmpz_sub_ui(e, e, 1);
    }

    fmpz_mpoly_clear(q, ctx);
}

/* Sort the bases of f according to fmpz_mpoly_cmp, combine equal bases
   and remove bases with exponent zero. The constant is not touched. */
void
_fmpz_mpoly_q_factored_sort_merge(fmpz_mpoly_factor_t f,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, n = f->num;

    /* insertion sort; the number of factors is typically small */
    for (i = 1; i < n; i++)
    {
        for (j = i; j > 0 && fmpz_mpoly_cmp(f->poly + j - 1, f->poly + j, ctx) > 0; j--)
        {
            fmpz_mpoly_swap(f->poly + j - 1, f->poly + j, ctx);
            fmpz_swap(f->exp + j - 1, f->exp + j);
        }
    }

    /* merge equal bases and remove zero exponents */
    j = 0;
    for (i = 0; i < n; i++)
    {
        if (j > 0 && fmpz_mpoly_equal(f->poly + j - 1, f->poly + i, ctx))
        {
            fmpz_add(f->exp + j - 1, f->exp + j - 1, f->exp + i);
        }
        else
        {
            if (j > 0 && fmpz_is_zero(f->exp + j - 1))
                j--;

            if (i != j)
            {
                fmpz_mpoly_swap(f->poly + j, f->poly + i, ctx);
                fmpz_swap(f->exp + j, f->exp + i);
            }
            j++;
        }
    }

    if (j > 0 && fmpz_is_zero(f->exp + j - 1))
        j--;

    f->num = j;
}

/* Factor the nonzero polynomial a as s * prod f_i^e_i where s is an
   integer and the f_i are irreducible and normalised as in the denominator
   of a canonical fmpz_mpoly_q_factored_t. The constant of f is set to one. */
void
_fmpz_mpoly_q_factored_factor_num(fmpz_t s, fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t a, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (fmpz_mpoly_is_zero(a, ctx))
        flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored: division by zero\n");

    if (fmpz_mpoly_is_fmpz(a, ctx))
    {
        fmpz_set(s, a->coeffs);
        fmpz_mpoly_factor_one(f, ctx);
        return;
    }

    if (!fmpz_mpoly_factor(f, a, ctx))
        flint_throw(FLINT_ERROR, "fmpz_mpoly_q_factored: factorization failed\n");

    fmpz_set(s, f->constant);
    fmpz_one(f->constant);
    fmpz_one(f->constant_den);

    for (i = 0; i < f->num; i++)
    {
        if (fmpz_sgn(f->poly[i].coeffs) < 0)
        {
            fmpz_mpoly_neg(f->poly + i, f->poly + i, ctx);
            if (fmpz_is_odd(f->exp + i))
                fmpz_neg(s, s);
        }
    }

    _fmpz_mpoly_q_factored_sort_merge(f, ctx);
}
