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

void
fmpz_mpoly_q_factored_canonicalise(fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_struct * num = fmpz_mpoly_q_factored_numref(x);
    fmpz_mpoly_factor_struct * den = fmpz_mpoly_q_factored_denref(x);
    fmpz_mpoly_factor_t F, G;
    fmpz_mpoly_t t;
    fmpz_t s, e;
    slong i, j;

    if (fmpz_is_zero(den->constant) || fmpz_is_zero(den->constant_den))
        flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored_canonicalise: division by zero\n");

    for (i = 0; i < den->num; i++)
        if (fmpz_sgn(den->exp + i) != 0 && fmpz_mpoly_is_zero(den->poly + i, ctx))
            flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored_canonicalise: division by zero\n");

    if (!fmpz_is_one(den->constant_den))
    {
        fmpz_mpoly_scalar_mul_fmpz(num, num, den->constant_den, ctx);
        fmpz_one(den->constant_den);
    }

    if (fmpz_mpoly_is_zero(num, ctx))
    {
        fmpz_mpoly_factor_one(den, ctx);
        return;
    }

    fmpz_mpoly_factor_init(F, ctx);
    fmpz_mpoly_factor_init(G, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_init(s);
    fmpz_init(e);

    fmpz_set(F->constant, den->constant);

    /* Factor every base of the denominator into irreducibles. Factors with
       negative exponent are moved to the numerator. */
    for (i = 0; i < den->num; i++)
    {
        if (fmpz_is_zero(den->exp + i))
            continue;

        if (fmpz_sgn(den->exp + i) < 0)
        {
            fmpz_neg(e, den->exp + i);
            if (!fmpz_mpoly_pow_fmpz(t, den->poly + i, e, ctx))
                flint_throw(FLINT_ERROR, "fmpz_mpoly_q_factored_canonicalise: exponent too large\n");
            fmpz_mpoly_mul(num, num, t, ctx);
            continue;
        }

        _fmpz_mpoly_q_factored_factor_num(s, G, den->poly + i, ctx);

        if (!fmpz_is_pm1(s) || (fmpz_sgn(s) < 0 && fmpz_is_odd(den->exp + i)))
        {
            if (!fmpz_pow_fmpz(s, s, den->exp + i))
                flint_throw(FLINT_ERROR, "fmpz_mpoly_q_factored_canonicalise: exponent too large\n");
            fmpz_mul(F->constant, F->constant, s);
        }

        for (j = 0; j < G->num; j++)
        {
            fmpz_mul(e, G->exp + j, den->exp + i);
            _fmpz_mpoly_q_factored_append(F, G->poly + j, e, ctx);
        }
    }

    _fmpz_mpoly_q_factored_sort_merge(F, ctx);

    /* Normalise the sign. */
    if (fmpz_sgn(F->constant) < 0)
    {
        fmpz_neg(F->constant, F->constant);
        fmpz_mpoly_neg(num, num, ctx);
    }

    /* Remove common polynomial factors. */
    for (i = 0; i < F->num; i++)
        _fmpz_mpoly_q_factored_divide_out(num, F->exp + i, F->poly + i, ctx);

    _fmpz_mpoly_q_factored_sort_merge(F, ctx);

    /* Remove common content. */
    _fmpz_vec_content2(s, num->coeffs, num->length, F->constant);
    if (!fmpz_is_one(s))
    {
        fmpz_mpoly_scalar_divexact_fmpz(num, num, s, ctx);
        fmpz_divexact(F->constant, F->constant, s);
    }

    fmpz_mpoly_factor_swap(den, F, ctx);

    fmpz_mpoly_factor_clear(F, ctx);
    fmpz_mpoly_factor_clear(G, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_clear(s);
    fmpz_clear(e);
}

int
fmpz_mpoly_q_factored_is_canonical(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
{
    const fmpz_mpoly_struct * num = fmpz_mpoly_q_factored_numref(x);
    const fmpz_mpoly_factor_struct * den = fmpz_mpoly_q_factored_denref(x);
    fmpz_mpoly_factor_t G;
    fmpz_mpoly_t t;
    fmpz_t g;
    slong i;
    int result = 1;

    if (!fmpz_mpoly_is_canonical(num, ctx))
        return 0;

    if (!fmpz_is_one(den->constant_den) || fmpz_sgn(den->constant) <= 0)
        return 0;

    if (fmpz_mpoly_is_zero(num, ctx))
        return fmpz_is_one(den->constant) && den->num == 0;

    fmpz_init(g);
    _fmpz_vec_content2(g, num->coeffs, num->length, den->constant);
    result = fmpz_is_one(g);
    fmpz_clear(g);

    if (!result)
        return 0;

    for (i = 0; i < den->num && result; i++)
    {
        if (!fmpz_mpoly_is_canonical(den->poly + i, ctx) ||
            fmpz_mpoly_is_fmpz(den->poly + i, ctx) ||
            fmpz_sgn(den->poly[i].coeffs) <= 0 ||
            fmpz_sgn(den->exp + i) <= 0)
            result = 0;

        if (i > 0 && fmpz_mpoly_cmp(den->poly + i - 1, den->poly + i, ctx) >= 0)
            result = 0;
    }

    if (!result)
        return 0;

    fmpz_mpoly_factor_init(G, ctx);
    fmpz_mpoly_init(t, ctx);

    for (i = 0; i < den->num && result; i++)
    {
        /* irreducible */
        if (!fmpz_mpoly_factor(G, den->poly + i, ctx) ||
            !fmpz_is_one(G->constant) || G->num != 1 || !fmpz_is_one(G->exp))
            result = 0;

        /* coprime to the numerator */
        if (result && fmpz_mpoly_divides(t, num, den->poly + i, ctx))
            result = 0;
    }

    fmpz_mpoly_factor_clear(G, ctx);
    fmpz_mpoly_clear(t, ctx);

    return result;
}
