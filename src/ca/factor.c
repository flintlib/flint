/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_factor.h"
#include "fmpz_mpoly_factor.h"
#include "ca.h"

#define HAVE_MPOLY_FAC 1

void
_ca_factor_fmpz(ca_factor_t res, const fmpz_t x, int inv, ulong flags, ca_ctx_t ctx)
{
    slong i;
    fmpz_factor_t fac;
    ca_t b, e;

    if (fmpz_is_one(x))
        return;

    fmpz_factor_init(fac);

    if (flags & CA_FACTOR_ZZ_FULL)
    {
        fmpz_factor(fac, x);
    }
    else if (flags & CA_FACTOR_ZZ_SMOOTH)
    {
        slong smooth_limit = ctx->options[CA_OPT_SMOOTH_LIMIT];

        fmpz_factor_smooth(fac, x, smooth_limit, -1); /* -1  =>  no primality test */
    }
    else
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    ca_init(b, ctx);
    ca_init(e, ctx);

    if (fac->sign != 1)
    {
        ca_set_si(b, fac->sign, ctx);
        ca_one(e, ctx);
        ca_factor_insert(res, b, e, ctx);
    }

    for (i = 0; i < fac->num; i++)
    {
        ca_set_fmpz(b, fac->p + i, ctx);
        ca_set_si(e, inv ? -(slong)(fac->exp[i]) : (fac->exp[i]), ctx);
        ca_factor_insert(res, b, e, ctx);
    }

    fmpz_factor_clear(fac);
    ca_clear(b, ctx);
    ca_clear(e, ctx);
}

void
_ca_factor_fmpq(ca_factor_t res, const fmpq_t x, ulong flags, ca_ctx_t ctx)
{
    if (flags & (CA_FACTOR_ZZ_SMOOTH | CA_FACTOR_ZZ_FULL))
    {
        _ca_factor_fmpz(res, fmpq_numref(x), 0, flags, ctx);
        _ca_factor_fmpz(res, fmpq_denref(x), 1, flags, ctx);
    }
    else if (!fmpq_is_one(x))
    {
        ca_t b, e;

        ca_init(b, ctx);
        ca_init(e, ctx);

        ca_set_fmpq(b, x, ctx);
        ca_one(e, ctx);

        ca_factor_insert(res, b, e, ctx);

        ca_clear(b, ctx);
        ca_clear(e, ctx);
    }
}

static int
_ca_fmpz_mpoly_factor(fmpz_mpoly_factor_t fac, const fmpz_mpoly_t poly, int full, const fmpz_mpoly_ctx_t ctx)
{
    if (full)
        return fmpz_mpoly_factor(fac, poly, ctx);
    else
        return fmpz_mpoly_factor_squarefree(fac, poly, ctx);
}

void
ca_factor(ca_factor_t res, const ca_t x, ulong flags, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        flint_throw(FLINT_ERROR, "ca_factor: expected a non-special value\n");
    }

    ca_factor_one(res, ctx);

    if (CA_IS_QQ(x, ctx))
    {
        _ca_factor_fmpq(res, CA_FMPQ(x), flags, ctx);
        return;
    }

    /* todo: factoring in number fields */
    if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        ca_t e;
        ca_init(e, ctx);
        ca_one(e, ctx);
        ca_factor_insert(res, x, e, ctx);
        ca_clear(e, ctx);
        return;
    }

    {
        if (flags & (CA_FACTOR_POLY_CONTENT | CA_FACTOR_POLY_SQF | CA_FACTOR_POLY_FULL))
        {
            const fmpz_mpoly_ctx_struct * mctx;
            ca_t b, e;
            fmpq_t content;

            mctx = CA_FIELD_MCTX(CA_FIELD(x, ctx), ctx);

            fmpq_init(content);
            ca_init(b, ctx);
            ca_init(e, ctx);

#if HAVE_MPOLY_FAC
            if (flags & (CA_FACTOR_POLY_SQF | CA_FACTOR_POLY_FULL))
            {
                int full;
                slong i;
                fmpz_mpoly_factor_t mfac;

                full = (flags & CA_FACTOR_POLY_FULL) ? 1 : 0;
                fmpz_mpoly_factor_init(mfac, mctx);

                if (!_ca_fmpz_mpoly_factor(mfac, fmpz_mpoly_q_numref(CA_MPOLY_Q(x)), full, mctx))
                {
                    flint_throw(FLINT_ERROR, "ca_factor: unable to factor numerator\n");
                }

                for (i = 0; i < mfac->num; i++)
                {
                    ca_set_fmpz(e, mfac->exp + i, ctx);
                    _ca_make_field_element(b, CA_FIELD(x, ctx), ctx);
                    fmpz_mpoly_swap(fmpz_mpoly_q_numref(CA_MPOLY_Q(b)), mfac->poly + i, mctx);
                    fmpz_mpoly_one(fmpz_mpoly_q_denref(CA_MPOLY_Q(b)), mctx);
                    ca_factor_insert(res, b, e, ctx);
                }

                fmpz_set(fmpq_numref(content), mfac->constant);

                fmpz_mpoly_factor_clear(mfac, mctx);
                fmpz_mpoly_factor_init(mfac, mctx);

                if (!_ca_fmpz_mpoly_factor(mfac, fmpz_mpoly_q_denref(CA_MPOLY_Q(x)), full, mctx))
                {
                    flint_throw(FLINT_ERROR, "ca_factor: unable to factor denominator\n");
                }

                for (i = 0; i < mfac->num; i++)
                {
                    ca_set_fmpz(e, mfac->exp + i, ctx);
                    ca_neg(e, e, ctx);
                    _ca_make_field_element(b, CA_FIELD(x, ctx), ctx);
                    fmpz_mpoly_swap(fmpz_mpoly_q_numref(CA_MPOLY_Q(b)), mfac->poly + i, mctx);
                    fmpz_mpoly_one(fmpz_mpoly_q_denref(CA_MPOLY_Q(b)), mctx);
                    ca_factor_insert(res, b, e, ctx);
                }

                fmpz_set(fmpq_denref(content), mfac->constant);

                fmpz_mpoly_factor_clear(mfac, mctx);
            }
            else
#endif
            {
                fmpz_mpoly_q_content(content, CA_MPOLY_Q(x), mctx);
                ca_div_fmpq(b, x, content, ctx);
                ca_one(e, ctx);
                ca_factor_insert(res, b, e, ctx);
            }

            if (fmpz_sgn(fmpq_denref(content)) < 0)
            {
                fmpz_neg(fmpq_numref(content), fmpq_numref(content));
                fmpz_neg(fmpq_denref(content), fmpq_denref(content));
            }

            _ca_factor_fmpq(res, content, flags, ctx);

            ca_clear(b, ctx);
            ca_clear(e, ctx);

            fmpq_clear(content);
        }
        else
        {
            ca_t e;
            ca_init(e, ctx);
            ca_one(e, ctx);
            ca_factor_insert(res, x, e, ctx);
            ca_clear(e, ctx);
        }
    }
}

