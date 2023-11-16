/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_sqrt_factor(ca_t res, const ca_t x, ulong flags, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_SIGNED_INF(x))
        {
            ca_sgn(res, x, ctx);
            ca_sqrt_factor(res, res, flags, ctx);
            if (!ca_is_unknown(res, ctx))
                res->field |= CA_INF;
        }
        else
        {
            ca_set(res, x, ctx);
        }
    }
    else
    {
        slong i;
        ca_factor_t fac;
        ca_t A, B, t;

        if (CA_IS_QQ(x, ctx))
        {
            qqbar_t t;
            qqbar_init(t);
            qqbar_fmpq_root_ui(t, CA_FMPQ(x), 2);
            ca_set_qqbar(res, t, ctx);
            qqbar_clear(t);
            return;
        }

        ca_factor_init(fac, ctx);
        ca_init(A, ctx);
        ca_init(B, ctx);
        ca_init(t, ctx);

        ca_factor(fac, x, flags, ctx);
        ca_one(A, ctx);
        ca_one(B, ctx);

        for (i = 0; i < fac->length; i++)
        {
            if (CA_IS_QQ(fac->exp + i, ctx) && fmpz_is_one(CA_FMPQ_DENREF(fac->exp + i)))
            {
                if (!fmpz_is_zero(CA_FMPQ_DENREF(fac->exp + i)))
                {
                    fmpz_t e;
                    ca_ext_ptr ext;

                    ext = ca_is_gen_as_ext(fac->base + i, ctx);

                    if (ext != NULL && CA_EXT_HEAD(ext) == CA_Exp)
                    {
                        /* sqrt(exp(a)^n) = +/- exp(a*n/2) */
                        ca_mul_fmpz(t, CA_EXT_FUNC_ARGS(ext), CA_FMPQ_NUMREF(fac->exp + i), ctx);
                        ca_div_ui(t, t, 2, ctx);
                        ca_exp(t, t, ctx);
                        ca_mul(A, A, t, ctx);
                        continue;
                    }

                    /* todo: assert a != 0? */
                    if (ext != NULL && CA_EXT_HEAD(ext) == CA_Sqrt)
                    {
                        /* sqrt(sqrt(a)^n) = +/- a^(n/4) */
                        ca_set_fmpz(t, CA_FMPQ_NUMREF(fac->exp + i), ctx);
                        ca_div_ui(t, t, 4, ctx);
                        ca_pow(t, CA_EXT_FUNC_ARGS(ext), t, ctx);
                        ca_mul(A, A, t, ctx);
                        continue;
                    }

                    /* todo: assert a != 0? */
                    if (ext != NULL && CA_EXT_HEAD(ext) == CA_Pow)
                    {
                        /* sqrt((a^b)^n) = +/- a^(n*b/2) */
                        ca_mul_fmpz(t, CA_EXT_FUNC_ARGS(ext) + 1, CA_FMPQ_NUMREF(fac->exp + i), ctx);
                        ca_div_ui(t, t, 2, ctx);
                        ca_pow(t, CA_EXT_FUNC_ARGS(ext), t, ctx);
                        ca_mul(A, A, t, ctx);
                        continue;
                    }

                    fmpz_init(e);

                    if (fmpz_is_odd(CA_FMPQ_NUMREF(fac->exp + i)))
                        ca_mul(B, B, fac->base + i, ctx);

                    fmpz_fdiv_q_2exp(e, CA_FMPQ_NUMREF(fac->exp + i), 1);

                    ca_pow_fmpz(t, fac->base + i, e, ctx);
                    ca_mul(A, A, t, ctx);

                    fmpz_clear(e);
                }
            }
            else
            {
                ca_pow(t, fac->base + i, fac->exp + i, ctx);
                ca_mul(B, B, t, ctx);
            }
        }

/*
        printf("factors:\n");
        ca_factor_print(fac, ctx); printf("\n");
        ca_print(A, ctx); printf("\n");
        ca_print(B, ctx); printf("\n");
*/

        ca_sqrt_nofactor(B, B, ctx);
        ca_mul(A, A, B, ctx);

        /* check sign (todo: much improvable) */
        {
            acb_t sA, sA2, sB;
            int success;
            slong prec, prec_limit, low_prec;

            low_prec = ctx->options[CA_OPT_LOW_PREC];
            prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
            prec_limit = FLINT_MAX(prec_limit, low_prec);

            ca_sqrt_inert(B, x, ctx);

            acb_init(sA);
            acb_init(sA2);
            acb_init(sB);

            success = 0;

            for (prec = low_prec; prec <= prec_limit; prec *= 2)
            {
                ca_get_acb_raw(sA, A, prec, ctx);
                ca_get_acb_raw(sB, B, prec, ctx);
                acb_neg(sA2, sA);

                if (acb_overlaps(sA, sB) && !acb_overlaps(sA2, sB))
                {
                    ca_set(res, A, ctx);
                    success = 1;
                    break;
                }

                if (acb_overlaps(sA2, sB) && !acb_overlaps(sA, sB))
                {
                    ca_neg(res, A, ctx);
                    success = 1;
                    break;
                }
            }

            if (!success)
            {
                if (ca_check_is_zero(A, ctx) == T_TRUE)
                    ca_zero(res, ctx);
                else
                    ca_set(res, B, ctx);
            }

            acb_clear(sA);
            acb_clear(sA2);
            acb_clear(sB);
        }

        ca_factor_clear(fac, ctx);
        ca_clear(A, ctx);
        ca_clear(B, ctx);
        ca_clear(t, ctx);
    }
}

