/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

truth_t
_ca_check_is_zero_qqbar(const ca_t x, ca_ctx_t ctx)
{
    qqbar_t t;
    truth_t res;
    qqbar_init(t);

    if (ca_get_qqbar(t, x, ctx))
        res = qqbar_is_zero(t) ? T_TRUE : T_FALSE;
    else
        res = T_UNKNOWN;

    qqbar_clear(t);
    return res;
}

void
ca_rewrite_complex_normal_form(ca_t res, const ca_t x, int deep, ca_ctx_t ctx);

truth_t
ca_is_zero_check_fast(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }

    if (CA_IS_QQ(x, ctx))
    {
        if (fmpq_is_zero(CA_FMPQ(x)))
            return T_TRUE;
        else
            return T_FALSE;
    }

    if (CA_IS_QQ_I(x, ctx))
    {
        const fmpz *n;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));

        if (fmpz_is_zero(n) && fmpz_is_zero(n + 1))
            return T_TRUE;

        return T_FALSE;
    }

    if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        if (nf_elem_is_zero(CA_NF_ELEM(x), CA_FIELD_NF(CA_FIELD(x, ctx))))
            return T_TRUE;
        else
            return T_FALSE;
    }

    return T_UNKNOWN;
}

int
_ca_generic_has_nontrivial_denominator(const ca_t x, ca_ctx_t ctx)
{
    return !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(CA_MPOLY_Q(x)), CA_FIELD_MCTX(CA_FIELD(x, ctx), ctx));
}

truth_t
ca_check_is_zero_no_factoring(const ca_t x, ca_ctx_t ctx)
{
    acb_t v;
    truth_t res;
    slong prec, prec_limit;

    res = ca_is_zero_check_fast(x, ctx);

    if (res != T_UNKNOWN || CA_IS_SPECIAL(x))
        return res;

    /* The denominator is irrelevant. */
    /* Todo: we should probably also factor out monomials. */
    if (_ca_generic_has_nontrivial_denominator(x, ctx))
    {
        ca_t t;
        ca_init(t, ctx);
        ca_set(t, x, ctx);
        /* Todo: could also remove content */
        fmpz_mpoly_one(fmpz_mpoly_q_denref(CA_MPOLY_Q(t)), CA_FIELD_MCTX(CA_FIELD(t, ctx), ctx));
        res = ca_check_is_zero_no_factoring(t, ctx);
        ca_clear(t, ctx);
        return res;
    }

    acb_init(v);

    prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
    prec_limit = FLINT_MAX(prec_limit, 64);

    for (prec = 64; (prec <= prec_limit) && (res == T_UNKNOWN); prec *= 2)
    {
        ca_get_acb_raw(v, x, prec, ctx);

        if (!acb_contains_zero(v))
        {
            res = T_FALSE;
            break;
        }

        /* try qqbar computation */
        /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
        if (prec == 64)
        {
            res = _ca_check_is_zero_qqbar(x, ctx);
        }
    }

    acb_clear(v);

    if (res == T_UNKNOWN)
    {
        ca_t tmp;
        ca_init(tmp, ctx);
        ca_rewrite_complex_normal_form(tmp, x, 1, ctx);
        res = ca_is_zero_check_fast(tmp, ctx);

        if (ctx->options[CA_OPT_VERBOSE])
        {
            flint_printf("is_zero: complex_normal form:\n");
            ca_print(x, ctx); flint_printf("\n");
            ca_print(tmp, ctx); flint_printf("\n");
            truth_print(res); flint_printf("\n");
        }

        ca_clear(tmp, ctx);
    }

    return res;
}

truth_t
ca_check_is_zero(const ca_t x, ca_ctx_t ctx)
{
    truth_t res;

    res = ca_check_is_zero_no_factoring(x, ctx);

    if (res == T_UNKNOWN && !CA_IS_SPECIAL(x))
    {
        ca_factor_t fac;
        ca_t t;
        truth_t factor_res;
        slong i, nontrivial_factors;

        /* the zero test will surely have succeeded over a number field */
        if (!CA_FIELD_IS_GENERIC(CA_FIELD(x, ctx)))
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        /* extract numerator */
        ca_init(t, ctx);
        ca_set(t, x, ctx);
        fmpz_mpoly_one(fmpz_mpoly_q_denref(CA_MPOLY_Q(t)), CA_FIELD_MCTX(CA_FIELD(t, ctx), ctx));

        ca_factor_init(fac, ctx);
        ca_factor(fac, t, CA_FACTOR_ZZ_NONE | CA_FACTOR_POLY_FULL, ctx);

        nontrivial_factors = 0;
        for (i = 0; i < fac->length; i++)
            nontrivial_factors += !CA_IS_QQ(fac->base + i, ctx);

        if (nontrivial_factors >= 2)
        {
            for (i = 0; i < fac->length; i++)
            {
                factor_res = ca_check_is_zero_no_factoring(fac->base + i, ctx);

                if (factor_res == T_TRUE)
                {
                    if (ctx->options[CA_OPT_VERBOSE])
                    {
                        flint_printf("is_zero: factoring:\n");
                        ca_print(x, ctx); flint_printf("\n");
                        ca_print(fac->base + i, ctx); flint_printf("\n");
                        truth_print(res); flint_printf("\n");
                    }

                    res = T_TRUE;
                    break;
                }
            }
        }

        ca_clear(t, ctx);
        ca_factor_clear(fac, ctx);
    }

    return res;
}
