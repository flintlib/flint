/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_vec.h"

/* todo */
void ca_rewrite_complex_normal_form(ca_t res, const ca_t x, int deep, ca_ctx_t ctx);

/* todo */
void ca_set_ext(ca_t res, ca_ext_srcptr ext, ca_ctx_t ctx);

/* todo */
ulong qqbar_try_as_cyclotomic(qqbar_t zeta, fmpq_poly_t poly, const qqbar_t x);

void
ca_atan_logarithmic(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u, v;
    acb_t z;
    arb_t one, minus_one;

    acb_init(z);
    arb_init(one);
    arb_init(minus_one);

    ca_init(t, ctx);
    ca_init(u, ctx);
    ca_init(v, ctx);

    ca_i(t, ctx);
    ca_mul(u, x, t, ctx);
    /* v = 1 + i x */
    ca_add_ui(v, u, 1, ctx);
    /* res = 1 - i x */
    ca_sub_ui(res, u, 1, ctx);
    ca_neg(res, res, ctx);

    ca_get_acb(z, x, ctx->options[CA_OPT_LOW_PREC], ctx);
    arb_set_si(one, 1);
    arb_set_si(minus_one, -1);

    if (arb_lt(acb_imagref(z), one))
    {
        /* atan(x) = i/2 log((1-ix)/(1+ix)) */
        ca_div(res, res, v, ctx);
        ca_log(res, res, ctx);
        ca_mul(res, res, t, ctx);
        ca_div_ui(res, res, 2, ctx);
    }
    else if (arb_gt(acb_imagref(z), minus_one))
    {
        /* atan(x) = -i/2 log((1+ix)/(1-ix)) */
        ca_div(res, v, res, ctx);
        ca_log(res, res, ctx);
        ca_mul(res, res, t, ctx);
        ca_div_ui(res, res, 2, ctx);
        ca_neg(res, res, ctx);
    }
    else
    {
        /* atan(x) = i/2 (log(1-ix) - log(1+ix)) */
        ca_log(res, res, ctx);
        ca_log(v, v, ctx);
        ca_sub(res, res, v, ctx);
        ca_mul(res, res, t, ctx);
        ca_div_ui(res, res, 2, ctx);
    }

    ca_clear(t, ctx);
    ca_clear(u, ctx);
    ca_clear(v, ctx);

    acb_clear(z);
    arb_clear(one);
    arb_clear(minus_one);
}

void
ca_sin_cos_exponential(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
{
    ca_t ix, y, t;

    ca_init(ix, ctx);
    ca_init(y, ctx);
    ca_init(t, ctx);

    ca_i(ix, ctx);
    ca_mul(ix, x, ix, ctx);
    ca_exp(y, ix, ctx);
    ca_inv(t, y, ctx);

    if (res2 != NULL)
    {
        ca_add(res2, y, t, ctx);
        ca_div_ui(res2, res2, 2, ctx);
    }

    if (res1 != NULL)
    {
        ca_sub(res1, y, t, ctx);
        ca_div_ui(res1, res1, 2, ctx);
        ca_neg_i(t, ctx);
        ca_mul(res1, res1, t, ctx);
    }

    ca_clear(ix, ctx);
    ca_clear(y, ctx);
    ca_clear(t, ctx);
}

/* todo: Re, Im, Abs, Sgn ... */

void
ca_rewrite_ext_complex_normal_form(ca_t res, ca_ext_ptr ext, int deep, ca_ctx_t ctx)
{
    switch (CA_EXT_HEAD(ext))
    {
        case CA_QQBar:

            if (qqbar_is_i(CA_EXT_QQBAR(ext)))
            {
                ca_set_ext(res, ext, ctx);
            }
            else if (qqbar_is_root_of_unity(NULL, NULL, CA_EXT_QQBAR(ext)))
            {
                ca_set_ext(res, ext, ctx);
            }
            else
            {
                fmpq_poly_t poly;
                qqbar_t zeta;
                ulong q;

                qqbar_init(zeta);
                fmpq_poly_init(poly);

                q = qqbar_try_as_cyclotomic(zeta, poly, CA_EXT_QQBAR(ext));

                if (q != 0)
                {
                    ca_set_qqbar(res, zeta, ctx);
                    ca_fmpq_poly_evaluate(res, poly, res, ctx);
                }
                else
                {
                    ca_set_ext(res, ext, ctx);
                }

                qqbar_clear(zeta);
                fmpq_poly_clear(poly);
            }
            break;

        case CA_Re:
        case CA_Im:
        case CA_Conjugate:
        case CA_Abs:
        case CA_Sign:
        case CA_Arg:
            {
                ca_t t, u;
                truth_t is_zero;

                ca_init(t, ctx);
                ca_init(u, ctx);

                if (deep)
                    ca_rewrite_complex_normal_form(t, CA_EXT_FUNC_ARGS(ext), deep, ctx);
                else
                    ca_set(t, CA_EXT_FUNC_ARGS(ext), ctx);

                switch (CA_EXT_HEAD(ext))
                {
                    case CA_Re:
                        ca_conj_deep(u, t, ctx);
                        ca_add(res, t, u, ctx);
                        ca_div_ui(res, res, 2, ctx);
                        break;

                    case CA_Im:
                        ca_conj_deep(u, t, ctx);
                        ca_sub(res, t, u, ctx);
                        ca_div_ui(res, res, 2, ctx);
                        ca_neg_i(t, ctx);
                        ca_mul(res, res, t, ctx);
                        break;

                    case CA_Conjugate:
                        ca_conj_deep(res, t, ctx);
                        break;

                    case CA_Abs:
                        ca_conj_deep(u, t, ctx);
                        ca_mul(t, t, u, ctx);
                        ca_sqrt(res, t, ctx);
                        break;

                    case CA_Sign:
                        is_zero = ca_check_is_zero(t, ctx);
                        if (is_zero == T_TRUE)
                        {
                            ca_zero(res, ctx);
                        }
                        else if (is_zero == T_FALSE)
                        {
                            ca_conj_deep(u, t, ctx);
                            ca_mul(u, t, u, ctx);
                            ca_sqrt(u, u, ctx);
                            /* todo: division by zero is impossible; optimizations here */
                            ca_div(res, t, u, ctx);
                        }
                        else
                        {
                            ca_sgn(res, t, ctx);
                        }
                        break;

                    case CA_Arg:
                        is_zero = ca_check_is_zero(t, ctx);
                        if (is_zero == T_TRUE)
                        {
                            ca_zero(res, ctx);
                        }
                        else if (is_zero == T_FALSE)
                        {
                            /* todo: better to create two logs? */
                            ca_conj_deep(u, t, ctx);
                            ca_mul(u, t, u, ctx);
                            ca_sqrt(u, u, ctx);
                            ca_div(u, t, u, ctx);
                            ca_log(u, u, ctx);
                            ca_neg_i(t, ctx);
                            ca_mul(res, t, u, ctx);
                        }
                        else
                        {
                            ca_arg(res, t, ctx);
                        }
                        break;

                    default:
                        flint_abort();
                }

                ca_clear(t, ctx);
                ca_clear(u, ctx);
            }
            break;

        /* sin(x) = (exp(ix) - 1/exp(ix))/(2i) */
        /* cos(x) = (exp(ix) + 1/exp(ix))/2 */
        /* tan(x) = 2i/(exp(2ix)+1) */
        case CA_Sin:
        case CA_Cos:
        case CA_Tan:
            {
                ca_t t, u;

                ca_init(t, ctx);
                ca_init(u, ctx);

                if (deep)
                    ca_rewrite_complex_normal_form(t, CA_EXT_FUNC_ARGS(ext), deep, ctx);
                else
                    ca_set(t, CA_EXT_FUNC_ARGS(ext), ctx);

                if (CA_EXT_HEAD(ext) == CA_Sin)
                    ca_sin_cos_exponential(res, NULL, t, ctx);
                else if (CA_EXT_HEAD(ext) == CA_Cos)
                    ca_sin_cos_exponential(NULL, res, t, ctx);
                else
                {
                    ca_sin_cos_exponential(t, u, t, ctx);
                    ca_div(res, t, u, ctx);
                }

                ca_clear(t, ctx);
                ca_clear(u, ctx);
            }
            break;

        case CA_Atan:
            {
                ca_t t;
                ca_init(t, ctx);
                if (deep)
                    ca_rewrite_complex_normal_form(t, CA_EXT_FUNC_ARGS(ext), deep, ctx);
                else
                    ca_set(t, CA_EXT_FUNC_ARGS(ext), ctx);
                ca_atan_logarithmic(res, t, ctx);
                ca_clear(t, ctx);
            }
            break;

        case CA_Sqrt:
        case CA_Exp:
        case CA_Log:
            if (deep)
            {
                ca_rewrite_complex_normal_form(res, CA_EXT_FUNC_ARGS(ext), deep, ctx);
                if (ca_equal_repr(res, CA_EXT_FUNC_ARGS(ext), ctx))
                {
                    ca_set_ext(res, ext, ctx);
                }
                else
                {
                    switch (CA_EXT_HEAD(ext))
                    {
                        case CA_Exp:
                            ca_exp(res, res, ctx);
                            break;
                        case CA_Log:
                            ca_log(res, res, ctx);
                            break;
                        case CA_Sqrt:
                            ca_sqrt(res, res, ctx);
                            break;
                        default:
                            flint_abort();
                    }
                }
            }
            else
            {
                ca_set_ext(res, ext, ctx);
            }
            break;

        default:
            /* todo: deep */
            ca_set_ext(res, ext, ctx);
    }

}

void
ca_rewrite_complex_normal_form(ca_t res, const ca_t x, int deep, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_SIGNED_INF(x))
        {
            ca_sgn(res, x, ctx);
            ca_rewrite_complex_normal_form(res, res, deep, ctx);
            if (!ca_is_unknown(res, ctx))
                res->field |= CA_INF;
        }
        else
        {
            ca_set(res, x, ctx);
        }
    }
    else if (CA_IS_QQ(x, ctx) || CA_IS_QQ_I(x, ctx))
    {
        ca_set(res, x, ctx);
    }
    else
    {
        ca_field_ptr K = CA_FIELD(x, ctx);

        if (CA_FIELD_IS_NF(K))
        {
            if (qqbar_is_root_of_unity(NULL, NULL, CA_FIELD_NF_QQBAR(K)))
            {
                ca_set(res, x, ctx);
            }
            else
            {
                fmpq_poly_t poly, xpoly;
                qqbar_t zeta;
                ulong q;

                qqbar_init(zeta);
                fmpq_poly_init(poly);

                q = qqbar_try_as_cyclotomic(zeta, poly, CA_FIELD_NF_QQBAR(K));

                if (q != 0)
                {
                    fmpq_poly_init(xpoly);
                    nf_elem_get_fmpq_poly(xpoly, CA_NF_ELEM(x), CA_FIELD_NF(K));
                    ca_set_qqbar(res, zeta, ctx);
                    ca_fmpq_poly_evaluate(res, poly, res, ctx);
                    ca_fmpq_poly_evaluate(res, xpoly, res, ctx);
                    fmpq_poly_clear(xpoly);
                }
                else
                {
                    ca_set(res, x, ctx);
                }

                qqbar_clear(zeta);
                fmpq_poly_clear(poly);
            }
        }
        else if (0)
        {
            ca_set(res, x, ctx);
        }
        else
        {
            slong nvars;
            int * used;
            slong i;
            ca_ptr cext;

            nvars = CA_FIELD_LENGTH(K);
            used = flint_calloc(nvars, sizeof(int));
            cext = _ca_vec_init(nvars, ctx);

            fmpz_mpoly_q_used_vars(used, CA_MPOLY_Q(x), CA_FIELD_MCTX(K, ctx));

            for (i = 0; i < nvars; i++)
            {
                if (used[i])
                {
                    ca_rewrite_ext_complex_normal_form(cext + i, CA_FIELD_EXT_ELEM(K, i), deep, ctx);
                }
            }

            ca_fmpz_mpoly_q_evaluate_zero_impossible(res, CA_MPOLY_Q(x), cext, CA_FIELD_MCTX(K, ctx), ctx);
            _ca_vec_clear(cext, nvars, ctx);
            flint_free(used);
        }
    }
}
