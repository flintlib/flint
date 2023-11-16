/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"

int
ca_get_qqbar(qqbar_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        return 0;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        qqbar_set_fmpq(res, CA_FMPQ(x));
        return 1;
    }
    else
    {
        if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
        {
            const fmpz * num;
            const fmpz * den;
            slong len;

            if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_LINEAR)
            {
                num = (fmpz *) LNF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = LNF_ELEM_DENREF(CA_NF_ELEM(x));
                len = 1;
            }
            else if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_QUADRATIC)
            {
                num = (fmpz *) QNF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = QNF_ELEM_DENREF(CA_NF_ELEM(x));
                len = 2;
            }
            else
            {
                num = (fmpz *) NF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = NF_ELEM_DENREF(CA_NF_ELEM(x));
                len = NF_ELEM(CA_NF_ELEM(x))->length;
            }

            _qqbar_evaluate_fmpq_poly(res, num, den, len, CA_FIELD_NF_QQBAR(CA_FIELD(x, ctx)));

            return 1;
        }
        else
        {
            slong i, len, deg_limit, bits_limit;
            qqbar_ptr xs;
            qqbar_t y, zero;
            int success;
            int * init_mask, * used;

            if (!ca_can_evaluate_qqbar(x, ctx))
                return 0;

            deg_limit = ctx->options[CA_OPT_QQBAR_DEG_LIMIT];
            bits_limit = 10 * ctx->options[CA_OPT_PREC_LIMIT]; /* xxx */

            len = CA_FIELD_LENGTH(CA_FIELD(x, ctx));

            success = 0;
            xs = (qqbar_struct *) flint_malloc(sizeof(qqbar_struct) * len);
            init_mask = (int *) flint_calloc(sizeof(int), len);
            used = (int *) flint_calloc(sizeof(int), len);
            qqbar_init(y);
            qqbar_init(zero);

            fmpz_mpoly_q_used_vars(used, CA_MPOLY_Q(x), CA_FIELD_MCTX(CA_FIELD(x, ctx), ctx));

            /* todo: allow non-qqbar extension elements to cache a qqbar value */
            for (i = 0; i < len; i++)
            {
                if (!used[i])
                {
                    xs[i] = *zero;
                    continue;
                }

                if (CA_EXT_IS_QQBAR(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)))
                {
                    xs[i] = *CA_EXT_QQBAR(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i));
                }
                else if (CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)->data.func_data.qqbar != NULL)
                {
                    xs[i] = *(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)->data.func_data.qqbar);
                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Sqrt)
                {
                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    if (!ca_get_qqbar(xs + i, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)), ctx))
                        goto cleanup;

                    /* todo: maybe do a x2 bounds check */
                    qqbar_sqrt(xs + i, xs + i);

                    /* todo: avoid copy here... */
                    CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)->data.func_data.qqbar = flint_malloc(sizeof(qqbar_struct));
                    qqbar_init(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)->data.func_data.qqbar);
                    qqbar_set(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)->data.func_data.qqbar, xs + i);

                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Abs)
                {
                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    if (!ca_get_qqbar(xs + i, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)), ctx))
                        goto cleanup;

                    /* todo: bounds check */
                    qqbar_abs(xs + i, xs + i);
                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Sign)
                {
                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    if (!ca_get_qqbar(xs + i, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)), ctx))
                        goto cleanup;

                    /* todo: bounds check */
                    qqbar_sgn(xs + i, xs + i);
                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Re)
                {
                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    if (!ca_get_qqbar(xs + i, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)), ctx))
                        goto cleanup;

                    /* todo: bounds check */
                    qqbar_re(xs + i, xs + i);
                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Im)
                {
                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    if (!ca_get_qqbar(xs + i, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)), ctx))
                        goto cleanup;

                    /* todo: bounds check */
                    qqbar_im(xs + i, xs + i);
                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Conjugate)
                {
                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    if (!ca_get_qqbar(xs + i, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)), ctx))
                        goto cleanup;

                    qqbar_conj(xs + i, xs + i);
                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Floor)
                {
                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    if (!ca_get_qqbar(xs + i, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)), ctx))
                        goto cleanup;

                    {
                        fmpz_t t;
                        fmpz_init(t);
                        qqbar_floor(t, xs + i);
                        qqbar_set_fmpz(xs + i, t);
                        fmpz_clear(t);
                    }
                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Ceil)
                {
                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    if (!ca_get_qqbar(xs + i, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)), ctx))
                        goto cleanup;

                    {
                        fmpz_t t;
                        fmpz_init(t);
                        qqbar_ceil(t, xs + i);
                        qqbar_set_fmpz(xs + i, t);
                        fmpz_clear(t);
                    }
                }
                else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) == CA_Pow)
                {
                    ca_srcptr base, exp;

                    qqbar_init(xs + i);
                    init_mask[i] = 1;

                    base = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i));
                    exp = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i)) + 1;

                    if (!CA_IS_QQ(exp, ctx))
                        goto cleanup;

                    if (!ca_get_qqbar(xs + i, base, ctx))
                        goto cleanup;

                    /* todo: better condition on the exponent numerator */
                    if (fmpz_bits(fmpq_denref(CA_FMPQ(exp))) > 20 ||
                        fmpz_bits(fmpq_numref(CA_FMPQ(exp))) > 20 ||
                        *fmpq_denref(CA_FMPQ(exp)) * qqbar_degree(xs + i) > ctx->options[CA_OPT_QQBAR_DEG_LIMIT] ||
                        (qqbar_height_bits(xs + i) + 1) * fmpz_bits(fmpq_numref(CA_FMPQ(exp))) > ctx->options[CA_OPT_PREC_LIMIT])
                        goto cleanup;

                    qqbar_root_ui(xs + i, xs + i, *fmpq_denref(CA_FMPQ(exp)));
                    if (*fmpq_numref(CA_FMPQ(exp)) >= 0)
                    {
                        qqbar_pow_ui(xs + i, xs + i, *fmpq_numref(CA_FMPQ(exp)));
                    }
                    else
                    {
                        qqbar_inv(xs + i, xs + i);
                        qqbar_pow_ui(xs + i, xs + i, -*fmpq_numref(CA_FMPQ(exp)));
                    }
                }
            }

            if (qqbar_evaluate_fmpz_mpoly(y, fmpz_mpoly_q_numref(CA_MPOLY_Q(x)), xs, deg_limit, bits_limit, CA_FIELD_MCTX(CA_FIELD(x, ctx), ctx)))
            {
                if (qqbar_evaluate_fmpz_mpoly(res, fmpz_mpoly_q_denref(CA_MPOLY_Q(x)), xs, deg_limit, bits_limit, CA_FIELD_MCTX(CA_FIELD(x, ctx), ctx)))
                {
                    if (qqbar_binop_within_limits(y, res, deg_limit, bits_limit))
                    {
                        qqbar_div(res, y, res);
                        success = 1;
                    }
                }
            }

cleanup:
            for (i = 0; i < len; i++)
            {
                if (init_mask[i])
                    qqbar_clear(xs + i);
            }

            flint_free(init_mask);
            flint_free(used);
            flint_free(xs);
            qqbar_clear(y);
            qqbar_clear(zero);

            return success;
        }
    }
}

int
ca_get_fmpq(fmpq_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        return 0;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        fmpq_set(res, CA_FMPQ(x));
        return 1;
    }
    else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        if (nf_elem_is_rational(CA_NF_ELEM(x), CA_FIELD_NF(CA_FIELD(x, ctx))))
        {
            const fmpz * num;
            const fmpz * den;
            slong len;

            if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_LINEAR)
            {
                num = (fmpz *) LNF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = LNF_ELEM_DENREF(CA_NF_ELEM(x));
                fmpz_set(fmpq_numref(res), num);
                fmpz_set(fmpq_denref(res), den);
            }
            else if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_QUADRATIC)
            {
                num = (fmpz *) QNF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = QNF_ELEM_DENREF(CA_NF_ELEM(x));
                fmpz_set(fmpq_numref(res), num);
                fmpz_set(fmpq_denref(res), den);
            }
            else
            {
                num = (fmpz *) NF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = NF_ELEM_DENREF(CA_NF_ELEM(x));
                len = NF_ELEM(CA_NF_ELEM(x))->length;
                if (len == 0)
                {
                    fmpq_zero(res);
                }
                else
                {
                    fmpz_set(fmpq_numref(res), num);
                    fmpz_set(fmpq_denref(res), den);
                }
            }

            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        /* todo: exclude complex numbers, obviously irrational numbers
           before evaluating */
        int success;
        qqbar_t t;
        qqbar_init(t);
        success = ca_get_qqbar(t, x, ctx);
        if (success && qqbar_is_rational(t))
        {
            qqbar_get_fmpq(res, t);
            success = 1;
        }
        else
        {
            success = 0;
        }
        qqbar_clear(t);
        return success;
    }
}

int
ca_get_fmpz(fmpz_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        return 0;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        if (fmpz_is_one(fmpq_denref(CA_FMPQ(x))))
        {
            fmpz_set(res, fmpq_numref(CA_FMPQ(x)));
            return 1;
        }

        return 0;
    }
    else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        if (nf_elem_is_integer(CA_NF_ELEM(x), CA_FIELD_NF(CA_FIELD(x, ctx))))
        {
            const fmpz * num;
            slong len;

            if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_LINEAR)
            {
                num = (fmpz *) LNF_ELEM_NUMREF(CA_NF_ELEM(x));
                fmpz_set(res, num);
            }
            else if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_QUADRATIC)
            {
                num = (fmpz *) QNF_ELEM_NUMREF(CA_NF_ELEM(x));
                fmpz_set(res, num);
            }
            else
            {
                num = (fmpz *) NF_ELEM_NUMREF(CA_NF_ELEM(x));
                len = NF_ELEM(CA_NF_ELEM(x))->length;
                if (len == 0)
                    fmpz_zero(res);
                else
                    fmpz_set(res, num);
            }

            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        /* todo: exclude (numerically) obvious non-integers before evaluating */
        int success;
        qqbar_t t;
        qqbar_init(t);
        success = ca_get_qqbar(t, x, ctx);
        if (success && qqbar_is_integer(t))
        {
            qqbar_get_fmpz(res, t);
            success = 1;
        }
        else
        {
            success = 0;
        }
        qqbar_clear(t);
        return success;
    }
}

