/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

ca_ext_ptr
ca_is_gen_pow_fmpz_as_ext(fmpz_t exp, const ca_t x, ca_ctx_t ctx)
{
    ca_field_ptr K;

    if (CA_IS_SPECIAL(x))
        return NULL;

    K = CA_FIELD(x, ctx);

    if (CA_FIELD_IS_QQ(K))
        return NULL;

    /* todo: detect powers of the generator in number fields */
    if (CA_FIELD_IS_NF(K))
    {
        if (!nf_elem_is_gen(CA_NF_ELEM(x), CA_FIELD_NF(K)))
            return NULL;

        fmpz_one(exp);
        return CA_FIELD_EXT_ELEM(K, 0);
    }

    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(CA_MPOLY_Q(x)), CA_FIELD_MCTX(K, ctx)))
    {
        if (fmpz_mpoly_length(fmpz_mpoly_q_numref(CA_MPOLY_Q(x)), CA_FIELD_MCTX(K, ctx)) == 1 &&
            fmpz_is_one(fmpz_mpoly_q_numref(CA_MPOLY_Q(x))->coeffs))
        {
            int * used;
            slong i, which, num_used;

            used = flint_calloc(CA_FIELD_LENGTH(K), sizeof(int));
            fmpz_mpoly_q_used_vars_num(used, CA_MPOLY_Q(x), CA_FIELD_MCTX(K, ctx));

            which = num_used = 0;
            for (i = 0; i < CA_FIELD_LENGTH(K); i++)
            {
                if (used[i])
                {
                    which = i;
                    num_used++;
                }
            }

            flint_free(used);

            if (num_used == 1)
            {
                fmpz_mpoly_total_degree_fmpz(exp, fmpz_mpoly_q_numref(CA_MPOLY_Q(x)),  CA_FIELD_MCTX(K, ctx));
                return CA_FIELD_EXT_ELEM(K, which);
            }
        }
    }

    if (fmpz_mpoly_is_one(fmpz_mpoly_q_numref(CA_MPOLY_Q(x)), CA_FIELD_MCTX(K, ctx)))
    {
        if (fmpz_mpoly_length(fmpz_mpoly_q_denref(CA_MPOLY_Q(x)), CA_FIELD_MCTX(K, ctx)) == 1 &&
            fmpz_is_one(fmpz_mpoly_q_denref(CA_MPOLY_Q(x))->coeffs))
        {
            int * used;
            slong i, which, num_used;

            used = flint_calloc(CA_FIELD_LENGTH(K), sizeof(int));
            fmpz_mpoly_q_used_vars_den(used, CA_MPOLY_Q(x), CA_FIELD_MCTX(K, ctx));

            which = num_used = 0;
            for (i = 0; i < CA_FIELD_LENGTH(K); i++)
            {
                if (used[i])
                {
                    which = i;
                    num_used++;
                }
            }

            flint_free(used);

            if (num_used == 1)
            {
                fmpz_mpoly_total_degree_fmpz(exp, fmpz_mpoly_q_denref(CA_MPOLY_Q(x)),  CA_FIELD_MCTX(K, ctx));
                fmpz_neg(exp, exp);
                return CA_FIELD_EXT_ELEM(K, which);
            }
        }
    }

    return NULL;
}


/* log(exp(z)) -- https://fungrim.org/entry/a3a253/ */
void
ca_log_exp(ca_t res, const ca_t z, ca_ctx_t ctx)
{
    ca_t t, pi;

    if (CA_IS_SPECIAL(z))
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    ca_init(t, ctx);
    ca_init(pi, ctx);

    ca_pi(pi, ctx);

    ca_im(t, z, ctx);
    ca_div(t, t, pi, ctx);
    ca_sub_ui(t, t, 1, ctx);
    ca_div_ui(t, t, 2, ctx);

    ca_ceil(t, t, ctx);

    if (ca_check_is_zero(t, ctx) == T_TRUE)
    {
        ca_set(res, z, ctx);
    }
    else
    {
        ca_t pi_i;

        ca_init(pi_i, ctx);
        ca_pi_i(pi_i, ctx);

        ca_mul(t, t, pi_i, ctx);
        ca_mul_ui(t, t, 2, ctx);

        ca_sub(res, z, t, ctx);
        ca_clear(pi_i, ctx);
    }

    ca_clear(t, ctx);
    ca_clear(pi, ctx);
}

/* log(z^a), assuming z != 0 */
void
ca_log_pow(ca_t res, const ca_t z, const ca_t a, ca_ctx_t ctx)
{
    ca_t t, u, pi;

    if (CA_IS_SPECIAL(z) || CA_IS_SPECIAL(a))
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    ca_init(t, ctx);
    ca_init(u, ctx);
    ca_init(pi, ctx);

    ca_log(u, z, ctx);
    ca_mul(u, u, a, ctx);

    ca_pi(pi, ctx);

    ca_im(t, u, ctx);
    ca_div(t, t, pi, ctx);
    ca_sub_ui(t, t, 1, ctx);
    ca_div_ui(t, t, 2, ctx);

    ca_ceil(t, t, ctx);

    if (ca_check_is_zero(t, ctx) == T_TRUE)
    {
        ca_set(res, u, ctx);
    }
    else
    {
        ca_t pi_i;

        ca_init(pi_i, ctx);
        ca_pi_i(pi_i, ctx);

        ca_mul(t, t, pi_i, ctx);
        ca_mul_ui(t, t, 2, ctx);

        ca_sub(res, u, t, ctx);
        ca_clear(pi_i, ctx);
    }

    ca_clear(t, ctx);
    ca_clear(u, ctx);
    ca_clear(pi, ctx);
}

void
ca_log(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    truth_t is_zero;
    ca_ext_ptr ext;

    if (CA_IS_SPECIAL(x))
    {
        if (ca_check_is_infinity(x, ctx) == T_TRUE)
            ca_pos_inf(res, ctx);
        else if (ca_check_is_undefined(x, ctx) == T_TRUE)
            ca_undefined(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }

    is_zero = ca_check_is_zero(x, ctx);

    if (is_zero == T_TRUE)
    {
        ca_neg_inf(res, ctx);
        return;
    }

    if (is_zero == T_UNKNOWN)
    {
        ca_unknown(res, ctx);
        return;
    }

    if (ca_check_is_one(x, ctx) == T_TRUE)
    {
        ca_zero(res, ctx);
        return;
    }

    /* log(i) and log(-i) */
    if (CA_IS_QQ_I(x, ctx))
    {
        if (ca_check_is_i(x, ctx) == T_TRUE)
        {
            ca_pi_i(res, ctx);
            ca_div_si(res, res, 2, ctx);
            return;
        }
        else if (ca_check_is_neg_i(x, ctx) == T_TRUE)
        {
            ca_pi_i(res, ctx);
            ca_div_si(res, res, -2, ctx);
            return;
        }
    }

    ext = ca_is_gen_as_ext(x, ctx);

    /* Fast detection of roots of unity. Todo: also detect roots
       of unity when in a number field, and in other situations. */
    if (ext != NULL && CA_EXT_HEAD(ext) == CA_QQBar)
    {
        slong p;
        ulong q;

        if (qqbar_log_pi_i(&p, &q, CA_EXT_QQBAR(ext)))
        {
            ca_pi_i(res, ctx);
            ca_mul_si(res, res, p, ctx);
            ca_div_ui(res, res, q, ctx);
            return;
        }
    }

    if (ext != NULL && CA_EXT_HEAD(ext) == CA_Exp)
    {
        /* log(exp(z)) */
        ca_log_exp(res, CA_EXT_FUNC_ARGS(ext), ctx);
        return;
    }

    if (ca_check_is_negative_real(x, ctx) == T_TRUE)
    {
        ca_t pi_i;
        ca_init(pi_i, ctx);
        ca_pi_i(pi_i, ctx);
        ca_neg(res, x, ctx);
        ca_log(res, res, ctx);
        ca_add(res, res, pi_i, ctx);
        ca_clear(pi_i, ctx);
        return;
    }

    if (ext != NULL && CA_EXT_HEAD(ext) == CA_Pow)
    {
        /* log(z^a) */
        if (ca_check_is_zero(CA_EXT_FUNC_ARGS(ext), ctx) == T_FALSE)
        {
            ca_log_pow(res, CA_EXT_FUNC_ARGS(ext), CA_EXT_FUNC_ARGS(ext) + 1, ctx);
            return;
        }
    }

    if (ext != NULL && CA_EXT_HEAD(ext) == CA_Sqrt)
    {
        /* log(sqrt(z)) */
        ca_t h;
        ca_init(h, ctx);
        ca_one(h, ctx);
        ca_div_ui(h, h, 2, ctx);
        ca_log_pow(res, CA_EXT_FUNC_ARGS(ext), h, ctx);
        ca_clear(h, ctx);
        return;
    }

    /* If x is not a generator, maybe it is a power of a generator */
    {
        fmpz_t n;
        ca_t t;
        int success = 0;

        fmpz_init(n);
        ca_init(t, ctx);

        ext = ca_is_gen_pow_fmpz_as_ext(n, x, ctx);

        if (ext != NULL && CA_EXT_HEAD(ext) == CA_Exp)
        {
            /* log(exp(z)^n) = log(exp(n*z)) */
            ca_mul_fmpz(t, CA_EXT_FUNC_ARGS(ext), n, ctx);
            ca_log_exp(res, t, ctx);
            success = 1;
        }

        if (ext != NULL && CA_EXT_HEAD(ext) == CA_Pow)
        {
            /* log((z^a)^n) = log(z^(a*n)) */
            if (ca_check_is_zero(CA_EXT_FUNC_ARGS(ext), ctx) == T_FALSE)
            {
                ca_mul_fmpz(t, CA_EXT_FUNC_ARGS(ext) + 1, n, ctx);
                ca_log_pow(res, CA_EXT_FUNC_ARGS(ext), t, ctx);
                success = 1;
            }
        }

        if (ext != NULL && CA_EXT_HEAD(ext) == CA_Sqrt)
        {
            /* log(sqrt(z)^n) = log(z^(n/2)) */
            ca_set_fmpz(t, n, ctx);
            ca_div_ui(t, t, 2, ctx);
            ca_log_pow(res, CA_EXT_FUNC_ARGS(ext), t, ctx);
            success = 1;
        }

        fmpz_clear(n);
        ca_clear(t, ctx);

        if (success)
            return;
    }

    _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_Log, x), ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
}
