/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

slong
nf_elem_bits(const nf_elem_t x, nf_t nf)
{
    slong b, c;

    if (nf->flag & NF_LINEAR)
    {
        flint_abort();
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        b = fmpz_bits(QNF_ELEM_NUMREF(x));
        c = fmpz_bits(QNF_ELEM_NUMREF(x) + 1);
        b = FLINT_MAX(b, c);
        c = fmpz_bits(QNF_ELEM_DENREF(x));
        b = FLINT_MAX(b, c);
        return b;
    }
    else
    {
        b = _fmpz_vec_max_bits(NF_ELEM_NUMREF(x), NF_ELEM(x)->length);
        b = FLINT_ABS(b);
        c = fmpz_bits(NF_ELEM_DENREF(x));
        b = FLINT_MAX(b, c);
        return b;
    }
}

void
_ca_pow_binexp(ca_t res, const ca_t x, slong n, ca_ctx_t ctx)
{
    if (n == 0)
    {
        ca_one(res, ctx);
    }
    else if (n == 1)
    {
        ca_set(res, x, ctx);
    }
    else if (n == 2)
    {
        ca_mul(res, x, x, ctx);
    }
    else if (n < 0)
    {
        ca_inv(res, x, ctx);
        _ca_pow_binexp(res, res, -n, ctx);
    }
    else
    {
        /* todo: is left-to-right or right-to-left better? */
        if (n % 2 == 0)
        {
#if 0
            _ca_pow_binexp(res, x, n / 2, ctx);
            _ca_pow_binexp(res, res, 2, ctx);
#else
            _ca_pow_binexp(res, x, 2, ctx);
            _ca_pow_binexp(res, res, n / 2, ctx);
#endif
        }
        else
        {
            if (res == x)
            {
                ca_t t;
                ca_init(t, ctx);
                _ca_pow_binexp(t, x, n - 1, ctx);
                ca_mul(res, t, x, ctx);
                ca_clear(t, ctx);
            }
            else
            {
                _ca_pow_binexp(res, x, n - 1, ctx);
                ca_mul(res, res, x, ctx);
            }
        }
    }
}

/* (z^a)^b, assuming z != 0 */
void
ca_pow_pow(ca_t res, const ca_t z, const ca_t a, const ca_t b, ca_ctx_t ctx)
{
    ca_t t, u, pi;

    if (CA_IS_SPECIAL(z) || CA_IS_SPECIAL(a) || CA_IS_SPECIAL(b))
        flint_abort();

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
        ca_mul(u, a, b, ctx);
        ca_pow(res, z, u, ctx);
    }
    else
    {
        ca_t pi_i;

        ca_init(pi_i, ctx);
        ca_pi_i(pi_i, ctx);

        ca_mul(t, t, pi_i, ctx);
        ca_mul_ui(t, t, 2, ctx);
        ca_mul(t, t, b, ctx);
        ca_neg(t, t, ctx);
        ca_exp(t, t, ctx);

        ca_mul(u, a, b, ctx);
        ca_pow(res, z, u, ctx);
        ca_mul(res, res, t, ctx);

        ca_clear(pi_i, ctx);
    }

    ca_clear(t, ctx);
    ca_clear(u, ctx);
    ca_clear(pi, ctx);
}

void
ca_pow(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    truth_t xzero;

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        ca_unknown(res, ctx);
        return;
    }

    xzero = ca_check_is_zero(x, ctx);

    if (xzero != T_FALSE)
    {
        if (xzero == T_TRUE)
        {
            if (CA_IS_QQ(y, ctx))
            {
                if (fmpq_is_zero(CA_FMPQ(y)))
                {
                    ca_one(res, ctx);
                    return;
                }
                else if (fmpq_sgn(CA_FMPQ(y)) > 0)
                {
                    ca_zero(res, ctx);
                    return;
                }
                else
                {
                    ca_uinf(res, ctx);
                    return;
                }
            }
        }

        ca_unknown(res, ctx);
    }
    else
    {
        if (CA_IS_QQ(y, ctx))
        {
            if (fmpz_is_one(CA_FMPQ_DENREF(y)))
            {
                if (fmpz_is_zero(CA_FMPQ_NUMREF(y)))
                {
                    ca_one(res, ctx);
                    return;
                }

                if (fmpz_is_one(CA_FMPQ_NUMREF(y)))
                {
                    ca_set(res, x, ctx);
                    return;
                }

                if (fmpz_equal_si(CA_FMPQ_NUMREF(y), -1))
                {
                    ca_inv(res, x, ctx);
                    return;
                }

                if (fmpz_equal_si(CA_FMPQ_NUMREF(y), 2))
                {
                    ca_mul(res, x, x, ctx);
                    return;
                }

                if (fmpz_equal_si(CA_FMPQ_NUMREF(y), -2))
                {
                    ca_inv(res, x, ctx);
                    ca_mul(res, res, res, ctx);
                    return;
                }

                if (CA_IS_QQ(x, ctx) && fmpz_bits(CA_FMPQ_NUMREF(y)) <= FLINT_BITS - 2)
                {
                    slong xbits1, xbits2;

                    xbits1 = fmpz_bits(CA_FMPQ_NUMREF(x));
                    xbits2 = fmpz_bits(CA_FMPQ_DENREF(x));
                    xbits1 = FLINT_MAX(xbits1, xbits2);

                    if (xbits1 * (double) FLINT_ABS(*CA_FMPQ_NUMREF(y)) < ctx->options[CA_OPT_PREC_LIMIT])
                    {
                        fmpq_t t;
                        fmpq_init(t);
                        fmpq_pow_si(t, CA_FMPQ(x), *CA_FMPQ_NUMREF(y));
                        ca_set_fmpq(res, t, ctx);
                        fmpq_clear(t);
                        return;
                    }
                }

                if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)) && fmpz_bits(CA_FMPQ_NUMREF(y)) <= FLINT_BITS - 2)
                {
                    slong xbits1;

                    xbits1 = nf_elem_bits(CA_NF_ELEM(x), CA_FIELD_NF(CA_FIELD(x, ctx)));

                    if (xbits1 * (double) FLINT_ABS(*CA_FMPQ_NUMREF(y)) < ctx->options[CA_OPT_PREC_LIMIT])
                    {
                        ca_t t;
                        ca_init(t, ctx);

                        if (fmpq_sgn(CA_FMPQ(y)) > 0)
                            ca_set(t, x, ctx);
                        else
                            ca_inv(t, x, ctx);

                        nf_elem_pow(CA_NF_ELEM(t), CA_NF_ELEM(t), FLINT_ABS(*CA_FMPQ_NUMREF(y)), CA_FIELD_NF(CA_FIELD(t, ctx)));
                        ca_condense_field(t, ctx);
                        ca_swap(res, t, ctx);

                        ca_clear(t, ctx);
                        return;
                    }
                }

                /* todo: evaluation limits */
                if (fmpz_cmp_si(CA_FMPQ_NUMREF(y), -ctx->options[CA_OPT_POW_LIMIT]) >= 0 && fmpz_cmp_si(CA_FMPQ_NUMREF(y), ctx->options[CA_OPT_POW_LIMIT]) <= 0)
                {
                    _ca_pow_binexp(res, x, *CA_FMPQ_NUMREF(y), ctx);
                    return;
                }
            }
            else if (fmpz_equal_ui(CA_FMPQ_DENREF(y), 2))
            {
                if (fmpz_equal_si(CA_FMPQ_NUMREF(y), 1))
                {
                    ca_sqrt(res, x, ctx);
                    return;
                }
                else if (fmpz_equal_si(CA_FMPQ_NUMREF(y), 3))
                {
                    ca_t t;
                    ca_init(t, ctx);
                    ca_sqrt(t, x, ctx);
                    ca_mul(res, t, x, ctx);
                    ca_clear(t, ctx);
                    return;
                }  /* todo: evaluation limits */
                else if (fmpz_cmp_si(CA_FMPQ_NUMREF(y), -ctx->options[CA_OPT_POW_LIMIT] / 2) >= 0 && fmpz_cmp_si(CA_FMPQ_NUMREF(y), ctx->options[CA_OPT_POW_LIMIT] / 2) <= 0)
                {
                    ca_t t;
                    ca_init(t, ctx);
                    ca_sqrt(t, x, ctx);
                    _ca_pow_binexp(res, t, *CA_FMPQ_NUMREF(y), ctx);
                    ca_clear(t, ctx);
                    return;
                }
            }
        }

        {
            ca_ext_ptr ext = ca_is_gen_as_ext(x, ctx);

            /* (z^a)^b */
            if (ext != NULL && CA_EXT_HEAD(ext) == CA_Pow)
            {
                if (ca_check_is_zero(CA_EXT_FUNC_ARGS(ext), ctx) == T_FALSE)
                {
                    ca_pow_pow(res, CA_EXT_FUNC_ARGS(ext), CA_EXT_FUNC_ARGS(ext) + 1, y, ctx);
                    return;
                }
            }
        }

        _ca_function_fxy(res, CA_Pow, x, y, ctx);
        /* todo: detect simple values instead of creating extension element in the first place */
        _ca_mpoly_q_reduce_ideal(CA_MPOLY_Q(res), CA_FIELD(res, ctx), ctx);
        ca_condense_field(res, ctx);
    }
}

void
ca_pow_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_fmpz(t, y, ctx);
    ca_pow(res, x, t, ctx);
    ca_clear(t, ctx);
}

void
ca_pow_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_fmpq(t, y, ctx);
    ca_pow(res, x, t, ctx);
    ca_clear(t, ctx);
}

void
ca_pow_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_si(t, y, ctx);
    ca_pow(res, x, t, ctx);
    ca_clear(t, ctx);
}

void
ca_pow_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_ui(t, y, ctx);
    ca_pow(res, x, t, ctx);
    ca_clear(t, ctx);
}
