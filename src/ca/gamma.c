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
ca_gamma_inert(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }
    else
    {
        _ca_function_fx(res, CA_Gamma, x, ctx);
    }
}

void
_ca_gamma_verbatim(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t y, tmp;

    ca_init(y, ctx);
    ca_init(tmp, ctx);

    _ca_function_fx(y, CA_Gamma, x, ctx);
    ca_merge_fields(tmp, res, x, y, ctx);

    ca_clear(y, ctx);
    ca_clear(tmp, ctx);
}

truth_t
ca_re_is_positive(const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u;
    truth_t res;

    ca_init(t, ctx);
    ca_init(u, ctx);
    ca_re(t, x, ctx);
    res = ca_check_gt(t, u, ctx);
    ca_clear(t, ctx);
    ca_clear(u, ctx);
    return res;
}

void
ca_gamma(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        /* todo: other signs */
        if (ca_check_is_pos_inf(x, ctx))
        {
            ca_pos_inf(res, ctx);
        }
        else if (ca_check_is_neg_inf(x, ctx) || ca_check_is_undefined(x, ctx) || ca_check_is_uinf(x, ctx))
        {
            ca_undefined(res, ctx);
        }
        else if (ca_check_is_pos_i_inf(x, ctx) || ca_check_is_neg_i_inf(x, ctx))
        {
            ca_zero(res, ctx);
        }
        else
        {
            ca_unknown(res, ctx);
        }
    }
    else
    {
        truth_t is_integer, is_positive;
        slong expand_limit;

        expand_limit = ctx->options[CA_OPT_PREC_LIMIT];

        is_integer = ca_check_is_integer(x, ctx);

        if (is_integer == T_TRUE)
        {
            is_positive = ca_re_is_positive(x, ctx);

            if (is_positive == T_TRUE)
            {
                fmpz_t t;
                fmpz_init(t);

                if (ca_get_fmpz(t, x, ctx) && fmpz_cmp_ui(t, expand_limit) < 0)
                {
                    fmpz_fac_ui(t, fmpz_get_ui(t) - 1);
                    ca_set_fmpz(res, t, ctx);
                }
                else
                {
                    /* todo: may want to use the potentially simplified x (to an integer) here... */
                    _ca_gamma_verbatim(res, x, ctx);
                }

                fmpz_clear(t);
            }
            else if (is_positive == T_FALSE)
            {
                ca_uinf(res, ctx);
            }
            else
            {
                ca_unknown(res, ctx);
            }
        }
        else if (is_integer == T_FALSE)
        {
            ca_t y;
            fmpz_t t, u;
            slong n;

            ca_init(y, ctx);
            fmpz_init(t);
            fmpz_init(u);

            ca_set_d(y, 0.5, ctx);
            ca_sub(y, x, y, ctx);

            if (ca_get_fmpz(t, y, ctx) && fmpz_cmp_si(t, expand_limit) < 0 &&
                                          fmpz_cmp_si(t, -expand_limit) > 0)
            {
                n = fmpz_get_si(t);

                ca_pi(res, ctx);
                ca_sqrt(res, res, ctx);

                if (n != 0)
                {
                    if (n >= 0)
                    {
                        fmpz_fac_ui(t, 2 * n);
                        fmpz_fac_ui(u, n);
                        fmpz_mul_2exp(u, u, 2 * n);
                    }
                    else
                    {
                        fmpz_fac_ui(t, -n);
                        fmpz_mul_2exp(t, t, -2 * n);
                        fmpz_fac_ui(u, -2 * n);
                        if (n & 1)
                            fmpz_neg(t, t);
                    }

                    ca_mul_fmpz(res, res, t, ctx);
                    ca_div_fmpz(res, res, u, ctx);
                }
            }
            else
            {
                _ca_gamma_verbatim(res, x, ctx);
            }

            ca_clear(y, ctx);
            fmpz_clear(t);
            fmpz_clear(u);
        }
        else
        {
            is_positive = ca_re_is_positive(x, ctx);

            if (is_positive == T_TRUE)
            {
                _ca_gamma_verbatim(res, x, ctx);
            }
            else
            {
                ca_unknown(res, ctx);
            }
        }
    }
}
