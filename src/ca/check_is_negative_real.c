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
ca_check_is_negative_real(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        return (fmpq_sgn(CA_FMPQ(x)) < 0) ? T_TRUE : T_FALSE;
    }
    else if (CA_IS_QQ_I(x, ctx))
    {
        const fmpz *n;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));

        if (fmpz_is_zero(n + 1))
            return (fmpz_sgn(n) < 0) ? T_TRUE : T_FALSE;

        return T_FALSE;
    }
    else
    {
        acb_t t;
        truth_t res, is_real;
        slong prec, prec_limit;

        res = T_UNKNOWN;

        acb_init(t);

        prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
        prec_limit = FLINT_MAX(prec_limit, 64);

        is_real = T_UNKNOWN;

        for (prec = 64; (prec <= prec_limit) && (res == T_UNKNOWN); prec *= 2)
        {
            ca_get_acb_raw(t, x, prec, ctx);

            if (is_real == T_UNKNOWN)
            {
                if (arb_is_zero(acb_imagref(t)))
                    is_real = T_TRUE;
                else if (!arb_contains_zero(acb_imagref(t)))
                    is_real = T_FALSE;
            }

            if ((is_real == T_TRUE) && arb_is_negative(acb_realref(t)))
            {
                res = T_TRUE;
                break;
            }

            if ((is_real == T_FALSE) || arb_is_nonnegative(acb_realref(t)))
            {
                res = T_FALSE;
                break;
            }

            if (prec == 64 && is_real == T_UNKNOWN)
            {
                ca_t t;
                ca_init(t, ctx);
                ca_conj_deep(t, x, ctx);
                is_real = ca_check_equal(t, x, ctx);
                ca_clear(t, ctx);
                if (is_real == T_FALSE)
                {
                    res = T_FALSE;
                    break;
                }
            }

            /* try qqbar computation */
            /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
            if (prec == 64)
            {
                qqbar_t a;
                qqbar_init(a);

                if (ca_get_qqbar(a, x, ctx))
                    res = (qqbar_sgn_im(a) == 0 && qqbar_sgn_re(a) < 0) ? T_TRUE : T_FALSE;

                qqbar_clear(a);
            }
        }

        acb_clear(t);

        return res;
    }
}

