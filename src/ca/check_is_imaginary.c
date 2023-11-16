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
ca_check_is_imaginary(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        if (fmpq_is_zero(CA_FMPQ(x)))
            return T_TRUE;
        else
            return T_FALSE;
    }
    else if (CA_IS_QQ_I(x, ctx))
    {
        const fmpz *n;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));

        if (fmpz_is_zero(n))
            return T_TRUE;

        return T_FALSE;
    }
    else
    {
        acb_t t;
        truth_t res;
        slong prec, prec_limit;

        res = T_UNKNOWN;

        acb_init(t);

        prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
        prec_limit = FLINT_MAX(prec_limit, 64);

        for (prec = 64; (prec <= prec_limit) && (res == T_UNKNOWN); prec *= 2)
        {
            ca_get_acb_raw(t, x, prec, ctx);

            if (arb_is_zero(acb_realref(t)))
            {
                res = T_TRUE;
                break;
            }

            if (!arb_contains_zero(acb_realref(t)))
            {
                res = T_FALSE;
                break;
            }

            /* try conjugation */
            /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
            if (prec == 64)
            {
                ca_t t;
                ca_init(t, ctx);
                ca_conj_deep(t, x, ctx);
                ca_neg(t, t, ctx);
                res = ca_check_equal(t, x, ctx);
                ca_clear(t, ctx);
                if (res != T_UNKNOWN)
                    break;
            }

            /* try qqbar computation */
            /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
            if (prec == 64)
            {
                qqbar_t a;
                qqbar_init(a);

                if (ca_get_qqbar(a, x, ctx))
                    res = (qqbar_sgn_re(a) == 0) ? T_TRUE : T_FALSE;

                qqbar_clear(a);
            }
        }

        acb_clear(t);

        return res;
    }
}
