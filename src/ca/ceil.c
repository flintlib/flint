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
ca_ceil(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_UNKNOWN(x))
            ca_unknown(res, ctx);
        else
            ca_undefined(res, ctx);
        return;
    }

    if (CA_IS_QQ(x, ctx))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_cdiv_q(t, CA_FMPQ_NUMREF(x), CA_FMPQ_DENREF(x));
        ca_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return;
    }

    {
        slong prec, prec_limit;
        acb_t v;
        mag_t m;
        fmpz_t n;
        int success = 0;

        acb_init(v);
        mag_init(m);
        fmpz_init(n);

        prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
        prec_limit = FLINT_MAX(prec_limit, 64);

        for (prec = 64; (prec <= prec_limit) && !success; prec *= 2)
        {
            ca_get_acb_raw(v, x, prec, ctx);
            arb_get_mag(m, acb_realref(v));

            if (arb_is_finite(acb_imagref(v)) && mag_cmp_2exp_si(m, prec_limit) <= 0)
            {
                arb_ceil(acb_realref(v), acb_realref(v), prec);

                if (arb_get_unique_fmpz(n, acb_realref(v)))
                {
                    ca_set_fmpz(res, n, ctx);
                    success = 1;
                    break;
                }
            }

            arb_get_mag_lower(m, acb_realref(v));

            if (mag_cmp_2exp_si(m, prec_limit) > 0)
                break;
        }

        acb_clear(v);
        mag_clear(m);
        fmpz_clear(n);

        if (success)
            return;
    }

    _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_Ceil, x), ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
}
