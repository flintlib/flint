/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_pow(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        ca_unknown(res, ctx);
    }
    else if (ca_check_is_zero(x, ctx) != T_FALSE)
    {
        /* todo */
        ca_unknown(res, ctx);
    }
    else
    {
        if (CA_FIELD_IS_QQ(y, ctx))
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
                    ca_mul(res, x, x, ctx);
                    return;
                }

                if (CA_FIELD_IS_QQ(x, ctx) && fmpz_bits(CA_FMPQ_NUMREF(y)) <= FLINT_BITS - 2)
                {
                    slong xbits1, xbits2;

                    xbits1 = fmpz_bits(CA_FMPQ_NUMREF(x));
                    xbits2 = fmpz_bits(CA_FMPQ_DENREF(x));
                    xbits1 = FLINT_MAX(xbits1, xbits2);

                    if (xbits1 * (double) *CA_FMPQ_NUMREF(y) < ctx->options[CA_OPT_PREC_LIMIT])
                    {
                        fmpq_t t;
                        fmpq_init(t);
                        fmpq_pow_si(t, CA_FMPQ(x), *CA_FMPQ_NUMREF(y));
                        ca_set_fmpq(res, t, ctx);
                        fmpq_clear(t);
                        return;
                    }
                }
            }
        }

        _ca_function_fxy(res, CA_Pow, x, y, ctx);
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

