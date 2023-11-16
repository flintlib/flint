/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "ca_poly.h"

ca_field_ptr
_ca_vec_same_field2(ca_srcptr A, slong Alen, ca_srcptr B, slong Blen, ca_ctx_t ctx);

void
_ca_poly_inv_series(ca_ptr Qinv,
    ca_srcptr Q, slong Qlen, slong len, ca_ctx_t ctx)
{
    Qlen = FLINT_MIN(Qlen, len);

    if (CA_IS_SPECIAL(Q))
    {
        if (ca_is_unknown(Q, ctx))
            _ca_vec_unknown(Qinv, len, ctx);
        else
            _ca_vec_undefined(Qinv, len, ctx);
        return;
    }

    /* todo: tuning */
    /* todo: shallow copy with integers */
    if (Qlen >= 4 && _ca_vec_is_fmpq_vec(Q, Qlen, ctx) && !fmpq_is_zero(CA_FMPQ(Q)))
    {
        fmpz *p, *r;
        fmpz_t pden, rden;

        p = _fmpz_vec_init(Qlen);
        r = _fmpz_vec_init(len);
        fmpz_init(pden);
        fmpz_init(rden);

        _ca_vec_fmpq_vec_get_fmpz_vec_den(p, pden, Q, Qlen, ctx);
        _fmpq_poly_inv_series(r, rden, p, pden, Qlen, len);
        _ca_vec_set_fmpz_vec_div_fmpz(Qinv, r, rden, len, ctx);

        fmpz_clear(pden);
        fmpz_clear(rden);
        _fmpz_vec_clear(p, Qlen);
        _fmpz_vec_clear(r, len);
        return;
    }

    ca_inv(Qinv, Q, ctx);

    if (CA_IS_SPECIAL(Qinv))
    {
        if (ca_is_unknown(Qinv, ctx))
            _ca_vec_unknown(Qinv + 1, len - 1, ctx);
        else
            _ca_vec_undefined(Qinv + 1, len - 1, ctx);
        return;
    }

    if (Qlen == 1)
    {
        _ca_vec_zero(Qinv + 1, len - 1, ctx);
    }
    else if (len == 2)
    {
        ca_mul(Qinv + 1, Qinv, Qinv, ctx);
        ca_mul(Qinv + 1, Qinv + 1, Q + 1, ctx);
        ca_neg(Qinv + 1, Qinv + 1, ctx);
    }
    else
    {
        int is_one;
        ca_field_ptr K;
        slong i, blen;

        if (Qlen <= 8)
        {
            blen = len;
        }
        else
        {
            K = _ca_vec_same_field2(Q, Qlen, NULL, 0, ctx);

            /* Newton iteration where we have fast multiplication */
            if (K != NULL && CA_FIELD_IS_NF(K))
            {
                blen = 2 * qqbar_degree(CA_FIELD_NF_QQBAR(K));
                blen = FLINT_MIN(blen, len);
            }
            else
            {
                blen = len;
            }
        }

        is_one = (ca_check_is_one(Qinv, ctx) == T_TRUE);

        for (i = 1; i < blen; i++)
        {
            ca_dot(Qinv + i, NULL, 1,
                Q + 1, 1, Qinv + i - 1, -1, FLINT_MIN(i, Qlen - 1), ctx);
            if (!is_one)
                ca_mul(Qinv + i, Qinv + i, Qinv, ctx);
        }

        if (len > blen)
        {
            slong Qnlen, Wlen, W2len;
            ca_ptr W;

            W = _ca_vec_init(len, ctx);

            NEWTON_INIT(blen, len)
            NEWTON_LOOP(m, n)

            Qnlen = FLINT_MIN(Qlen, n);
            Wlen = FLINT_MIN(Qnlen + m - 1, n);
            W2len = Wlen - m;
            _ca_poly_mullow(W, Q, Qnlen, Qinv, m, Wlen, ctx);
            _ca_poly_mullow(Qinv + m, Qinv, m, W + m, W2len, n - m, ctx);
            _ca_vec_neg(Qinv + m, Qinv + m, n - m, ctx);

            NEWTON_END_LOOP
            NEWTON_END

            _ca_vec_clear(W, len, ctx);
        }
    }
}

void
ca_poly_inv_series(ca_poly_t Qinv, const ca_poly_t Q, slong len, ca_ctx_t ctx)
{
    if (len == 0)
    {
        ca_poly_zero(Qinv, ctx);
        return;
    }

    if (Q->length == 0)
    {
        ca_poly_fit_length(Qinv, len, ctx);
        ca_uinf(Qinv->coeffs, ctx);
        _ca_vec_undefined(Qinv->coeffs + 1, len - 1, ctx);
        _ca_poly_set_length(Qinv, len, ctx);
        return;
    }

    if (Qinv == Q)
    {
        ca_poly_t t;
        ca_poly_init(t, ctx);
        ca_poly_inv_series(t, Q, len, ctx);
        ca_poly_swap(Qinv, t, ctx);
        ca_poly_clear(t, ctx);
        return;
    }

    ca_poly_fit_length(Qinv, len, ctx);
    _ca_poly_inv_series(Qinv->coeffs, Q->coeffs, Q->length, len, ctx);
    _ca_poly_set_length(Qinv, len, ctx);
    _ca_poly_normalise(Qinv, ctx);
}
