/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

ca_field_ptr
_ca_vec_same_field2(ca_srcptr A, slong Alen, ca_srcptr B, slong Blen, ca_ctx_t ctx);

void
_ca_poly_div_series(ca_ptr Q,
    ca_srcptr A, slong Alen,
    ca_srcptr B, slong Blen, slong len, ca_ctx_t ctx)
{
    Alen = FLINT_MIN(Alen, len);
    Blen = FLINT_MIN(Blen, len);

    if (CA_IS_SPECIAL(A) || CA_IS_SPECIAL(B))
    {
        if (ca_is_unknown(A, ctx) || ca_is_unknown(B, ctx))
            _ca_vec_unknown(Q, len, ctx);
        else
            _ca_vec_undefined(Q, len, ctx);
        return;
    }

    /* todo: tuning */
    /* todo: shallow copy with integers */
    if (Blen >= 4 && _ca_vec_is_fmpq_vec(A, Alen, ctx)
        && _ca_vec_is_fmpq_vec(B, Blen, ctx) && !fmpq_is_zero(CA_FMPQ(B)))
    {
        fmpz *p, *r, *s;
        fmpz_t pden, rden, sden;

        p = _fmpz_vec_init(Alen);
        r = _fmpz_vec_init(Blen);
        s = _fmpz_vec_init(len);
        fmpz_init(pden);
        fmpz_init(rden);
        fmpz_init(sden);

        _ca_vec_fmpq_vec_get_fmpz_vec_den(p, pden, A, Alen, ctx);
        _ca_vec_fmpq_vec_get_fmpz_vec_den(r, rden, B, Blen, ctx);
        _fmpq_poly_div_series(s, sden, p, pden, Alen, r, rden, Blen, len);
        _ca_vec_set_fmpz_vec_div_fmpz(Q, s, sden, len, ctx);

        fmpz_clear(pden);
        fmpz_clear(rden);
        fmpz_clear(sden);
        _fmpz_vec_clear(p, Alen);
        _fmpz_vec_clear(r, Blen);
        _fmpz_vec_clear(s, len);
        return;
    }

    if (Blen == 1)
    {
        _ca_vec_scalar_div_ca(Q, A, Alen, B, ctx);
        _ca_vec_zero(Q + Alen, len - Alen, ctx);
    }
    else
    {
        int is_one;
        ca_field_ptr K;
        ca_t q;
        slong i;

        if (Blen > 8)
        {
            K = _ca_vec_same_field2(A, Alen, B, Blen, ctx);

            /* If we have fast multiplication and inversion */
            /* Todo: also want this if inversion alone is fast and Alen is small */
            if (K != NULL && CA_FIELD_IS_NF(K))
            {
                if (len > 2 * qqbar_degree(CA_FIELD_NF_QQBAR(K)))
                {
                    ca_ptr Binv;
                    Binv = _ca_vec_init(len, ctx);
                    _ca_poly_inv_series(Binv, B, Blen, len, ctx);
                    _ca_poly_mullow(Q, Binv, len, A, Alen, len, ctx);
                    _ca_vec_clear(Binv, len, ctx);
                    return;
                }
            }
        }

        ca_init(q, ctx);

        ca_inv(q, B + 0, ctx);
        ca_mul(Q, A + 0, q, ctx);

        is_one = (ca_check_is_one(q, ctx) == T_TRUE);

        for (i = 1; i < len; i++)
        {
            ca_dot(Q + i, (i < Alen) ? A + i : NULL, 1,
                B + 1, 1, Q + i - 1, -1, FLINT_MIN(i, Blen - 1), ctx);
            if (!is_one)
                ca_mul(Q + i, Q + i, q, ctx);
        }

        ca_clear(q, ctx);
    }
}

void
ca_poly_div_series(ca_poly_t Q, const ca_poly_t A, const ca_poly_t B, slong len, ca_ctx_t ctx)
{
    if (len == 0)
    {
        ca_poly_zero(Q, ctx);
        return;
    }

    if (B->length == 0)
    {
        ca_poly_fit_length(Q, len, ctx);
        ca_unknown(Q->coeffs, ctx);  /* todo: uinf when A != 0? */
        _ca_vec_undefined(Q->coeffs + 1, len - 1, ctx);
        _ca_poly_set_length(Q, len, ctx);
        return;
    }

    if (A->length == 0)
    {
        if (ca_check_is_zero(B->coeffs, ctx) == T_FALSE)
        {
            ca_poly_zero(Q, ctx);
        }
        else
        {
            ca_poly_fit_length(Q, len, ctx);
            _ca_vec_unknown(Q->coeffs, len, ctx);
            _ca_poly_set_length(Q, len, ctx);
        }
        return;
    }

    if (Q == A || Q == B)
    {
        ca_poly_t t;
        ca_poly_init(t, ctx);
        ca_poly_div_series(t, A, B, len, ctx);
        ca_poly_swap(Q, t, ctx);
        ca_poly_clear(t, ctx);
        return;
    }

    ca_poly_fit_length(Q, len, ctx);
    _ca_poly_div_series(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, len, ctx);
    _ca_poly_set_length(Q, len, ctx);
    _ca_poly_normalise(Q, ctx);
}
