/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

int
_ca_vec_is_fmpq_vec(ca_srcptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!CA_IS_QQ(vec + i, ctx))
                return 0;
    return 1;
}

int
_ca_vec_fmpq_vec_is_fmpz_vec(ca_srcptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!fmpz_is_one(CA_FMPQ_DENREF(vec + i)))
                return 0;
    return 1;
}

int
_ca_poly_roots(ca_ptr roots, ca_srcptr poly, slong len, ulong flags, ca_ctx_t ctx)
{
    slong deg;
    truth_t leading_zero;

    if (len == 0)
        return 0;

    deg = len - 1;

    leading_zero = ca_check_is_zero(poly + deg, ctx);

    if (leading_zero != T_FALSE)
        return 0;

    if (deg == 0)
        return 1;

    if (deg == 1)
    {
        ca_div(roots, poly, poly + 1, ctx);
        ca_neg(roots, roots, ctx);
        return 1;
    }

    if (deg == 2)
    {
        ca_srcptr a, b, c;
        ca_t d, t;

        a = poly + 2;
        b = poly + 1;
        c = poly + 0;

        ca_init(d, ctx);
        ca_init(t, ctx);

        ca_mul(t, a, c, ctx);
        ca_mul_ui(t, t, 4, ctx);
        ca_sqr(d, b, ctx);
        ca_sub(d, d, t, ctx);
        ca_sqrt(d, d, ctx);

        ca_inv(t, a, ctx);
        ca_div_ui(t, t, 2, ctx);

        ca_sub(roots, d, b, ctx);
        ca_add(roots + 1, b, d, ctx);
        ca_neg(roots + 1, roots + 1, ctx);
        ca_mul(roots, roots, t, ctx);
        ca_mul(roots + 1, roots + 1, t, ctx);

        ca_clear(d, ctx);
        ca_clear(t, ctx);
        return 1;
    }

    if (_ca_vec_is_fmpq_vec(poly, len, ctx))
    {
        fmpz_poly_t f;
        qqbar_ptr r;
        slong i;

        f->coeffs = _fmpz_vec_init(len);
        f->length = f->alloc = len;
        r = qqbar_vec_init(len - 1);

        if (_ca_vec_fmpq_vec_is_fmpz_vec(poly, len, ctx))
        {
            for (i = 0; i < len; i++)
                fmpz_set(f->coeffs + i, CA_FMPQ_NUMREF(poly + i));
        }
        else
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_one(t);
            for (i = 0; i < len; i++)
                fmpz_lcm(t, t, CA_FMPQ_DENREF(poly + i));

            for (i = 0; i < len; i++)
            {
                fmpz_divexact(f->coeffs + i, t, CA_FMPQ_DENREF(poly + i));
                fmpz_mul(f->coeffs + i, f->coeffs + i, CA_FMPQ_NUMREF(poly + i));
            }

            fmpz_clear(t);
        }

        qqbar_roots_fmpz_poly(r, f, 0);

        for (i = 0; i < deg; i++)
            ca_set_qqbar(roots + i, r + i, ctx);

        _fmpz_vec_clear(f->coeffs, len);
        qqbar_vec_clear(r, len - 1);

        return 1;
    }

    return 0;
}

int
ca_poly_roots(ca_vec_t roots, const ca_poly_t poly, ulong flags, ca_ctx_t ctx)
{
    if (poly->length == 0)
        return 0;

    ca_vec_set_length(roots, poly->length - 1, ctx);

    if (_ca_poly_roots(roots->coeffs, poly->coeffs, poly->length, 0, ctx))
        return 1;

    return 0;
}
