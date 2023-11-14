/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

static inline
int _fmpz_mod_poly_sqrt_2(fmpz * s, const fmpz * p, slong len)
{
   slong i;

   for (i = 1; i < len; i += 2)
       if (!fmpz_is_zero(p + i))
           return 0;

   for (i = 0; i < len; i += 2)
       fmpz_set(s + i / 2,  p + i);

   return 1;
}

int
_fmpz_mod_poly_sqrt(fmpz * s, const fmpz * p, slong len, const fmpz_mod_ctx_t mod)
{
    slong slen, i;
    int result;
    fmpz * t;
    fmpz_t c, d;

    if (len % 2 == 0)
        return len == 0;

    if (fmpz_cmp_ui(mod->n, 2) == 0)
        return _fmpz_mod_poly_sqrt_2(s, p, len);

    /* valuation must be even, and then can be reduced to 0 */
    while (fmpz_is_zero(p + 0))
    {
        if (!fmpz_is_zero(p + 1))
            return 0;

        fmpz_zero(s + 0);
        p += 2;
        len -= 2;
        s++;
    }

    fmpz_init(c);
    fmpz_init(d);

    fmpz_set(d, p + 0);
    fmpz_set(c, d);

    if (!fmpz_is_one(c))
    {
        if (!fmpz_sqrtmod(c, c, mod->n))
        {
           result = 0;
           goto cleanup;
        }
    }

    if (len == 1)
    {
        fmpz_set(s + 0, c);

        result = 1;
        goto cleanup;
    }

    slen = len / 2 + 1;

    t = _fmpz_vec_init(len);

    if (fmpz_is_one(c))
        _fmpz_mod_poly_sqrt_series(s, p, slen, slen, mod);
    else
    {
        fmpz_invmod(d, d, mod->n);

        _fmpz_mod_vec_scalar_mul_fmpz_mod(t, p, slen, d, mod);
        _fmpz_mod_poly_sqrt_series(s, t, slen, slen, mod);
    }

    if (!fmpz_is_one(c))
        _fmpz_mod_vec_scalar_mul_fmpz_mod(s, s, slen, c, mod);

    _fmpz_poly_mulhigh(t, s, slen, s, slen, slen);

    for (i = 0; i < slen; i++)
        fmpz_zero(t + i);

    _fmpz_mod_vec_set_fmpz_vec(t + slen, t + slen, slen - 1, mod);

    result = _fmpz_vec_equal(t + slen, p + slen, len - slen);

    _fmpz_vec_clear(t, len);

cleanup:

    fmpz_clear(c);
    fmpz_clear(d);

    return result;
}

int
fmpz_mod_poly_sqrt(fmpz_mod_poly_t b, const fmpz_mod_poly_t a, const fmpz_mod_ctx_t ctx)
{
    slong blen, len = a->length;
    int result;

    if (len % 2 == 0)
    {
        fmpz_mod_poly_zero(b, ctx);

        return len == 0;
    }

    if (b == a)
    {
        fmpz_mod_poly_t tmp;

        fmpz_mod_poly_init(tmp, ctx);

        result = fmpz_mod_poly_sqrt(tmp, a, ctx);
        fmpz_mod_poly_swap(b, tmp, ctx);

        fmpz_mod_poly_clear(tmp, ctx);

        return result;
    }

    blen = len / 2 + 1;
    fmpz_mod_poly_fit_length(b, blen, ctx);

    result = _fmpz_mod_poly_sqrt(b->coeffs, a->coeffs, len, ctx);

    if (!result)
        blen = 0;

    _fmpz_mod_poly_set_length(b, blen);
    _fmpz_mod_poly_normalise(b);

    return result;
}
