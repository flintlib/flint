/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

static inline
int _nmod_poly_sqrt_2(mp_ptr s, mp_srcptr p, slong len)
{
   slong i;

   for (i = 1; i < len; i += 2)
       if (p[i] != 0)
           return 0;

   for (i = 0; i < len; i += 2)
       s[i / 2] = p[i];

   return 1;
}

int
_nmod_poly_sqrt(mp_ptr s, mp_srcptr p, slong len, nmod_t mod)
{
    slong slen;
    int result;
    mp_ptr t;
    mp_limb_t c, d;

    if (len % 2 == 0)
        return len == 0;

    if (mod.n == 2)
        return _nmod_poly_sqrt_2(s, p, len);

    /* valuation must be even, and then can be reduced to 0 */
    while (p[0] == 0)
    {
        if (p[1] != 0)
            return 0;

        s[0] = 0;
        p += 2;
        len -= 2;
        s++;
    }

    c = d = p[0];
    if (c != 1)
    {
        c = n_sqrtmod(c, mod.n);
        if (c == 0)
            return 0;
    }

    if (len == 1)
    {
        s[0] = c;
        return 1;
    }

    slen = len / 2 + 1;

    t = _nmod_vec_init(len);

    if (c == 1)
        _nmod_poly_sqrt_series(s, p, slen, slen, mod);
    else
    {
       _nmod_vec_scalar_mul_nmod(t, p, slen, n_invmod(d, mod.n), mod);
        _nmod_poly_sqrt_series(s, t, slen, slen, mod);
    }

    if (c != 1)
        _nmod_vec_scalar_mul_nmod(s, s, slen, c, mod);

    _nmod_poly_mulhigh(t, s, slen, s, slen, slen, mod);


    result = _nmod_vec_equal(t + slen, p + slen, len - slen);
    _nmod_vec_clear(t);
    return result;
}

int
nmod_poly_sqrt(nmod_poly_t b, const nmod_poly_t a)
{
    slong blen, len = a->length;
    int result;

    if (len % 2 == 0)
    {
        nmod_poly_zero(b);
        return len == 0;
    }

    if (b == a)
    {
        nmod_poly_t tmp;
        nmod_poly_init_preinv(tmp, a->mod.n, a->mod.ninv);
        result = nmod_poly_sqrt(tmp, a);
        nmod_poly_swap(b, tmp);
        nmod_poly_clear(tmp);
        return result;
    }

    blen = len / 2 + 1;
    nmod_poly_fit_length(b, blen);
    b->length = blen;
    result = _nmod_poly_sqrt(b->coeffs, a->coeffs, len, a->mod);
    if (!result)
        b->length = 0;
    return result;
}
