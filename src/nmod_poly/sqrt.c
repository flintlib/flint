/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

static inline
int _nmod_poly_sqrt_2(nn_ptr s, nn_srcptr p, slong len)
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
_nmod_poly_sqrt(nn_ptr s, nn_srcptr p, slong len, nmod_t mod)
{
    slong slen;
    int result;
    nn_ptr t;
    ulong c, d;

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
        if (nmod_mul(c, c, mod) != d)
            flint_throw(FLINT_ERROR, "Exception (nmod_poly_sqrt). "
                "Modulus must be prime.\n");
    }

    if (len == 1)
    {
        s[0] = c;
        return 1;
    }

    slen = len / 2 + 1;

    /* Even modulus: _nmod_poly_sqrt_series requires prime modulus */
    if (mod.n % 2 == 0)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_sqrt). "
            "Modulus must be prime.\n");

    t = _nmod_vec_init(len);

    if (c == 1)
        _nmod_poly_sqrt_series(s, p, slen, slen, mod);
    else
    {
        ulong dinv, g;
        g = n_gcdinv(&dinv, d, mod.n);
        if (g != 1)
            flint_throw(FLINT_ERROR, "Exception (nmod_poly_sqrt). "
                "Modulus must be prime.\n");
       _nmod_vec_scalar_mul_nmod(t, p, slen, dinv, mod);
        _nmod_poly_sqrt_series(s, t, slen, slen, mod);
    }

    if (c != 1)
        _nmod_vec_scalar_mul_nmod(s, s, slen, c, mod);

    /* Full verification: s^2 must equal p (all coefficients).
       This catches incorrect results from Newton iteration
       when the modulus is not prime. */
    _nmod_poly_mul(t, s, slen, s, slen, mod);

    result = _nmod_vec_equal(t, p, len);
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
