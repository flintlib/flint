/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod.h"

ulong
_nmod_poly_evaluate_nmod(nn_srcptr poly, slong len, ulong c, nmod_t mod)
{
    slong m;
    ulong val;

    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    m = len - 1;

    val = poly[m];
    m--;

    for ( ; m >= 0; m--)
    {
        val = nmod_mul(val, c, mod);
        val = n_addmod(val, poly[m], mod.n);
    }

    return val;
}

ulong
_nmod_poly_evaluate_nmod_precomp(nn_srcptr poly, slong len, ulong c, ulong c_precomp, nmod_t mod)
{
    slong m;
    ulong val;

    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    m = len - 1;

    val = poly[m];
    m--;

    for ( ; m >= 0; m--)
    {
        val = n_mulmod_shoup(c, val, c_precomp, mod.n);
        val = n_addmod(val, poly[m], mod.n);
    }

    return val;
}

ulong
_nmod_poly_evaluate_nmod_precomp_lazy(nn_srcptr poly, slong len, ulong c, ulong c_precomp, nmod_t mod)
{
    slong m;
    ulong val, p_hi, p_lo;

    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    m = len - 1;

    val = poly[m];
    m--;

    for ( ; m >= 0; m--)
    {
        // computes either val = (c*val mod n) or val = (c*val mod n) + n
        // see documentation of ulong_extras / n_mulmod_shoup for details
        umul_ppmm(p_hi, p_lo, c_precomp, val);
        val = c * val - p_hi * mod.n;
        // lazy addition, yields val in [0..3n-1)
        // --> requires 3n-2 < 2**FLINT_BITS
        val += poly[m];
    }

    // correct excess
    if (val >= 2*mod.n)
        return val - 2*mod.n;
    if (val >= mod.n)
        return val - mod.n;
    return val;
}

ulong
nmod_poly_evaluate_nmod(const nmod_poly_t poly, ulong c)
{
    if (poly->length == 0)
        return 0;

    if (poly->length == 1 || c == 0)
        return poly->coeffs[0];

    // if degree below the n_mulmod_shoup threshold
    // or modulus forbids n_mulmod_shoup usage, use nmod_mul
#if FLINT_MULMOD_SHOUP_THRESHOLD <= 2
    if (poly->mod.norm == 0)  // here poly->length >= threshold
#else
    if ((poly->length < FLINT_MULMOD_SHOUP_THRESHOLD)
           || (poly->mod.norm == 0))
#endif
    {
        return _nmod_poly_evaluate_nmod(poly->coeffs, poly->length, c, poly->mod);
    }

    // if 3*mod.n - 2 < 2**FLINT_BITS, use n_mulmod_shoup, lazy variant
    else
#if FLINT_BITS == 64
    if (poly->mod.n <= UWORD(6148914691236517205))
#else // FLINT_BITS == 32
    if (poly->mod.n <= UWORD(1431655765))
#endif
    {
        const ulong c_precomp = n_mulmod_precomp_shoup(c, poly->mod.n);
        return _nmod_poly_evaluate_nmod_precomp_lazy(poly->coeffs, poly->length, c, c_precomp, poly->mod);
    }

    // use n_mulmod_shoup, non-lazy variant
    else
    {
        const ulong c_precomp = n_mulmod_precomp_shoup(c, poly->mod.n);
        return _nmod_poly_evaluate_nmod_precomp(poly->coeffs, poly->length, c, c_precomp, poly->mod);
    }
}
