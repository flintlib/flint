/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_poly.h"

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
        val = n_mulmod2_preinv(val, c, mod.n, mod.ninv);
        val = n_addmod(val, poly[m], mod.n);
    }

    return val;
}

ulong
_nmod_poly_evaluate_nmod_shoup(nn_srcptr poly, slong len, ulong c, ulong c_precomp, nmod_t mod)
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
_nmod_poly_evaluate_nmod_shoup_lazy(nn_srcptr poly, slong len, ulong c, ulong c_precomp, nmod_t mod)
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
        // see FIXME n_mulmod_shoup explanations for details
        umul_ppmm(p_hi, p_lo, c_precomp, val);
        val = c * val - p_hi * mod.n;
        // lazy addition in [0..3n-1) --> requires 3n-2 < 2**FLINT_BITS
        val += poly[m];
    }

    // correct excess
    if (val >= 2*mod.n)
        val -= 2*mod.n;
    if (val >= mod.n)
        val -= mod.n;

    return val;
}

ulong
nmod_poly_evaluate_nmod(const nmod_poly_t poly, ulong c)
{
    // if degree below the n_mulmod_shoup threshold (FIXME issue #2059)
    // or modulus forbids n_mulmod_shoup usage, use straightforward code
    if ((poly->length <= 10)         
           || (poly->mod.norm == 0)) // 
    {
        return _nmod_poly_evaluate_nmod(poly->coeffs, poly->length, c, poly->mod);
    }

    // if 3*mod.n - 2 < 2**FLINT_BITS, apply Shoup's precomputation, lazy variant
    else
#if FLINT_BITS == 64
    if (poly->mod.n <= UWORD(6148914691236517205))
#else // FLINT_BITS == 32
    if (poly->mod.n <= UWORD(1431655765))
#endif
    {
        const ulong c_precomp = n_mulmod_precomp_shoup(c, poly->mod.n);
        return _nmod_poly_evaluate_nmod_shoup_lazy(poly->coeffs, poly->length, c, c_precomp, poly->mod);
    }

    // apply Shoup's precomputation, non-lazy variant
    else
    {
        const ulong c_precomp = n_mulmod_precomp_shoup(c, poly->mod.n);
        return _nmod_poly_evaluate_nmod_shoup(poly->coeffs, poly->length, c, c_precomp, poly->mod);
    }
}
