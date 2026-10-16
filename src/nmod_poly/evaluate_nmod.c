/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2024 Vincent Neiger
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod.h"

static void
_nmod_vec_set_powers(nn_ptr res, ulong x, slong len, nmod_t mod)
{
    slong i;
    res[0] = 1;
    res[1] = x;

    for (i = 2; i < len; i++)
        res[i] = nmod_mul(res[i - 1], x, mod);
}

static void
_nmod_vec_set_powers_precomp(nn_ptr res, ulong x, slong len, nmod_t mod)
{
    slong i;
    res[0] = 1;
    res[1] = x;
    ulong xpre = n_mulmod_precomp_shoup(x, mod.n);

    for (i = 2; i < len; i++)
        res[i] = n_mulmod_shoup(x, res[i - 1], xpre, mod.n);
}

ulong
_nmod_poly_evaluate_nmod_rectangular(nn_srcptr poly,
    slong len, ulong x, nmod_t mod)
{
    slong i, m, r;
    nn_ptr xs;
    ulong s, y, xmpre;
    TMP_INIT;

    if (len == 0)
        return 0;
    if (len == 1 || x == 0)
        return poly[0];
    if (len == 2)
        return nmod_addmul(poly[0], poly[1], x, mod);

    m = n_sqrt(len) + 1;
    m = FLINT_MIN(m, 1024);
    r = (len + m - 1) / m;

    dot_params_t params = _nmod_vec_dot_params(m - 1, mod);

    TMP_START;
    xs = TMP_ALLOC((m + 1) * sizeof(ulong));

    if (NMOD_CAN_USE_SHOUP(mod))
    {
        _nmod_vec_set_powers_precomp(xs, x, m + 1, mod);
        xmpre = n_mulmod_precomp_shoup(xs[m], mod.n);

        y = _nmod_vec_dot(poly + (r - 1) * m + 1, xs + 1, len - (r - 1) * m - 1, mod, params);
        y = nmod_add(y, poly[(r - 1) * m], mod);

        for (i = r - 2; i >= 0; i--)
        {
            s = _nmod_vec_dot(poly + i * m + 1, xs + 1, m - 1, mod, params);
            s = nmod_add(s, poly[i * m], mod);
            y = n_mulmod_shoup(xs[m], y, xmpre, mod.n);
            y = nmod_add(y, s, mod);
        }
    }
    else
    {
        _nmod_vec_set_powers(xs, x, m + 1, mod);
        y = _nmod_vec_dot(poly + (r - 1) * m + 1, xs + 1, len - (r - 1) * m - 1, mod, params);
        y = nmod_add(y, poly[(r - 1) * m], mod);

        for (i = r - 2; i >= 0; i--)
        {
            s = _nmod_vec_dot(poly + i * m + 1, xs + 1, m - 1, mod, params);
            s = nmod_add(s, poly[i * m], mod);
            y = nmod_mul(xs[m], y, mod);
            y = nmod_add(y, s, mod);
        }
    }

    TMP_END;

    return y;
}

ulong
_nmod_poly_evaluate_nmod_horner(nn_srcptr poly, slong len, ulong c, nmod_t mod)
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
_nmod_poly_evaluate_nmod_precomp(nn_srcptr poly, slong len, ulong c, ulong c_precomp, ulong modn)
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
        val = n_mulmod_shoup(c, val, c_precomp, modn);
        val = n_addmod(val, poly[m], modn);
    }

    return val;
}

ulong
_nmod_poly_evaluate_nmod_precomp_lazy(nn_srcptr poly, slong len, ulong c, ulong c_precomp, ulong modn)
{
    slong m;
    ulong val, p_hi;

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
        p_hi = n_mulhi(c_precomp, val);
        val = c * val - p_hi * modn;
        // lazy addition, yields val in [0..k+2n-1), where max(poly) < k
        // --> if k == n (poly is reduced mod n), constraint: 3n-1 <= 2**FLINT_BITS
        val += poly[m];
    }

    return val;
}

/* if 3*mod.n - 1 <= 2**FLINT_BITS, can use the lazy variant */
#if FLINT_BITS == 64
#define LAZY_MAX UWORD(6148914691236517205)
#else
#define LAZY_MAX UWORD(1431655765)
#endif

ulong
_nmod_poly_evaluate_nmod(nn_srcptr poly, slong len, ulong c, nmod_t mod)
{
    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    if (!NMOD_CAN_USE_SHOUP(mod))
    {
        if (len <= 13)
            return _nmod_poly_evaluate_nmod_horner(poly, len, c, mod);
        else
            return _nmod_poly_evaluate_nmod_rectangular(poly, len, c, mod);
    }

    if (len < FLINT_MULMOD_SHOUP_THRESHOLD)
        return _nmod_poly_evaluate_nmod_horner(poly, len, c, mod);

    if (len > 50)
        return _nmod_poly_evaluate_nmod_rectangular(poly, len, c, mod);

    const ulong modn = mod.n;

    if (modn <= LAZY_MAX)
    {
        const ulong c_precomp = n_mulmod_precomp_shoup(c, modn);
        ulong val = _nmod_poly_evaluate_nmod_precomp_lazy(poly, len, c, c_precomp, modn);
        /* Correct excess. */
        if (val >= 2*modn)
            val -= 2*modn;
        else if (val >= modn)
            val -= modn;
        return val;
    }
    else
    {
        const ulong c_precomp = n_mulmod_precomp_shoup(c, modn);
        return _nmod_poly_evaluate_nmod_precomp(poly, len, c, c_precomp, modn);
    }
}

ulong
nmod_poly_evaluate_nmod(const nmod_poly_t poly, ulong c)
{
    return _nmod_poly_evaluate_nmod(poly->coeffs, poly->length, c, poly->mod);
}

