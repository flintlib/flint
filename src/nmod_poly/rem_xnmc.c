/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

/* division by x**n - 1 */
void _nmod_poly_rem_xnm1(nn_ptr R, nn_srcptr A, slong len, ulong n, ulong modn)
{
    /* assumes len >= n */
    slong i;
    ulong j, r;

    r = len % n;
    i = len - r - n;  /* multiple of n, >= 0 by assumption */

    for (j = 0; j < r; j++)
        R[j] = n_addmod(A[i+n+j], A[i+j], modn);
    for (j = r; j < n; j++)
        R[j] = A[i+j];

    i -= n;
    while (i >= 0)
    {
        for (j = 0; j < n; j++)
            R[j] = n_addmod(R[j], A[i+j], modn);
        i -= n;
    }
}

/* division by x**n + 1 */
void _nmod_poly_rem_xnp1(nn_ptr R, nn_srcptr A, slong len, ulong n, ulong modn)
{
    /* assumes len >= n */
    slong i;
    ulong j, r;

    r = len % n;
    i = len - r - n;  /* multiple of n, >= 0 by assumption */

    for (j = 0; j < r; j++)
        R[j] = n_submod(A[i+j], A[i+n+j], modn);
    for (j = r; j < n; j++)
        R[j] = A[i+j];

    i -= n;
    while (i >= 0)
    {
        for (j = 0; j < n; j++)
            R[j] = n_submod(A[i+j], R[j], modn);
        i -= n;
    }
}

/* division by x**n - c, general variant */
void _nmod_poly_rem_xnmc(nn_ptr R, nn_srcptr A, slong len, ulong n, ulong c, nmod_t mod)
{
    /* assumes len >= n */
    slong i;
    ulong j, r, val;

    r = len % n;
    i = len - r - n;  /* multiple of n, >= 0 by assumption */

    for (j = 0; j < r; j++)
    {
        val = nmod_mul(A[i+n+j], c, mod);
        R[j] = n_addmod(val, A[i+j], mod.n);
    }
    for (j = r; j < n; j++)
        R[j] = A[i+j];

    i -= n;
    while (i >= 0)
    {
        for (j = 0; j < n; j++)
        {
            val = nmod_mul(R[j], c, mod);
            R[j] = n_addmod(val, A[i+j], mod.n);
        }
        i -= n;
    }
}

/* division by x**n - c, with precomputation on c */
/* constraint: modn < 2**(FLINT_BITS-1) */
void _nmod_poly_rem_xnmc_precomp(nn_ptr R, nn_srcptr A, slong len, ulong n, ulong c, ulong c_precomp, ulong modn)
{
    /* assumes len >= n */
    slong i;
    ulong j, r, val;

    r = len % n;
    i = len - r - n;  /* multiple of n, >= 0 by assumption */

    for (j = 0; j < r; j++)
    {
        val = n_mulmod_shoup(c, A[i+n+j], c_precomp, modn);
        R[j] = n_addmod(val, A[i+j], modn);
    }
    for (j = r; j < n; j++)
        R[j] = A[i+j];

    i -= n;
    while (i >= 0)
    {
        for (j = 0; j < n; j++)
        {
            val = n_mulmod_shoup(c, R[j], c_precomp, modn);
            R[j] = n_addmod(val, A[i+j], modn);
        }
        i -= n;
    }
}

/* division by x**n - c, lazy with precomputation */
/* constraint: max(A) + 2*modn <= 2**FLINT_BITS */
/* coeff bounds: in [0, max(A)] | out [0, max(A) + 2*modn) */
void _nmod_poly_rem_xnmc_precomp_lazy(nn_ptr R, nn_srcptr A, slong len, ulong n, ulong c, ulong c_precomp, ulong modn)
{
    /* assumes len >= n */
    slong i;
    ulong j, r, val, p_hi, p_lo;

    r = len % n;
    i = len - r - n;  /* multiple of n, >= 0 by assumption */

    for (j = 0; j < r; j++)
    {
        /* computes either val = (c*val mod n) or val = (c*val mod n) + n */
        val = A[i+n+j];
        umul_ppmm(p_hi, p_lo, c_precomp, val);
        val = c * val - p_hi * modn;
        /* lazy addition, yields R[i+j] in [0..k+2n), where max(R) <= k */
        R[j] = val + A[i+j];
    }
    for (j = r; j < n; j++)
        R[j] = A[i+j];

    i -= n;
    while (i >= 0)
    {
        for (j = 0; j < n; j++)
        {
            /* computes either val = (c*val mod n) or val = (c*val mod n) + n */
            val = R[j];
            umul_ppmm(p_hi, p_lo, c_precomp, val);
            val = c * val - p_hi * modn;
            /* lazy addition, yields R[i+j] in [0..k+2n), where max(R) <= k */
            R[j] = val + A[i+j];
        }
        i -= n;
    }
}

void nmod_poly_rem_xnmc(nmod_poly_t R, nmod_poly_t A, ulong n, ulong c)
{
    const ulong len = A->length;

    if (len <= n)
    {
        nmod_poly_set(R, A);
        return;
    }

    if (c == 0)
    {
        nmod_poly_set_trunc(R, A, n);
        return;
    }

    int lazy = 0;
    nn_ptr Rvec = _nmod_vec_init(n);

    /* perform division */
    if (c == 1)
        _nmod_poly_rem_xnm1(Rvec, A->coeffs, len, n, A->mod.n);

    else if (c == A->mod.n - 1)
        _nmod_poly_rem_xnp1(Rvec, A->coeffs, len, n, A->mod.n);

    /* if degree below the n_mulmod_shoup threshold, */
    /* or if modulus forbids n_mulmod_shoup usage, use general */
#if FLINT_MULMOD_SHOUP_THRESHOLD <= 2
    else if (A->mod.norm == 0)  /* here A->length >= threshold */
#else
    else if ((A->length < FLINT_MULMOD_SHOUP_THRESHOLD)
           || (A->mod.norm == 0))
#endif
    {
        _nmod_poly_rem_xnmc(Rvec, A->coeffs, len, n, c, A->mod);
    }

    else
    {
        const ulong modn = A->mod.n;
        const ulong c_precomp = n_mulmod_precomp_shoup(c, modn);

        /* if 3*mod.n - 1 <= 2**FLINT_BITS, use precomp+lazy variant */
#if FLINT_BITS == 64
        if (modn <= UWORD(6148914691236517205))
#else /* FLINT_BITS == 32 */
        if (modn <= UWORD(1431655765))
#endif
        {
            lazy = 1;
            _nmod_poly_rem_xnmc_precomp_lazy(Rvec, A->coeffs, len, n, c, c_precomp, modn);
        }

        /* use n_mulmod_shoup, non-lazy variant */
        else
        {
            _nmod_poly_rem_xnmc_precomp(Rvec, A->coeffs, len, n, c, c_precomp, modn);
        }
    }

    /* copy remainder R */
    nmod_poly_fit_length(R, n);
    if (lazy)
    {
        /* correct excess */
        const ulong modn = A->mod.n;
        for (ulong i = 0; i < n; i++)
        {
            ulong v = Rvec[i];
            if (v >= 2*modn)
                v -= 2*modn;
            else if (v >= modn)
                v -= modn;
            R->coeffs[i] = v;
        }
    }
    else
    {
        _nmod_vec_set(R->coeffs, Rvec, n);
    }
    _nmod_poly_set_length(R, n);
    _nmod_poly_normalise(R);

    _nmod_vec_clear(Rvec);
}
