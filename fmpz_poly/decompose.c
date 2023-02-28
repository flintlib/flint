/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz_poly.h"
#include "fmpq.h"
#include "long_extras.h"

#define FMPZ_PTR_SWAP(a, b) \
    do {                    \
        fmpz * __tt__ = b;  \
        b = a;              \
        a = __tt__;         \
    } while (0)

/* try f(x) = g(h(x)), deg g = r, deg h = s, Kozen-Landau */
int _fmpz_poly_decompose(
    fmpz * g,           /* length r + 1 */
    slong r,
    fmpz * h,           /* length s + 1 */
    slong s,
    const fmpz * f)    /* length r*s + 1 */
{
    int success = 1;
    slong hshift, k, j, i;
    fmpz_t cnum, cden, gg, hden, tden, hdenbar;
    fmpz * t, * tq, * tr;

    FLINT_ASSERT(g != h && g != f && h != f);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(s > 1);

    t = _fmpz_vec_init(r*s + 1);
    tq = _fmpz_vec_init(r*s + 1);
    tr = _fmpz_vec_init(r*s + 1);

    fmpz_init(cnum);
    fmpz_init(cden);
    fmpz_init(gg);
    fmpz_init(hden);
    fmpz_init(tden);
    fmpz_init(hdenbar);

    /* h = x^s, hshift is highest power of x dividing h */
    k = s; 
    hshift = s;
    fmpz_one(h + hshift);
    fmpz_one(hden);
    for (i = 0; i < hshift; i++)
        fmpz_zero(h + i);

    for (k = 1; k < s; k++)
    {
        _fmpz_poly_pow(t, h + hshift, 1 + s - hshift, r);
        fmpz_pow_ui(tden, hden, r);

        if (fmpz_sgn(f + r*s) < 0)
        {
            fmpz_neg(cnum, f + r*s - k);
            fmpz_neg(cden, f + r*s);
        }
        else
        {
            fmpz_set(cnum, f + r*s - k);
            fmpz_set(cden, f + r*s);
        }

        if (k <= r*(s - hshift))
        {
            _fmpq_sub(cnum, cden, cnum, cden, t + r*(s - hshift) - k, tden);
        }

        fmpz_mul_ui(cden, cden, r);


        if (fmpz_is_zero(cnum))
            continue;

        FLINT_ASSERT(fmpz_sgn(cden) > 0);
        fmpz_gcd(gg, cnum, cden);
        fmpz_divexact(cnum, cnum, gg);
        fmpz_divexact(cden, cden, gg);

        fmpz_gcd(gg, hden, cden);
        fmpz_divexact(hdenbar, hden, gg);
        fmpz_divexact(cden, cden, gg);

        for (i = hshift; i <= s; i++)
            fmpz_mul(h + i, h + i, cden);

        hshift = FLINT_MIN(hshift, s - k);
        fmpz_mul(h + s - k, cnum, hdenbar);
        fmpz_mul(hden, hden, cden);
    }

    j = 0;

    fmpz_set(g + 0, f + 0);

    for (k = 1; k < hshift; k++)
        if (!fmpz_is_zero(f + k))
            goto fail;

    if (!_fmpz_poly_divrem(t, tr, f + hshift, r*s + 1 - s*j - hshift,
                                             h + hshift, s + 1 - hshift, 1))
        goto fail;

    for (k = 0; k < r*s + 1 - s*j - hshift; k++)
        if (!fmpz_is_zero(tr + k))
            goto fail;

    for (j = 1; j < r; j++)
    {
        /* t has degree r*s - s*j */
        fmpz_swap(g + j, t + 0);

        for (k = 1; k < hshift; k++)
            if (!fmpz_is_zero(t + k))
                goto fail;

        if (!_fmpz_poly_divrem(tq, tr, t + hshift, r*s + 1 - hshift - s*j,
                                       h + hshift, s + 1 - hshift, 1))
            goto fail;

        for (k = 0; k < r*s + 1 - hshift - s*j; k++)
            if (!fmpz_is_zero(tr + k))
                goto fail;

        FMPZ_PTR_SWAP(t, tq);
    }

    /* t has degree 0 */
    fmpz_swap(g + r, t + 0);

cleanup:

    _fmpz_vec_clear(t, r*s + 1);
    _fmpz_vec_clear(tq, r*s + 1);
    _fmpz_vec_clear(tr, r*s + 1);

    fmpz_clear(cnum);
    fmpz_clear(cden);
    fmpz_clear(gg);
    fmpz_clear(hden);
    fmpz_clear(tden);
    fmpz_clear(hdenbar);

    return success;

fail:

    success = 0;
    goto cleanup;
}


int fmpz_poly_decompose(fmpz_poly_t g, slong r, fmpz_poly_t h, slong s,
                                                           const fmpz_poly_t f)
{
    int success;
    slong rs;

    if (z_mul_checked(&rs, r, s) || r <= 1 || s <= 1 ||
        rs != fmpz_poly_degree(f))
    {
        return 0;
    }

    if (g != f && h != f)
    {
        fmpz_poly_fit_length(g, r + 1);
        fmpz_poly_fit_length(h, s + 1);
        success = _fmpz_poly_decompose(g->coeffs, r, h->coeffs, s, f->coeffs);
    }
    else
    {
        fmpz_poly_t gt, ht;
        fmpz_poly_init2(gt, r + 1);
        fmpz_poly_init2(ht, s + 1);
        success = _fmpz_poly_decompose(gt->coeffs, r, ht->coeffs, s, f->coeffs);
        fmpz_poly_swap(g, gt);
        fmpz_poly_swap(h, ht);
        fmpz_poly_clear(gt);
        fmpz_poly_clear(ht);
    }

    _fmpz_poly_set_length(g, r + 1);
    _fmpz_poly_set_length(h, s + 1);

    /* retain valid (but junk) outputs in case of failure */
    _fmpz_poly_normalise(g);
    _fmpz_poly_normalise(h);

    return success;
}
