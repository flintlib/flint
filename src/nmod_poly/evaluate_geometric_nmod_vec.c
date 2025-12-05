/*
    Copyright (C) 2025, Vincent Neiger, Éric Schost
    Copyright (C) 2025, Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_poly.h"
#include "nmod_vec.h"
#if FLINT_HAVE_FFT_SMALL
#  include "fft_small.h"
#endif

void _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(nn_ptr vs, nn_srcptr poly, slong plen,
                                                         const nmod_geometric_progression_t G, slong len,
                                                         nmod_t mod)
{
    FLINT_ASSERT(len <= G->len);

    if (plen == 0)
    {
        _nmod_vec_zero(vs, len);
        return;
    }

    /* val = valuation of poly */
    slong val;
    for (val = 0; val < plen; val++)
    {
        if (poly[val] != 0)
        {
            break;
        }
    }

    /* FIXME below are 3 different versions (one requires fft_small) */
    /** goal is to compute [rev(p) * G->f]_{plen - 1}^{len}  (that is, coeffs [plen - 1, plen - 1 + len))
     * if p has valuation val, define a = p / x**val of length alen = plen - val, and this becomes
     * [rev(a) * (G->f >> x**val)]_{alen - 1}^{len}  (i.e. coeffs [alen - 1, alen - 1 + len))
     **/

    const slong alen = plen - val;

    /* version 1:  (2025-12-04: fastest in small lengths, waiting for optimized middle product) */
    /** this uses a short product: write rev(a) = x**(alen-1) a(1/x), and write F = (G->f(x) >> x**val) rem x**(alen - 1 + len)
     * rev(a) * F has length L = 2 * alen - 2 + len, reverse it: we get a * rev(F),
     * we want its coefficients from L - 1 - (alen - 1) = alen - 1 + len - 1
     * down to, included, L - 1 - (alen - 1 + len - 1) = alen - 1
     **/
    if (2 * (plen - val) - 2 + len <= 192  || !FLINT_HAVE_FFT_SMALL)
    {
        slong blen = alen - 1 + len;
        nn_ptr a = _nmod_vec_init(alen);
        nn_ptr b = _nmod_vec_init(blen);

        for (slong i = val; i < plen; i++)
            a[i - val] = nmod_mul(G->x[i], poly[i], mod);

        nn_ptr Frev = _nmod_vec_init(blen);
        _nmod_poly_reverse(Frev, G->f->coeffs + val, blen, blen);
        _nmod_poly_mullow(b, Frev, blen, a, alen, blen, mod);

        for (slong i = 0; i < len; i++)
            vs[i] = nmod_mul(G->x[i], b[alen - 1 + len - 1 - i], mod);

        _nmod_vec_clear(Frev);
        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
    }
    else
    {
    /* version 2 */
    /* this uses a middle product to compute [rev(p) * G->f]_{plen - 1}^{len}  (i.e. coeffs [plen - 1, plen - 1 + len)) */
#if FLINT_HAVE_FFT_SMALL
        /* version 2.a uses fft_small directly  (2025-12-04: fastest in medium and large lengths, like 100 and more) */
        nn_ptr b = _nmod_vec_init(alen + len - 1);

        for (slong i = val; i < plen; i++)
            b[plen - 1 - i] = nmod_mul(G->x[i], poly[i], mod);

        _nmod_poly_mul_mid_default_mpn_ctx(b, alen - 1, alen - 1 + len, G->f->coeffs + val, alen - 1 + len, b, alen, mod);

        for (slong i = 0; i < len; i++)
            vs[i] = nmod_mul(G->x[i], b[i], mod);

        _nmod_vec_clear(b);
#else
        /* version 2.b uses nmod_poly_mulhigh  (2025-12-04: tested correct, but disabled: nmod_poly_mulhigh not yet optimized) */
        nn_ptr a = _nmod_vec_init(alen);
        nn_ptr b = _nmod_vec_init(alen + (alen - 1 + len));

        for (slong i = val; i < plen; i++)
            a[plen - 1 - i] = nmod_mul(G->x[i], poly[i], mod);

        _nmod_poly_mulhigh(b, G->f->coeffs + val, alen - 1 + len, a, alen, alen - 1, mod);

        for (slong i = 0; i < len; i++)
            vs[i] = nmod_mul(G->x[i], b[alen - 1 + i], mod);

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
#endif
    }
}

void _nmod_poly_evaluate_geometric_nmod_vec_fast(nn_ptr ys, nn_srcptr poly, slong plen, ulong r, slong n, nmod_t mod)
{
    if (n == 0)
        return;

    nmod_geometric_progression_t G;
    nmod_geometric_progression_init(G, r, FLINT_MAX(n, plen), mod);
    _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(ys, poly, plen, G, n, mod);
    nmod_geometric_progression_clear(G);
}

void nmod_poly_evaluate_geometric_nmod_vec_fast(nn_ptr ys, const nmod_poly_t poly, ulong r, slong n)
{
    _nmod_poly_evaluate_geometric_nmod_vec_fast(ys, poly->coeffs,
                                                poly->length, r, n, poly->mod);
}

void _nmod_poly_evaluate_geometric_nmod_vec_iter(nn_ptr ys, nn_srcptr coeffs, slong len, ulong r, slong n, nmod_t mod)
{
    slong i;
    ulong rpow = 1;

    ulong r2 = nmod_mul(r, r, mod);

    for (i = 0; i < n; i++)
    {
        ys[i] = _nmod_poly_evaluate_nmod(coeffs, len, rpow, mod);
        rpow = nmod_mul(rpow, r2, mod);
    }
}

void nmod_poly_evaluate_geometric_nmod_vec_iter(nn_ptr ys, const nmod_poly_t poly, ulong r, slong n)
{
    _nmod_poly_evaluate_geometric_nmod_vec_iter(ys, poly->coeffs,
                                        poly->length, r, n, poly->mod);
}
