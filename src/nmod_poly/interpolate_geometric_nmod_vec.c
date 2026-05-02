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

void _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nn_ptr poly,
            nn_srcptr v, const nmod_geometric_progression_t G, slong len, nmod_t mod)
{
    FLINT_ASSERT(len <= G->len);
    FLINT_ASSERT((G->function >> 1) & 1);

    slong i, N, f_len, h_len;
    nn_ptr f, h;

    N = G->len;

    if (len == 1)
    {
        poly[0] = v[0];
        return;
    }

    f = _nmod_vec_init(N);
    h = _nmod_vec_init(N);

    /** step1: Newton interpolation
     * [Bostan - Schost, J.Complexity 2005, Section 5.1]
     * -> The coefficients of the interpolant, in the Newton basis associated
     *  to the geometric progression 1, q, q**2, q**3, etc., are obtained as
     *  c_0 / q_0, ..., c_{len-1} / q_{len-1}, where c_0, ..., c_{len-1} are
     *  the first coefficients of the product
     *     (sum_{i=0}^{len-1} v[i]/u_i x**i) (sum_{i=0}^{len-1} (-1)**i q_i/u_i x**i)
     *  where v[i] is the element at index `i` in the input values `v`,
     *  and u_i = prod_{1 <= k <= i} (q**k - 1),
     *  and q_i = q**(i*(i-1)/2)
     * -> With the precomputed data, these are the `len` coefficients of
     *    f * G->int_f1  mod x**len
     * where f = sum_{i=0}^{len-1} v[i] * G->int_s1[i] x**i
     */

    /* val = valuation of output poly in Newton basis */
    slong val = 0;
    for (val = 0; val < len; val++)
        if (v[val] != 0)
            break;
    if (val == len)
    {
        _nmod_vec_zero(poly, len);
        return;
    }

    /* i = number of trailing zero coefficients */
    for (f_len = len; f_len > val; f_len--)
        if (v[f_len-1] != 0)
            break;
    f_len = f_len - val;

    /* f = sum_{i=val}^{len-1} v[i] * G->int_s1[i] x**{i-val} */
    /*   == sum_{i=0}^{f_len-1} v[i+val] * G->int_s1[i+val] x**i */
    for (i = 0; i < f_len; i++)
        f[i] = nmod_mul(v[i+val], G->int_s1[i+val], mod);

    /* h = (x**val * f) * G->int_f1  mod x**len                     */
    /*   == x**val (f * G->int_f1  mod x**(len-val))                */
    /* note: len - val is <= G->int_f1->length, since G->intf_1 has */
    /* length G->len >= len (all its coefficients are nonzero)      */
    _nmod_vec_zero(h, val);
    _nmod_poly_mullow(h+val, G->int_f1->coeffs, len - val, f, f_len, len - val, mod);

    for (h_len = len; h_len > val; h_len--)
        if (h[h_len-1] != 0)
            break;

    /* FIXME division by q_i ?? */

    /* TODO remove debug */
    /* if (f_len != len) */
    /*     flint_printf("len = %ld, G->len = %ld, flen = %ld\n", len, G->len, f_len); */

    /** step2: Newton basis -> monomial basis             
     * [Bostan - Schost, J.Complexity 2005, Section 5.2]
     * in the Newton basis, poly has coefficients h[i] / q_i, i=0..len-1
     * -> Convert it to monomial basis, through the transposed product
     *       mul_t(len-1, sum_{i=0}^{len-1} uu_i x**i,
     *                    sum_{i=0}^{len-1} (-1)**i * (h[i]/q_i)*q_i/uu_i x**i)
     *  where q_i = q**(i*(i-1)/2) as above,
     *  and uu_i =  prod_{1 <= k <= i} q**k / (1 - q**k)
     *           == (-1)**i * q_i / u_i, for u_i as above
     *  This gives `len` coefficients and it just remains to scale the i-th one
     *  by (-1)**i * uu_i / q_i == 1 / u_i
     * -> so, with the precomputed data, the two polynomials have coefficients
     *     uu_i == G->int_f1[i] and (-1)**i * h[i]/uu_i == (-1)**i * h[i] * G->int_s2[i]
     *  and the final scaling factors are 1/u_i == G->int_s1[i]
     *   meaning that we are looking to compute
     *          mul_t(len-1, G->int_f1, f)
     *   where f = sum_{i=0}^{len-1} h[i] * G->int_s2[i] x**i
     * -> from that, it remains to scale coefficients by 1/u_i == G->int_s1[i]
     */

    for (val = 0; val < h_len; val++)
        if (h[val] != 0)
            break;

    for (i = val; i < h_len; i++)
        f[N - 1 - i] = nmod_mul(h[i], G->int_s2[i], mod);

    f_len = N - val;
    _nmod_vec_zero(f, N - h_len);
    _nmod_poly_mullow(h, G->int_f2->coeffs, G->int_f2->length, f, f_len, N, mod);

    for (i = 0; i < len; i++)
        poly[i] = nmod_mul(h[N - 1 - i], G->int_s3[i], mod);

    _nmod_vec_clear(f);
    _nmod_vec_clear(h);
}

void nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nmod_poly_t poly,
                nn_srcptr v, const nmod_geometric_progression_t G, slong len)
{
    FLINT_ASSERT((G->function >> 1) & 1);

    nmod_poly_fit_length(poly, len);
    _nmod_poly_set_length(poly, len);
    _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly->coeffs, v, G, len, G->mod);
    _nmod_poly_normalise(poly);
}

void nmod_poly_interpolate_geometric_nmod_vec_fast(nmod_poly_t poly,
                                                   ulong r, nn_srcptr ys, slong n)
{
    nmod_geometric_progression_t G;
    _nmod_geometric_progression_init_function(G, r, n, poly->mod, UWORD(2));
    nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly, ys, G, n);
    _nmod_geometric_progression_clear_function(G, UWORD(2));
}
