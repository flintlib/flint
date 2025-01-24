/*
    Copyright (C) 2007 Marco Bodrato
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/*
   Toom33 (interpolation in 5 points) using Bodrato scheme
   http://marco.bodrato.it/papers/Bodrato2007-OptimalToomCookMultiplicationForBinaryFieldAndIntegers.pdf

   Assumes commutativity, division by 3.
   Todo: squaring version.
   Todo: skip unnecessary zero-extensions of vectors and tighten
         allocations.
*/
int
_gr_poly_mul_toom33(gr_ptr res, gr_srcptr f, slong flen, gr_srcptr g, slong glen, gr_ctx_t ctx)
{
    gr_srcptr U0, U1, U2, V0, V1, V2;
    gr_ptr tmp, W0, W1, W2, W3, W4;
    slong m, U2len, V2len, U1len, V1len, U0len, V0len, rlen, len;
    slong W4len;
    slong sz = ctx->sizeof_elem;
    slong alloc;
    int status = GR_SUCCESS;

    /* TODO: should explicitly call basecase mul. */
    if (flen <= 1 || glen <= 1)
        return _gr_poly_mullow_generic(res, f, flen, g, glen, flen + glen - 1, ctx);

    /* U = U2*x^(2m) + U1*x^m + U0 */
    /* V = V2*x^(2m) + V1*x^m + V0 */
    /* Each block has length m */
    m = FLINT_MAX(flen, glen);
    m = (m + 3 - 1) / 3;
    U0 = f;
    U1 = GR_ENTRY(f, m, sz);
    U2 = GR_ENTRY(f, 2 * m, sz);
    V0 = g;
    V1 = GR_ENTRY(g, m, sz);
    V2 = GR_ENTRY(g, 2 * m, sz);

    U2len = FLINT_MAX(flen - 2 * m, 0);
    V2len = FLINT_MAX(glen - 2 * m, 0);
    U1len = FLINT_MIN(FLINT_MAX(flen - m, 0), m);
    V1len = FLINT_MIN(FLINT_MAX(glen - m, 0), m);
    U0len = FLINT_MIN(flen, m);
    V0len = FLINT_MIN(glen, m);
 
    alloc = 10 * m;
    GR_TMP_INIT_VEC(tmp, alloc, ctx);
    W0 = tmp;
    W1 = GR_ENTRY(W0, 2 * m, sz);
    W2 = GR_ENTRY(W1, 2 * m, sz);
    W3 = GR_ENTRY(W2, 2 * m, sz);
    W4 = GR_ENTRY(W3, 2 * m, sz);

    /* Evaluation: 5*2 add, 2 shift; 5mul */
    /* W0 = U2 + U0 */
    /* if max(U2len,U0len) < m, assumes top coefficients are already zeroed from the initialization */
    status |= _gr_poly_add(W0, U2, U2len, U0, U0len, ctx);
    /* W4 = V2 + V0 */
    /* if max(V2len,V0len) < m, assumes top coefficients are already zeroed from the initialization */
    status |= _gr_poly_add(W4, V2, V2len, V0, V0len, ctx);
    /* W2 = W0 - U1 */
    status |= _gr_poly_sub(W2, W0, m, U1, U1len, ctx);
    /* W1 = W4 - V1 */
    status |= _gr_poly_sub(W1, W4, m, V1, V1len, ctx);
    /* W0 = W0 + U1 */
    status |= _gr_poly_add(W0, W0, m, U1, U1len, ctx);
    /* W4 = W4 + V1 */
    status |= _gr_poly_add(W4, W4, m, V1, V1len, ctx);
    /* W3 = W2 * W1 */
    status |= _gr_poly_mul(W3, W2, m, W1, m, ctx);
    /* W1 = W0 * W4 */
    status |= _gr_poly_mul(W1, W0, m, W4, m, ctx);
    /* W0 = ((W0 + U2) << 1) - U0 */
    status |= _gr_poly_add(W0, W0, m, U2, U2len, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(W0, W0, m, 1, ctx);
    status |= _gr_poly_sub(W0, W0, m, U0, U0len, ctx);
    /* W4 = ((W4 + V2) << 1) - V0 */
    status |= _gr_poly_add(W4, W4, m, V2, V2len, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(W4, W4, m, 1, ctx);
    status |= _gr_poly_sub(W4, W4, m, V0, V0len, ctx);
    /* W2 = W0 * W4 */
    status |= _gr_poly_mul(W2, W0, m, W4, m, ctx);
    /* W0 = U0 * V0 */
    if (U0len > 0 && V0len > 0)
    {
        status |= _gr_poly_mul(W0, U0, U0len, V0, V0len, ctx);
        status |= _gr_vec_zero(GR_ENTRY(W0, U0len + V0len - 1, sz), 2 * m - (U0len + V0len - 1), ctx);
    }
    else
        status |= _gr_vec_zero(W0, 2 * m, ctx);
    /* W4 = U2 * V2 */
    /* We compute this length accurately instead of zero-extending. */
    if (U2len > 0 && V2len > 0)
    {
        W4len = U2len + V2len - 1;
        status |= _gr_poly_mul(W4, U2, U2len, V2, V2len, ctx);
    }
    else
    {
        W4len = 0;
    }

    /* toom42 variant */
    /* U = U3*x^(3m) + U2*x^(2m) + U1*x^m + U0 */
    /* V = V1*x^m + V0 */
    /* Evaluation: 7+3 add, 3 shift; 5mul */
    /*
    W0 = U1 + U3;
    W4 = U0 + U2;
    W3 = W4 + W0;
    W4 = W4 - W0;
    W0 = V0 + V1;
    W2 = V0 - V1;
    W1 = W3 * W0;
    W3 = W4 * W2;
    W4 = (((((U3<<1) + U2) << 1) + U1) << 1) + U0;
    W0 = W0 + V1;
    W2 = W4 * W0;
    W0 = U0 * V0;
    W4 = U3 * V1;
    */

    /* Interpolation: 8 add, 3 shift, 1 Sdiv */
    len = 2 * m - 1;
    /* W2 = (W2 - W3) / 3 */
    status |= _gr_vec_sub(W2, W2, W3, len, ctx);
    status |= _gr_vec_divexact_scalar_ui(W2, W2, len, 3, ctx);
    /* W3 = (W1 - W3) >> 1 */
    status |= _gr_vec_sub(W3, W1, W3, len, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(W3, W3, len, -1, ctx);
    /* W1 = W1 - W0 */
    status |= _gr_vec_sub(W1, W1, W0, len, ctx);
    /* W2 = ((W2 - W1) >> 1) - (W4 << 1) */
    status |= _gr_vec_sub(W2, W2, W1, len, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(W2, W2, len, -1, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(res, W4, W4len, 1, ctx);
    status |= _gr_vec_sub(W2, W2, res, W4len, ctx);
    /* W1 = W1 - W3 - W4 */
    status |= _gr_vec_sub(W1, W1, W3, len, ctx);
    status |= _gr_poly_sub(W1, W1, len, W4, W4len, ctx);
    /* W3 = W3 - W2 */
    status |= _gr_vec_sub(W3, W3, W2, len, ctx);

    /* Recomposition: */
    /* W = W4 * x^(4m) + W2*x^(3m) + W1*x^(2m) + W*x^m + W0 */

    rlen = flen + glen - 1;
    len = FLINT_MIN(rlen, m);
    status |= _gr_vec_set(res, W0, FLINT_MIN(rlen, m), ctx);
    len = FLINT_MIN(rlen - m, m);
    status |= _gr_vec_add(GR_ENTRY(res, m, sz), W3, GR_ENTRY(W0, m, sz), len, ctx);
    len = FLINT_MIN(rlen - 2 * m, m);
    status |= _gr_vec_add(GR_ENTRY(res, 2 * m, sz), W1, GR_ENTRY(W3, m, sz), len, ctx);
    len = FLINT_MIN(rlen - 3 * m, m);
    status |= _gr_vec_add(GR_ENTRY(res, 3 * m, sz), W2, GR_ENTRY(W1, m, sz), len, ctx);
    len = FLINT_MIN(rlen - 4 * m, m);
    status |= _gr_poly_add(GR_ENTRY(res, 4 * m, sz), W4, FLINT_MIN(W4len, len), GR_ENTRY(W2, m, sz), len, ctx);
    len = FLINT_MIN(rlen - 5 * m, m);
    status |= _gr_vec_set(GR_ENTRY(res, 5 * m, sz), GR_ENTRY(W4, m, sz), len, ctx);

    GR_TMP_CLEAR_VEC(tmp, alloc, ctx);

    return status;
}

int
gr_poly_mul_toom33(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    slong len_out;
    int status;

    if (poly1->length == 0 || poly2->length == 0)
        return gr_poly_zero(res, ctx);

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, len_out, ctx);
        status = _gr_poly_mul_toom33(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len_out, ctx);
        status = _gr_poly_mul_toom33(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, ctx);
    }

    _gr_poly_set_length(res, len_out, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
