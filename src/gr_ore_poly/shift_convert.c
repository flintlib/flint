/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "fmpz.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_vec.h"
#include "gr_ore_poly.h"

static int
shift_alg_is_shift(ore_algebra_t a)
{
    return a == ORE_ALGEBRA_FORWARD_SHIFT || a == ORE_ALGEBRA_BACKWARD_SHIFT
        || a == ORE_ALGEBRA_FORWARD_DIFFERENCE || a == ORE_ALGEBRA_BACKWARD_DIFFERENCE;
}

static int
shift_alg_is_backward(ore_algebra_t a)
{
    return a == ORE_ALGEBRA_BACKWARD_SHIFT || a == ORE_ALGEBRA_BACKWARD_DIFFERENCE;
}

static int
shift_alg_is_difference(ore_algebra_t a)
{
    return a == ORE_ALGEBRA_FORWARD_DIFFERENCE || a == ORE_ALGEBRA_BACKWARD_DIFFERENCE;
}

/* In-place change of basis between a shift generator and the associated
   finite-difference generator, computing
       vec[m] <- sum_{k >= m} w(k, m) * binomial(k, m) * vec[k],
   where the sign w is 1 for sign == 0, (-1)^(k-m) for sign == +1, and (-1)^m
   for sign == -1. The binomial coefficients along each row are built up
   incrementally, as in _arb_poly_binomial_transform_basecase. The update reads
   only entries of index >= m, so writing vec[m] in place is safe. */
static int
binomial_transform(gr_ptr vec, slong len, int sign, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem, m, k;
    gr_ptr acc;
    fmpz_t t;

    GR_TMP_INIT(acc, ctx);
    fmpz_init(t);

    for (m = 0; m < len; m++)
    {
        status |= gr_zero(acc, ctx);
        for (k = m; k < len; k++)
        {
            if (k == m)
                fmpz_set_si(t, (sign == -1 && (m & 1)) ? -1 : 1);
            else
            {
                fmpz_mul_si(t, t, (sign == 1) ? -k : k);
                fmpz_divexact_ui(t, t, (ulong) (k - m));
            }
            status |= gr_addmul_fmpz(acc, GR_ENTRY(vec, k, sz), t, ctx);
        }
        status |= gr_set(GR_ENTRY(vec, m, sz), acc, ctx);
    }

    fmpz_clear(t);
    GR_TMP_CLEAR(acc, ctx);

    return status;
}

static void
reverse(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, k;

    for (k = 0; k < len / 2; k++)
        gr_swap(GR_ENTRY(vec, k, sz), GR_ENTRY(vec, len - 1 - k, sz), ctx);
}

/* Taylor-shift every coefficient by p; requires a polynomial base ring. */
static int
taylor_shift_all(gr_ptr vec, slong len, slong p, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem, k;
    gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(ctx);
    gr_ptr c;

    GR_TMP_INIT(c, sctx);
    status |= gr_set_si(c, p, sctx);
    for (k = 0; k < len; k++)
    {
        gr_poly_struct * vk = (gr_poly_struct *) GR_ENTRY(vec, k, sz);
        status |= gr_poly_taylor_shift(vk, vk, c, sctx);
    }
    GR_TMP_CLEAR(c, sctx);

    return status;
}

int
_gr_ore_poly_shift_convert(gr_ptr res, slong * p, gr_srcptr op, slong len,
                           ore_algebra_t src_alg, ore_algebra_t dst_alg, slong var, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    int src_back, dst_back, crossing;
    slong s;

    *p = 0;

    if (!shift_alg_is_shift(src_alg) || !shift_alg_is_shift(dst_alg))
        return GR_DOMAIN;

    if (len <= 0)
        return GR_SUCCESS;

    /* The implementation is limited to a univariate polynomial base ring, whose
       only generator has index 0. */
    if (var != 0)
        return GR_UNABLE;

    if (src_alg == dst_alg)
        return _gr_vec_set(res, op, len, ctx);

    src_back = shift_alg_is_backward(src_alg);
    dst_back = shift_alg_is_backward(dst_alg);
    crossing = (src_back != dst_back);
    /* res is computed as S^s * op, so the reported p (satisfying S^p * res = op)
       is -s. */
    s = dst_back ? -(len - 1) : (len - 1);

    /* Crossing between the forward side S, S - 1 and the backward side S^{-1},
       1 - S^{-1} conjugates by a power of S, realized as a Taylor shift of the
       coefficients; this needs a polynomial base ring once len >= 2. */
    if (crossing && len >= 2 && ctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    status |= _gr_vec_set(res, op, len, ctx);

    /* The two difference generators sit on opposite sides, so converting between
       them always crosses. Because S - 1 = S (1 - S^{-1}), a single binomial
       transform suffices: reversing the coefficients turns the combination into
       the standard one handled by binomial_transform. */
    if (shift_alg_is_difference(src_alg) && shift_alg_is_difference(dst_alg))
    {
        reverse(res, len, ctx);
        status |= binomial_transform(res, len, dst_back ? 1 : 0, ctx);
        reverse(res, len, ctx);
        if (s != 0)
            status |= taylor_shift_all(res, len, s, ctx);
        *p = -s;
        return status;
    }

    /* General pipeline: source difference -> its shift, cross sides if needed,
       then shift -> destination difference. At most one binomial transform runs
       (the difference <-> difference case above is the only one needing two). */
    if (shift_alg_is_difference(src_alg))
        status |= binomial_transform(res, len, src_back ? -1 : 1, ctx);

    if (crossing)
    {
        if (s != 0)
            status |= taylor_shift_all(res, len, s, ctx);
        reverse(res, len, ctx);
        *p = -s;
    }

    if (shift_alg_is_difference(dst_alg))
        status |= binomial_transform(res, len, dst_back ? -1 : 0, ctx);

    return status;
}
