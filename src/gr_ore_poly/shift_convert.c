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
is_shift_case(ore_algebra_t a)
{
    return a == ORE_ALGEBRA_FORWARD_SHIFT || a == ORE_ALGEBRA_BACKWARD_SHIFT
        || a == ORE_ALGEBRA_FORWARD_DIFFERENCE || a == ORE_ALGEBRA_BACKWARD_DIFFERENCE;
}

static int
is_backward(ore_algebra_t a)
{
    return a == ORE_ALGEBRA_BACKWARD_SHIFT || a == ORE_ALGEBRA_BACKWARD_DIFFERENCE;
}

static int
is_difference(ore_algebra_t a)
{
    return a == ORE_ALGEBRA_FORWARD_DIFFERENCE || a == ORE_ALGEBRA_BACKWARD_DIFFERENCE;
}

/* vec[m] <- sum_{k >= m} w(k, m) * binomial(k, m) * vec[k], where the sign w is
   1 for sign == 0, (-1)^(k-m) for sign == +1, and (-1)^m for sign == -1 */
static int
binomial_transform(gr_ptr vec, slong len, int sign, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem, m, k;
    fmpz_t t;

    fmpz_init(t);

    for (m = 0; m < len; m++)
    {
        fmpz_set_si(t, (sign == -1 && (m & 1)) ? -1 : 1);
        status |= gr_mul_fmpz(GR_ENTRY(vec, m, sz), GR_ENTRY(vec, m, sz), t, ctx);
        for (k = m + 1; k < len; k++)
        {
            fmpz_mul_si(t, t, (sign == 1) ? -k : k);
            fmpz_divexact_ui(t, t, (ulong) (k - m));
            status |= gr_addmul_fmpz(GR_ENTRY(vec, m, sz), GR_ENTRY(vec, k, sz), t, ctx);
        }
    }

    fmpz_clear(t);

    return status;
}

static void
reverse(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, k;

    for (k = 0; k < len / 2; k++)
        gr_swap(GR_ENTRY(vec, k, sz), GR_ENTRY(vec, len - 1 - k, sz), ctx);
}

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

    if (!is_shift_case(src_alg) || !is_shift_case(dst_alg))
        return GR_UNABLE;

    if (len <= 0)
        return GR_SUCCESS;

    status |= _gr_vec_set(res, op, len, ctx);

    if (len == 1 || src_alg == dst_alg)
        return status;

    src_back = is_backward(src_alg);
    dst_back = is_backward(dst_alg);
    crossing = (src_back != dst_back);
    s = dst_back ? -(len - 1) : (len - 1);

    if (crossing && len >= 2)  /* cases requiring a Taylor shift */
    {
        if (ctx->which_ring != GR_CTX_GR_POLY)
            return GR_UNABLE;
        else
            FLINT_ASSERT(var == 0);
    }

    if (is_difference(src_alg))
        status |= binomial_transform(res, len, src_back ? -1 : 1, ctx);

    if (crossing)
    {
        status |= taylor_shift_all(res, len, s, ctx);
        reverse(res, len, ctx);
        *p = -s;
    }

    if (is_difference(dst_alg))
        status |= binomial_transform(res, len, dst_back ? -1 : 0, ctx);

    return status;
}


int
_gr_ore_poly_shift_convert_difference(gr_ptr res, slong * p, gr_srcptr op,
                                      slong len, int to_backward, slong var,
                                      gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    *p = 0;

    if (len <= 0)
        return GR_SUCCESS;

    status |= _gr_vec_set(res, op, len, ctx);

    if (len == 1)
        return status;

    if (len >= 2)
    {
        if (ctx->which_ring != GR_CTX_GR_POLY)
            return GR_UNABLE;
        else
            FLINT_ASSERT(var == 0);
    }

    *p = to_backward ? (len - 1) : -(len - 1);

    reverse(res, len, ctx);
    status |= binomial_transform(res, len, to_backward ? 1 : 0, ctx);
    reverse(res, len, ctx);
    status |= taylor_shift_all(res, len, -*p, ctx);

    return status;
}
