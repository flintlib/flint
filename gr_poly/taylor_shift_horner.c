#include "gr_poly.h"

/* todo: use vector functions when overlap behavior is guaranteed? */
int
_gr_poly_taylor_shift_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    int status;
    slong i, j;
    slong sz = ctx->sizeof_elem;

    status = GR_SUCCESS;

    if (res != poly)
        status |= _gr_vec_set(res, poly, len, ctx);

    if (gr_is_one(c, ctx) == T_TRUE)
    {
        gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);

        for (i = len - 2; i >= 0; i--)
            for (j = i; j < len - 1; j++)
                status |= add(GR_ENTRY(res, j, sz), GR_ENTRY(res, j, sz), GR_ENTRY(res, j + 1, sz), ctx);
    }
    else if (gr_is_neg_one(c, ctx) == T_TRUE)
    {
        gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);

        for (i = len - 2; i >= 0; i--)
            for (j = i; j < len - 1; j++)
                status |= sub(GR_ENTRY(res, j, sz), GR_ENTRY(res, j, sz), GR_ENTRY(res, j + 1, sz), ctx);
    }
    else if (gr_is_zero(c, ctx) != T_TRUE)
    {
        /* todo: when we have generic addmul, use a temporary */
        gr_method_binary_op addmul = GR_BINARY_OP(ctx, ADDMUL);

        for (i = len - 2; i >= 0; i--)
            for (j = i; j < len - 1; j++)
                status |= addmul(GR_ENTRY(res, j, sz), GR_ENTRY(res, j + 1, sz), c, ctx);
    }

    return status;
}

int
gr_poly_taylor_shift_horner(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (res != f)
        status |= gr_poly_set(res, f, ctx);

    status |= _gr_poly_taylor_shift_horner(res->coeffs, res->coeffs, res->length, c, ctx);
    return status;
}
