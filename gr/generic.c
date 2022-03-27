#include "gr.h"
#include "gr_mat.h"
#include "gr_poly.h"

/* Generic arithmetic functions */

truth_t gr_generic_is_zero(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t eq;
    GR_TMP_START;

    GR_TMP_INIT1(t, ctx);
    eq = gr_equal(x, t, ctx);
    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;

    return eq;
}

truth_t gr_generic_is_one(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t eq;
    GR_TMP_START;

    GR_TMP_INIT1(t, ctx);

    if (gr_one(t, ctx) == GR_SUCCESS)
        eq = gr_equal(x, t, ctx);
    else
        eq = T_UNKNOWN;

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;

    return eq;
}

truth_t gr_generic_is_neg_one(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t eq;
    GR_TMP_START;

    GR_TMP_INIT1(t, ctx);

    if (gr_neg_one(t, ctx) == GR_SUCCESS)
        eq = gr_equal(x, t, ctx);
    else
        eq = T_UNKNOWN;

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;

    return eq;
}

int gr_generic_neg_one(gr_ptr res, gr_ctx_t ctx)
{
    int status;
    status = gr_one(res, ctx);
    status |= gr_neg(res, res, ctx);
    return status;
}

int gr_generic_set_fmpq(gr_ptr res, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t, u;
    int status;
    GR_TMP_START;

    status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);

    status |= gr_set_fmpz(t, fmpq_numref(y), ctx);
    status |= gr_set_fmpz(u, fmpq_denref(y), ctx);

    if (status == GR_SUCCESS)
        status = gr_inv(u, u, ctx);

    if (status == GR_SUCCESS)
        status = gr_mul(res, t, u, ctx);

    GR_TMP_CLEAR2(t, u, ctx);
    GR_TMP_END;
    return status;
}

int gr_generic_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_START;

    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    status |= gr_set_fmpz(t, y, ctx);

    if (status == GR_SUCCESS)
        status = gr_add(res, x, t, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int gr_generic_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_ui(t, y);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_si(t, y);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_START;

    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    status |= gr_set_fmpq(t, y, ctx);
    if (status == GR_SUCCESS)
        status = gr_add(res, x, t, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int gr_generic_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init(t);
    fmpz_neg_ui(t, y);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, y);
    fmpz_neg(t, t);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init(t);
    fmpz_neg(t, y);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    int status;
    fmpq_t t;
    fmpq_init(t);
    fmpq_neg(t, y);
    status = gr_add_fmpq(res, x, t, ctx);
    fmpq_clear(t);
    return status;
}

int gr_generic_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_START;

    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    status |= gr_set_fmpz(t, y, ctx);

    if (status == GR_SUCCESS)
        status = gr_mul(res, x, t, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int gr_generic_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_ui(t, y);
    status = gr_mul_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_si(t, y);
    status = gr_mul_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_START;

    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    status |= gr_set_fmpq(t, y, ctx);
    if (status == GR_SUCCESS)
        status = gr_mul(res, x, t, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int gr_generic_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    return gr_add(res, x, x, ctx);
}

int gr_generic_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    return gr_mul(res, x, x, ctx);
}

int gr_generic_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_one(x, ctx) == T_TRUE)
        return gr_one(res, ctx);

    if (gr_is_neg_one(x, ctx) == T_TRUE)
        return gr_neg_one(res, ctx);

    /* todo: dubious in the zero ring, if comparing with 1 above
       somehow failed */
    if (gr_is_zero(x, ctx) == T_TRUE)
        return GR_DOMAIN;

    return GR_UNABLE;
}

truth_t gr_generic_is_invertible(gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_one(x, ctx) == T_TRUE)
        return T_TRUE;

    if (gr_is_neg_one(x, ctx) == T_TRUE)
        return T_TRUE;

    /* todo: dubious in the zero ring, if comparing with 1 above
       somehow failed */
    if (gr_is_zero(x, ctx) == T_TRUE)
        return T_FALSE;

    return T_UNKNOWN;
}


int gr_generic_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_START;

    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    status |= gr_set_fmpz(t, y, ctx);

    if (status == GR_SUCCESS)
        status = gr_div(res, x, t, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int gr_generic_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_ui(t, y);
    status = gr_div_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_si(t, y);
    status = gr_div_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_START;

    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    status |= gr_set_fmpq(t, y, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(res, x, t, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}


/* Generic powering */

/* Assumes exp >= 2; res and tmp not not aliased with x. */
static int
gr_generic_pow_ui_binexp(gr_ptr res, gr_ptr tmp, gr_srcptr x, ulong exp, gr_ctx_t ctx)
{
    gr_ptr R, S, T;
    gr_method_unary_op sqr = GR_UNARY_OP(ctx, SQR);
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    int status;
    int zeros;
    ulong bit;

    status = GR_SUCCESS;

    /* Determine parity due to swaps */
    zeros = 0;
    bit = exp;
    while (bit > 1)
    {
        zeros += !(bit & 1);
        bit >>= 1;
    }

    if (zeros % 2)
    {
        R = res;
        S = tmp;
    }
    else
    {
        R = tmp;
        S = res;
    }

    bit = UWORD(1) << (FLINT_BIT_COUNT(exp) - 2);

    status |= sqr(R, x, ctx);

    if (bit & exp)
    {
        status |= mul(S, R, x, ctx);
        T = R;
        R = S;
        S = T;
    }

    while (bit >>= 1)
    {
        status |= sqr(S, R, ctx);

        if (bit & exp)
        {
            status |= mul(R, S, x, ctx);
        }
        else
        {
            T = R;
            R = S;
            S = T;
        }
    }

    return status;
}

/* todo: optimize swaps */
static int
gr_generic_pow_fmpz_binexp(gr_ptr res, gr_srcptr x, const fmpz_t exp, gr_ctx_t ctx)
{
    gr_ptr t, u;
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_unary_op sqr = GR_UNARY_OP(ctx, SQR);
    gr_method_swap_op swap = GR_SWAP_OP(ctx, SWAP);
    int status;
    slong i;
    GR_TMP_START;

    status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);

    status |= gr_set(t, x, ctx);

    for (i = fmpz_bits(exp) - 2; i >= 0; i--)
    {
        status |= sqr(u, t, ctx);

        if (fmpz_tstbit(exp, i))
            status |= mul(t, u, x, ctx);
        else
            swap(t, u, ctx);
    }

    swap(res, t, ctx);

    GR_TMP_CLEAR2(t, u, ctx);
    GR_TMP_END;

    return status;
}

int
gr_generic_pow_ui(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
{
    int status;

    if (e > (ulong) ctx->size_limit)
        return GR_UNABLE;

    if (e == 0)
    {
        return gr_one(res, ctx);
    }
    else if (e == 1)
    {
        return gr_set(res, x, ctx);
    }
    else if (e == 2)
    {
        return gr_sqr(res, x, ctx);
    }
    else if (e == 4)
    {
        status = gr_sqr(res, x, ctx);
        status |= gr_sqr(res, res, ctx);
        return status;
    }
    else
    {
        gr_ptr t, u;
        GR_TMP_START;

        if (res == x)
        {
            GR_TMP_INIT2(t, u, ctx);
            gr_set(u, x, ctx);
            status = gr_generic_pow_ui_binexp(res, t, u, e, ctx);
            GR_TMP_CLEAR2(t, u, ctx);
        }
        else
        {
            GR_TMP_INIT1(t, ctx);
            status = gr_generic_pow_ui_binexp(res, t, x, e, ctx);
            GR_TMP_CLEAR1(t, ctx);
        }

        GR_TMP_END;
        return status;
    }
}

/* todo: call gr_pow_ui instead of gr_generic_pow_ui? */
int
gr_generic_pow_si(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    if (e >= 0)
    {
        return gr_generic_pow_ui(res, x, e, ctx);
    }
    else
    {
        int status;

        /* todo: some heuristic for when we want to invert before/after powering */
        status = gr_inv(res, x, ctx);
        if (status == GR_SUCCESS)
            status = gr_generic_pow_ui(res, res, -e, ctx);

        return status;
    }
}

int
gr_generic_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t e, gr_ctx_t ctx)
{
    int status;

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);

        /* todo: some heuristic for when we want to invert before/after powering */
        status = gr_inv(res, x, ctx);
        if (status == GR_SUCCESS)
            status = gr_generic_pow_fmpz(res, res, f, ctx);

        fmpz_clear(f);
        return status;
    }

    if (*e == 0)
        return gr_one(res, ctx);
    else if (*e == 1)
        return gr_set(res, x, ctx);
    else if (*e == 2)
        return gr_sqr(res, x, ctx);
    else
        return gr_generic_pow_fmpz_binexp(res, x, e, ctx);
}

truth_t
gr_generic_is_square(gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return T_TRUE;

    if (gr_is_one(x, ctx) == T_TRUE)
        return T_TRUE;

    return GR_UNABLE;
}

int
gr_generic_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
    {
        gr_zero(res, ctx);
        return GR_SUCCESS;
    }

    if (gr_is_one(x, ctx) == T_TRUE)
    {
        gr_one(res, ctx);
        return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int
gr_generic_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return GR_DOMAIN;

    if (gr_is_one(x, ctx) == T_TRUE)
    {
        gr_one(res, ctx);
        return GR_SUCCESS;
    }

    return GR_UNABLE;
}


/* Generic vector functions */

void
gr_generic_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT);
    slong i, sz;

    sz = ctx->sizeof_elem;

    for (i = 0; i < len; i++)
        init(GR_ENTRY(vec, i, sz), ctx);
}

void
gr_generic_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR);
    slong i, sz;

    sz = ctx->sizeof_elem;

    for (i = 0; i < len; i++)
        clear(GR_ENTRY(vec, i, sz), ctx);
}

void
gr_generic_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx)
{
    gr_method_swap_op swap = GR_SWAP_OP(ctx, SWAP);
    slong i, sz;

    sz = ctx->sizeof_elem;

    for (i = 0; i < len; i++)
        swap(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);
}

int
gr_generic_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_constant_op zero = GR_CONSTANT_OP(ctx, ZERO);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= zero(GR_ENTRY(vec, i, sz), ctx);

    return status;
}

int
gr_generic_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)
{
    gr_method_unary_op set = GR_UNARY_OP(ctx, SET);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= set(GR_ENTRY(res, i, sz), GR_ENTRY(src, i, sz), ctx);

    return status;
}

int
gr_generic_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)
{
    gr_method_unary_op neg = GR_UNARY_OP(ctx, NEG);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= neg(GR_ENTRY(res, i, sz), GR_ENTRY(src, i, sz), ctx);

    return status;
}

int
gr_generic_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= add(GR_ENTRY(res, i, sz), GR_ENTRY(src1, i, sz), GR_ENTRY(src2, i, sz), ctx);

    return status;
}

int
gr_generic_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= sub(GR_ENTRY(res, i, sz), GR_ENTRY(src1, i, sz), GR_ENTRY(src2, i, sz), ctx);

    return status;
}

int
gr_generic_vec_scalar_mul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= mul(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), c, ctx);

    return status;
}

int
gr_generic_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;
    GR_TMP_START;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= mul(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= add(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);
    int status;
    slong i, sz;
    gr_ptr t;
    GR_TMP_START;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= mul(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= sub(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
{
    gr_method_binary_op_si mul_si = GR_BINARY_OP_SI(ctx, MUL_SI);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;
    GR_TMP_START;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= mul_si(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= add(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
{
    gr_method_binary_op_si mul_si = GR_BINARY_OP_SI(ctx, MUL_SI);
    gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);
    int status;
    slong i, sz;
    gr_ptr t;
    GR_TMP_START;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= mul_si(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= sub(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

truth_t
gr_generic_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_predicate equal = GR_BINARY_PREDICATE(ctx, EQUAL);
    truth_t eq, this_eq;
    slong i, sz;

    sz = ctx->sizeof_elem;

    eq = T_TRUE;

    for (i = 0; i < len; i++)
    {
        this_eq = equal(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);

        if (this_eq == T_FALSE)
            return T_FALSE;

        if (this_eq == T_UNKNOWN)
            eq = T_UNKNOWN;
    }

    return eq;
}

int
gr_generic_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);
    truth_t eq, this_eq;
    slong i, sz;

    sz = ctx->sizeof_elem;

    eq = T_TRUE;

    for (i = 0; i < len; i++)
    {
        this_eq = is_zero(GR_ENTRY(vec, i, sz), ctx);

        if (this_eq == T_FALSE)
            return T_FALSE;

        if (this_eq == T_UNKNOWN)
            eq = T_UNKNOWN;
    }

    return eq;
}

int
gr_generic_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;
    GR_TMP_START;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    if (initial == NULL)
    {
        status |= mul(res, vec1, vec2, ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= mul(t, vec1, vec2, ctx);
        status |= add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= mul(t, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);
        status |= add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;
    GR_TMP_START;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    if (initial == NULL)
    {
        status |= mul(res, vec1, GR_ENTRY(vec2, len - 1, sz), ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= mul(t, vec1, GR_ENTRY(vec2, len - 1, sz), ctx);
        status |= add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= mul(t, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, len - 1 - i, sz), ctx);
        status |= add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_ctx_clear(gr_ctx_t ctx)
{
    return GR_SUCCESS;
}

/* Generic method implementations */
const gr_method_tab_input _gr_generic_methods[] =
{
    {GR_METHOD_CTX_CLEAR,               (gr_funcptr) gr_generic_ctx_clear},

    {GR_METHOD_NEG_ONE,                 (gr_funcptr) gr_generic_neg_one},

    {GR_METHOD_IS_ZERO,                 (gr_funcptr) gr_generic_is_zero},
    {GR_METHOD_IS_ONE,                  (gr_funcptr) gr_generic_is_one},
    {GR_METHOD_IS_NEG_ONE,              (gr_funcptr) gr_generic_is_neg_one},

    {GR_METHOD_SET_FMPQ,                (gr_funcptr) gr_generic_set_fmpq},

    {GR_METHOD_ADD_UI,                  (gr_funcptr) gr_generic_add_ui},
    {GR_METHOD_ADD_SI,                  (gr_funcptr) gr_generic_add_si},
    {GR_METHOD_ADD_FMPZ,                (gr_funcptr) gr_generic_add_fmpz},
    {GR_METHOD_ADD_FMPQ,                (gr_funcptr) gr_generic_add_fmpq},

    {GR_METHOD_SUB_UI,                  (gr_funcptr) gr_generic_sub_ui},
    {GR_METHOD_SUB_SI,                  (gr_funcptr) gr_generic_sub_si},
    {GR_METHOD_SUB_FMPZ,                (gr_funcptr) gr_generic_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,                (gr_funcptr) gr_generic_sub_fmpq},

    {GR_METHOD_MUL_UI,                  (gr_funcptr) gr_generic_mul_ui},
    {GR_METHOD_MUL_SI,                  (gr_funcptr) gr_generic_mul_si},
    {GR_METHOD_MUL_FMPZ,                (gr_funcptr) gr_generic_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,                (gr_funcptr) gr_generic_mul_fmpq},

    {GR_METHOD_DIV_UI,                  (gr_funcptr) gr_generic_div_ui},
    {GR_METHOD_DIV_SI,                  (gr_funcptr) gr_generic_div_si},
    {GR_METHOD_DIV_FMPZ,                (gr_funcptr) gr_generic_div_fmpz},
    {GR_METHOD_DIV_FMPQ,                (gr_funcptr) gr_generic_div_fmpq},

    {GR_METHOD_IS_INVERTIBLE,           (gr_funcptr) gr_generic_is_invertible},
    {GR_METHOD_INV,                     (gr_funcptr) gr_generic_inv},

    {GR_METHOD_POW_SI,                  (gr_funcptr) gr_generic_pow_si},
    {GR_METHOD_POW_UI,                  (gr_funcptr) gr_generic_pow_ui},
    {GR_METHOD_POW_FMPZ,                (gr_funcptr) gr_generic_pow_fmpz},

    {GR_METHOD_IS_SQUARE,               (gr_funcptr) gr_generic_is_square},
    {GR_METHOD_SQRT,                    (gr_funcptr) gr_generic_sqrt},
    {GR_METHOD_RSQRT,                   (gr_funcptr) gr_generic_rsqrt},

    {GR_METHOD_VEC_INIT,                (gr_funcptr) gr_generic_vec_init},
    {GR_METHOD_VEC_CLEAR,               (gr_funcptr) gr_generic_vec_clear},
    {GR_METHOD_VEC_SWAP,                (gr_funcptr) gr_generic_vec_swap},
    {GR_METHOD_VEC_ZERO,                (gr_funcptr) gr_generic_vec_zero},
    {GR_METHOD_VEC_SET,                 (gr_funcptr) gr_generic_vec_set},
    {GR_METHOD_VEC_NEG,                 (gr_funcptr) gr_generic_vec_neg},
    {GR_METHOD_VEC_ADD,                 (gr_funcptr) gr_generic_vec_add},
    {GR_METHOD_VEC_SUB,                 (gr_funcptr) gr_generic_vec_sub},
    {GR_METHOD_VEC_SCALAR_MUL,          (gr_funcptr) gr_generic_vec_scalar_mul},
    {GR_METHOD_VEC_SCALAR_ADDMUL,       (gr_funcptr) gr_generic_vec_scalar_addmul},
    {GR_METHOD_VEC_SCALAR_SUBMUL,       (gr_funcptr) gr_generic_vec_scalar_submul},
    {GR_METHOD_VEC_SCALAR_ADDMUL_SI,    (gr_funcptr) gr_generic_vec_scalar_addmul_si},
    {GR_METHOD_VEC_SCALAR_SUBMUL_SI,    (gr_funcptr) gr_generic_vec_scalar_submul_si},
    {GR_METHOD_VEC_EQUAL,               (gr_funcptr) gr_generic_vec_equal},
    {GR_METHOD_VEC_IS_ZERO,             (gr_funcptr) gr_generic_vec_is_zero},
    {GR_METHOD_VEC_DOT,                 (gr_funcptr) gr_generic_vec_dot},
    {GR_METHOD_VEC_DOT_REV,             (gr_funcptr) gr_generic_vec_dot_rev},

    {GR_METHOD_POLY_MULLOW,             (gr_funcptr) _gr_poly_mullow_generic},

    {0,                                 (gr_funcptr) NULL},
};

void
gr_method_tab_init(gr_funcptr * methods, gr_method_tab_input * tab)
{
    slong i;

    for (i = 0; i < GR_METHOD_TAB_SIZE; i++)
        methods[i] = (gr_funcptr) gr_not_implemented;

    /* Assign generic methods as fallbacks */
    for (i = 0; ; i++)
    {
        if (_gr_generic_methods[i].function == NULL)
            break;

        if (_gr_generic_methods[i].index >= GR_METHOD_TAB_SIZE)
            abort();

        methods[_gr_generic_methods[i].index] = _gr_generic_methods[i].function;
    }

    for (i = 0; ; i++)
    {
        if (tab[i].function == NULL)
            break;

        if (tab[i].index >= GR_METHOD_TAB_SIZE)
            abort();

        methods[tab[i].index] = tab[i].function;
    }
}
