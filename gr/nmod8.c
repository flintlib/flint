#include "gr.h"

#define NMOD8_CTX(ring_ctx) (((nmod_t *)((ring_ctx)->elem_ctx))[0])

typedef unsigned char nmod8_struct;
typedef nmod8_struct nmod8_t[1];

void
nmod8_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_ui(out, NMOD8_CTX(ctx).n);
    gr_stream_write(out, " (nmod8)");
}

void
nmod8_init(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = 0;
}

void
nmod8_clear(nmod8_t x, const gr_ctx_t ctx)
{
}

void
nmod8_swap(nmod8_t x, nmod8_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

int
nmod8_randtest(nmod8_t res, flint_rand_t state, const void * options, const gr_ctx_t ctx)
{
    res[0] = n_randtest(state) % NMOD8_CTX(ctx).n;
    return GR_SUCCESS;
}

int
nmod8_write(gr_stream_t out, const nmod8_t x, const gr_ctx_t ctx)
{
    gr_stream_write_si(out, x[0]);
    return GR_SUCCESS;
}

int
nmod8_zero(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = 0;
    return GR_SUCCESS;
}

int
nmod8_one(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = (NMOD8_CTX(ctx).n != 0);
    return GR_SUCCESS;
}

int
nmod8_set_si(nmod8_t res, slong v, const gr_ctx_t ctx)
{
    ulong t;
    nmod_t mod = NMOD8_CTX(ctx);
    t = FLINT_ABS(v);
    NMOD_RED(t, t, mod);
    if (v < 0)
        t = nmod_neg(t, mod);
    res[0] = t;
    return GR_SUCCESS;
}

int
nmod8_set_ui(nmod8_t res, ulong v, const gr_ctx_t ctx)
{
    ulong t;
    nmod_t mod = NMOD8_CTX(ctx);
    NMOD_RED(t, v, mod);
    res[0] = t;
    return GR_SUCCESS;
}

int
nmod8_set_fmpz(nmod8_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    nmod_t mod = NMOD8_CTX(ctx);
    res[0] = fmpz_get_nmod(v, mod);
    return GR_SUCCESS;
}

truth_t
nmod8_is_zero(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == 0) ? T_TRUE : T_FALSE;
}

truth_t
nmod8_is_one(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == (NMOD8_CTX(ctx).n != 1)) ? T_TRUE : T_FALSE;
}

truth_t
nmod8_is_neg_one(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == NMOD8_CTX(ctx).n - 1) ? T_TRUE : T_FALSE;
}

truth_t
nmod8_equal(const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    return (x[0] == y[0]) ? T_TRUE : T_FALSE;
}

int
nmod8_set(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    res[0] = x[0];
    return GR_SUCCESS;
}

int
nmod8_neg(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    res[0] = nmod_neg(x[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

int
nmod8_add(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_add(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

int
nmod8_add_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_add(res, x, t, ctx);
}

int
nmod8_sub(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_sub(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

int
nmod8_mul(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_mul(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

int
nmod8_mul_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_mul(res, x, t, ctx);
}

int
nmod8_inv(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    ulong r, g;

    g = n_gcdinv(&r, x[0], NMOD8_CTX(ctx).n);
    if (g == 1)
    {
        res[0] = r;
        return GR_SUCCESS;
    }
    else
    {
        res[0] = 0;
        return GR_DOMAIN;
    }
}

int
nmod8_div(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    int status;

    status = nmod8_inv(t, y, ctx);
    nmod8_mul(res, x, t, ctx);
    return status;
}

int
nmod8_div_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

truth_t
nmod8_is_invertible(const nmod8_t x, const gr_ctx_t ctx)
{
    ulong r, g;
    g = n_gcdinv(&r, x[0], NMOD8_CTX(ctx).n);
    return (g == 1) ? T_TRUE : T_FALSE;
}

int
_nmod8_vec_dot(nmod8_t res, const nmod8_t initial, int subtract, const nmod8_struct * vec1, const nmod8_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n, s;

    if (len <= 0)
    {
        if (initial == NULL)
            nmod8_zero(res, ctx);
        else
            nmod8_set(res, initial, ctx);
        return GR_SUCCESS;
    }

    n = NMOD8_CTX(ctx).n;

    if (initial == NULL)
    {
        s = 0;
    }
    else
    {
        if (subtract)
            s = (n - initial[0]);
        else
            s = initial[0];
    }

    for (i = 0; i + 4 < len; i += 4)
    {
        s += (ulong) vec1[i    ] * (ulong) vec2[i];
        s += (ulong) vec1[i + 1] * (ulong) vec2[i + 1];
        s += (ulong) vec1[i + 2] * (ulong) vec2[i + 2];
        s += (ulong) vec1[i + 3] * (ulong) vec2[i + 3];
    }

    for ( ; i < len; i++)
        s += (ulong) vec1[i] * (ulong) vec2[i];

    nmod8_set_ui(res, s, ctx);

    if (subtract && res[0] != 0)
        res[0] = (n - res[0]);

    return GR_SUCCESS;
}

int
_nmod8_vec_dot_rev(nmod8_t res, const nmod8_t initial, int subtract, const nmod8_struct * vec1, const nmod8_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n, s;

    if (len <= 0)
    {
        if (initial == NULL)
            nmod8_zero(res, ctx);
        else
            nmod8_set(res, initial, ctx);
        return GR_SUCCESS;
    }

    n = NMOD8_CTX(ctx).n;

    if (initial == NULL)
    {
        s = 0;
    }
    else
    {
        if (subtract)
            s = (n - initial[0]);
        else
            s = initial[0];
    }

    for (i = 0; i + 4 < len; i += 4)
    {
        s += (ulong) vec1[i    ] * (ulong) vec2[len - 1 - i];
        s += (ulong) vec1[i + 1] * (ulong) vec2[len - 1 - i - 1];
        s += (ulong) vec1[i + 2] * (ulong) vec2[len - 1 - i - 2];
        s += (ulong) vec1[i + 3] * (ulong) vec2[len - 1 - i - 3];
    }

    for ( ; i < len; i++)
        s += (ulong) vec1[i] * (ulong) vec2[len - 1 - i];

    nmod8_set_ui(res, s, ctx);

    if (subtract && res[0] != 0)
        res[0] = (n - res[0]);

    return GR_SUCCESS;
}

int
nmod8_ctx_clear(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
    return GR_SUCCESS;
}

int _nmod8_methods2_initialized = 0;
gr_static_method_table _nmod8_static_table;
gr_method_tab_t _nmod8_methods2;

gr_method_tab_input nmod8_methods2[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) nmod8_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) nmod8_ctx_clear},
    {GR_METHOD_INIT,            (gr_funcptr) nmod8_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) nmod8_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) nmod8_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) nmod8_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) nmod8_write},
    {GR_METHOD_ZERO,            (gr_funcptr) nmod8_zero},
    {GR_METHOD_ONE,             (gr_funcptr) nmod8_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) nmod8_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) nmod8_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) nmod8_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) nmod8_equal},
    {GR_METHOD_SET,             (gr_funcptr) nmod8_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) nmod8_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) nmod8_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) nmod8_set_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) nmod8_neg},
    {GR_METHOD_ADD,             (gr_funcptr) nmod8_add},
    {GR_METHOD_ADD_SI,          (gr_funcptr) nmod8_add_si},
    {GR_METHOD_SUB,             (gr_funcptr) nmod8_sub},
    {GR_METHOD_MUL,             (gr_funcptr) nmod8_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) nmod8_mul_si},
    {GR_METHOD_DIV,             (gr_funcptr) nmod8_div},
    {GR_METHOD_DIV_SI,          (gr_funcptr) nmod8_div_si},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) nmod8_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) nmod8_inv},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _nmod8_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _nmod8_vec_dot_rev},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n)
{
    ctx->flags = GR_COMMUTATIVE_RING;
    ctx->which_ring = GR_WHICH_RING_ZZ_MOD;
    ctx->sizeof_elem = sizeof(nmod8_struct);
    ctx->size_limit = WORD_MAX;

    ctx->elem_ctx = flint_malloc(sizeof(nmod_t));  /* This could be something more interesting */
    nmod_init(ctx->elem_ctx, n);

    if (!_nmod8_methods2_initialized)
    {
        gr_method_tab_init_static(&_nmod8_methods2, _nmod8_static_table, nmod8_methods2);
        _nmod8_methods2_initialized = 1;
    }

    ctx->methods2 = &_nmod8_methods2;
}
