#include "gr.h"
#include "flint/fmpz_mod.h"
#include "flint/fmpz_mod_poly.h"
#include "flint/fmpz_mod_mat.h"

typedef struct
{
    fmpz_mod_ctx_struct ctx;
    truth_t is_prime;
}
fmpz_mod_ctx_extended_struct;

#define FMPZ_MOD_CTX(ring_ctx) (&(((fmpz_mod_ctx_extended_struct *)((ring_ctx)->elem_ctx))->ctx))
#define FMPZ_MOD_IS_PRIME(ring_ctx) ((((fmpz_mod_ctx_extended_struct *)((ring_ctx)->elem_ctx))->is_prime))

int
_gr_fmpz_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_fmpz(out, FMPZ_MOD_CTX(ctx)->n);
    gr_stream_write(out, " (fmpz)");
    return GR_SUCCESS;
}

void
_gr_fmpz_mod_ctx_clear(gr_ctx_t ctx)
{
    fmpz_mod_ctx_clear(FMPZ_MOD_CTX(ctx));
    flint_free(ctx->elem_ctx);
}

truth_t
_gr_fmpz_mod_ctx_is_field(gr_ctx_t ctx)
{
/*
    if (FMPZ_MOD_IS_PRIME(ctx) == T_UNKNOWN)
        FMPZ_MOD_IS_PRIME(ctx) = fmpz_is_prime(FMPZ_MOD_CTX(ctx)->n) ? T_TRUE : T_FALSE;
*/

    return FMPZ_MOD_IS_PRIME(ctx);
}

void
_gr_fmpz_mod_init(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_init(x);
}

void
_gr_fmpz_mod_clear(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_clear(x);
}

void
_gr_fmpz_mod_swap(fmpz_t x, fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

int
_gr_fmpz_mod_randtest(fmpz_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fmpz_mod_rand(res, state, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_write(gr_stream_t out, const fmpz_t x, const gr_ctx_t ctx)
{
    gr_stream_write_fmpz(out, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_zero(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_zero(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_one(fmpz_t x, const gr_ctx_t ctx)
{
    if (fmpz_is_one(FMPZ_MOD_CTX(ctx)->n))
        fmpz_zero(x);
    else
        fmpz_one(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_set_si(fmpz_t res, slong v, const gr_ctx_t ctx)
{
    fmpz_mod_set_si(res, v, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_set_ui(fmpz_t res, ulong v, const gr_ctx_t ctx)
{
    fmpz_mod_set_ui(res, v, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_set_fmpz(fmpz_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fmpz_mod_set_fmpz(res, v, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mod_is_zero(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mod_is_one(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_mod_is_one(x, FMPZ_MOD_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mod_is_neg_one(const fmpz_t x, const gr_ctx_t ctx)
{
    truth_t res;
    fmpz_t t;
    fmpz_init(t);
    fmpz_mod_set_si(t, -1, FMPZ_MOD_CTX(ctx));
    res = fmpz_equal(t, x) ? T_TRUE : T_FALSE;
    fmpz_clear(t);
    return res;
}

truth_t
_gr_fmpz_mod_equal(const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    return fmpz_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mod_set(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_set(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_neg(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_mod_neg(res, x, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_add(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_mod_add(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_add_si(fmpz_t res, const fmpz_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_mod_add_si(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_add_ui(fmpz_t res, const fmpz_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_mod_add_ui(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_sub(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_mod_sub(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
#if 1
    fmpz_mod_mul(res, x, y, FMPZ_MOD_CTX(ctx));
#else
    fmpz_mul(res, x, y);
    fmpz_mod(res, res, FMPZ_MOD_CTX(ctx)->n);
#endif
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mul_si(fmpz_t res, const fmpz_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_mod_mul_si(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_addmul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul(t, x, y);
    fmpz_add(t, t, res);
    fmpz_mod(res, t, FMPZ_MOD_CTX(ctx)->n);
    fmpz_clear(t);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_submul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul(t, x, y);
    fmpz_sub(t, res, t);
    fmpz_mod(res, t, FMPZ_MOD_CTX(ctx)->n);
    fmpz_clear(t);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mul_two(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_mod_add(res, x, x, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_sqr(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_mod_mul(res, x, x, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_inv(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    /* todo: also check for -1 when fast? */
    if (fmpz_is_one(x))
    {
        fmpz_one(res);
        return GR_SUCCESS;
    }
    else if (fmpz_is_zero(x))
    {
        fmpz_zero(res);
        return fmpz_is_one(FMPZ_MOD_CTX(ctx)->n) ? GR_SUCCESS : GR_DOMAIN;
    }
    else
    {
        int status;

        fmpz_t d;
        fmpz_init(d);
        fmpz_gcdinv(d, res, x, FMPZ_MOD_CTX(ctx)->n);

        if (fmpz_is_one(d))
            status = GR_SUCCESS;
        else
            status = GR_DOMAIN;

        fmpz_clear(d);
        return status;
    }
}

int
_gr_fmpz_mod_div(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    int status;
    fmpz_t t;

    fmpz_init(t);
    status = _gr_fmpz_mod_inv(t, y, ctx);
    if (status == GR_SUCCESS)
        fmpz_mod_mul(res, x, t, FMPZ_MOD_CTX(ctx));

    fmpz_clear(t);
    return status;
}

truth_t
_gr_fmpz_mod_is_invertible(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_mod_is_invertible(x, FMPZ_MOD_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mod_pow_ui(fmpz_t res, const fmpz_t x, ulong exp, const gr_ctx_t ctx)
{
    fmpz_mod_pow_ui(res, x, exp, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_pow_fmpz(fmpz_t res, const fmpz_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    return fmpz_mod_pow_fmpz(res, x, exp, FMPZ_MOD_CTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

/* todo: len 1 */
int
_gr_fmpz_mod_vec_dot(fmpz_t res, const fmpz_t initial, int subtract, const fmpz * vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            fmpz_zero(res);
        else
            fmpz_set(res, initial);
        return GR_SUCCESS;
    }

    if (initial == NULL)
    {
        fmpz_mul(res, vec1, vec2);
    }
    else
    {
        if (subtract)
            fmpz_neg(res, initial);
        else
            fmpz_set(res, initial);

        fmpz_addmul(res, vec1, vec2);
    }

    for (i = 1; i < len; i++)
        fmpz_addmul(res, vec1 + i, vec2 + i);

    if (subtract)
        fmpz_neg(res, res);

    fmpz_mod(res, res, FMPZ_MOD_CTX(ctx)->n);

    return GR_SUCCESS;
}

/* todo: len 1 */
int
_gr_fmpz_mod_vec_dot_rev(fmpz_t res, const fmpz_t initial, int subtract, const fmpz * vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            fmpz_zero(res);
        else
            fmpz_set(res, initial);
        return GR_SUCCESS;
    }

    if (initial == NULL)
    {
        fmpz_mul(res, vec1, vec2 + len - 1);
    }
    else
    {
        if (subtract)
            fmpz_neg(res, initial);
        else
            fmpz_set(res, initial);

        fmpz_addmul(res, vec1, vec2 + len - 1);
    }

    for (i = 1; i < len; i++)
        fmpz_addmul(res, vec1 + i, vec2 + len - 1 - i);

    if (subtract)
        fmpz_neg(res, res);

    fmpz_mod(res, res, FMPZ_MOD_CTX(ctx)->n);

    return GR_SUCCESS;
}

int
_gr_fmpz_mod_poly_mullow(fmpz * res,
    const fmpz * poly1, slong len1,
    const fmpz * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    /* weird argument order? */
    if (len1 >= len2)
        _fmpz_mod_poly_mullow(res, poly1, len1, poly2, len2, FMPZ_MOD_CTX(ctx)->n, n);
    else
        _fmpz_mod_poly_mullow(res, poly2, len2, poly1, len1, FMPZ_MOD_CTX(ctx)->n, n);

    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mat_mul(fmpz_mat_t res, const fmpz_mat_t x, const fmpz_mat_t y, gr_ctx_t ctx)
{
    fmpz_mod_mat_t MM;

    fmpz_mat_mul(res, x, y);

    *MM->mat = *res;
    *MM->mod = *(FMPZ_MOD_CTX(ctx)->n);

    _fmpz_mod_mat_reduce(MM);

    return GR_SUCCESS;
}

int _fmpz_mod_methods_initialized = 0;

gr_static_method_table _fmpz_mod_methods;

gr_method_tab_input _fmpz_mod_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fmpz_mod_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fmpz_mod_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_fmpz_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_fmpz_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) _gr_fmpz_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fmpz_mod_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fmpz_mod_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fmpz_mod_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fmpz_mod_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fmpz_mod_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fmpz_mod_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fmpz_mod_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fmpz_mod_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fmpz_mod_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_fmpz_mod_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fmpz_mod_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fmpz_mod_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fmpz_mod_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fmpz_mod_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fmpz_mod_set_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fmpz_mod_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fmpz_mod_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_fmpz_mod_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_fmpz_mod_add_si},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fmpz_mod_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fmpz_mod_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fmpz_mod_mul_si},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_fmpz_mod_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_fmpz_mod_submul},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_fmpz_mod_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fmpz_mod_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fmpz_mod_div},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpz_mod_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpz_mod_inv},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_fmpz_mod_pow_ui},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_fmpz_mod_pow_fmpz},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_fmpz_mod_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_fmpz_mod_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fmpz_mod_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fmpz_mod_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpz_mod(gr_ctx_t ctx, const fmpz_t n)
{
    ctx->which_ring = GR_CTX_FMPZ_MOD;
    ctx->sizeof_elem = sizeof(fmpz);

    ctx->elem_ctx = flint_malloc(sizeof(fmpz_mod_ctx_extended_struct));
    fmpz_mod_ctx_init(FMPZ_MOD_CTX(ctx), n);
    FMPZ_MOD_IS_PRIME(ctx) = T_UNKNOWN;

    ctx->size_limit = WORD_MAX;

    ctx->methods = _fmpz_mod_methods;

    if (!_fmpz_mod_methods_initialized)
    {
        gr_method_tab_init(_fmpz_mod_methods, _fmpz_mod_methods_input);
        _fmpz_mod_methods_initialized = 1;
    }
}

void
gr_ctx_fmpz_mod_set_primality(gr_ctx_t ctx, truth_t is_prime)
{
    FMPZ_MOD_IS_PRIME(ctx) = is_prime;
}
