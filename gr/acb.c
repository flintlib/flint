#include "acb.h"
#include "acb_poly.h"
#include "gr.h"

typedef struct
{
    slong prec;
}
gr_acb_ctx;

#define ACB_CTX_PREC(ring_ctx) (((gr_acb_ctx *)((ring_ctx)->elem_ctx))->prec)

int
_gr_acb_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Complex numbers (acb, prec = ");
    gr_stream_write_si(out, ACB_CTX_PREC(ctx));
    gr_stream_write(out, ")");
    return GR_SUCCESS;
}

void
_gr_acb_init(acb_t x, const gr_ctx_t ctx)
{
    acb_init(x);
}

void
_gr_acb_clear(acb_t x, const gr_ctx_t ctx)
{
    acb_clear(x);
}

void
_gr_acb_swap(acb_t x, acb_t y, const gr_ctx_t ctx)
{
    acb_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
int
_gr_acb_randtest(acb_t res, flint_rand_t state, const void * options, const gr_ctx_t ctx)
{
    acb_randtest(res, state, ACB_CTX_PREC(ctx), 10);
    return GR_SUCCESS;
}

/* todo */
int
_gr_acb_write(gr_stream_t out, const acb_t x, const gr_ctx_t ctx)
{
    slong digits = ACB_CTX_PREC(ctx) * 0.30102999566398 + 1;
    int flags = 0;

    if (arb_is_zero(acb_imagref(x)))
    {
        gr_stream_write_free(out, arb_get_str(acb_realref(x), digits, flags));
    }
    else if (arb_is_zero(acb_realref(x)))
    {
        gr_stream_write_free(out, arb_get_str(acb_imagref(x), digits, flags));
        gr_stream_write(out, "*I");
    }
    else
    {
        gr_stream_write_free(out, arb_get_str(acb_realref(x), digits, flags));

        gr_stream_write(out, "(");

        if ((arb_is_exact(acb_imagref(x)) || (flags & ARB_STR_NO_RADIUS))
                && arf_sgn(arb_midref(acb_imagref(x))) < 0)
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, acb_imagref(x));
            gr_stream_write(out, " - ");
            gr_stream_write_free(out, arb_get_str(t, digits, flags));
            arb_clear(t);
        }
        else
        {
            gr_stream_write(out, " + ");
            gr_stream_write_free(out, arb_get_str(acb_imagref(x), digits, flags));
        }

        gr_stream_write(out, "*I)");
    }

    return GR_SUCCESS;
}

int
_gr_acb_zero(acb_t x, const gr_ctx_t ctx)
{
    acb_zero(x);
    return GR_SUCCESS;
}

int
_gr_acb_one(acb_t x, const gr_ctx_t ctx)
{
    acb_one(x);
    return GR_SUCCESS;
}

int
_gr_acb_set_si(acb_t res, slong v, const gr_ctx_t ctx)
{
    acb_set_si(res, v);
    acb_set_round(res, res, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_ui(acb_t res, ulong v, const gr_ctx_t ctx)
{
    acb_set_ui(res, v);
    acb_set_round(res, res, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_fmpz(acb_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    acb_set_round_fmpz(res, v, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_fmpq(acb_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    acb_set_fmpq(res, v, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_acb_is_zero(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
        return T_TRUE;

    if (acb_contains_zero(x))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_is_one(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_one(x))
        return T_TRUE;

    if (arb_contains_zero(acb_imagref(x)) && arb_contains_si(acb_realref(x), 1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_is_neg_one(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_equal_si(x, -1))
        return T_TRUE;

    if (arb_contains_zero(acb_imagref(x)) && arb_contains_si(acb_realref(x), -1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_equal(const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    if (acb_is_exact(x) && acb_equal(x, y))
        return T_TRUE;

    if (acb_overlaps(x, y))
        return T_UNKNOWN;

    return T_FALSE;
}

int
_gr_acb_set(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_set(res, x);
    return GR_SUCCESS;
}

int
_gr_acb_neg(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_acb_add(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_add(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_add_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_add_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_add_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_sub(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_sub_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_sub_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_sub_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_mul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_mul_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_mul_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_mul_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_two(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_mul_2exp_si(res, x, 1);
    return GR_SUCCESS;
}

int
_gr_acb_sqr(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_sqr(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_inv(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_inv(res, x, ACB_CTX_PREC(ctx));
        if (acb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_acb_div(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    if (acb_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div(res, x, y, ACB_CTX_PREC(ctx));

        if (acb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_acb_div_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_si(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_div_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_ui(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_div_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_fmpz(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

truth_t
_gr_acb_is_invertible(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
        return T_FALSE;

    if (acb_contains_zero(x))
        return T_UNKNOWN;

    return T_TRUE;
}

int
_gr_acb_pow_ui(acb_t res, const acb_t x, ulong exp, const gr_ctx_t ctx)
{
    acb_pow_ui(res, x, exp, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_pow_si(acb_t res, const acb_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0 && acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (exp < 0 && acb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        fmpz_t t;
        fmpz_init_set_si(t, exp);
        acb_pow_fmpz(res, x, t, ACB_CTX_PREC(ctx));
        fmpz_clear(t);
        return GR_SUCCESS;
    }
}

int
_gr_acb_pow_fmpz(acb_t res, const acb_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (fmpz_sgn(exp) < 0 && acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (fmpz_sgn(exp) < 0 && acb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        acb_pow_fmpz(res, x, exp, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_pow(acb_t res, const acb_t x, const acb_t exp, const gr_ctx_t ctx)
{
    if (acb_is_int(exp))
    {
        if (arf_sgn(arb_midref(acb_realref(exp))) < 0)
        {
            if (acb_is_zero(x))
                return GR_DOMAIN;

            if (acb_contains_zero(x))
                return GR_UNABLE;
        }

        acb_pow(res, x, exp, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (!acb_contains_zero(x) || arb_is_positive(acb_realref(exp)))
    {
        acb_pow(res, x, exp, ACB_CTX_PREC(ctx));

        if (!acb_is_finite(res))
            return GR_UNABLE;

        return GR_SUCCESS;
    }
    else if (acb_is_zero(x) && arb_is_negative(acb_realref(exp)))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_acb_pow_fmpq(acb_t res, const acb_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    int status;
    acb_t t;
    acb_init(t);
    acb_set_fmpq(t, exp, ACB_CTX_PREC(ctx) + 20);
    status = _gr_acb_pow(res, x, t, ctx);
    acb_clear(t);
    return status;
}

truth_t
_gr_acb_is_square(const acb_t x, const gr_ctx_t ctx)
{
    return T_TRUE;
}

int
_gr_acb_sqrt(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_sqrt(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_rsqrt(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_contains_zero(x))
    {
        acb_rsqrt(res, x, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_acb_vec_dot(acb_t res, const acb_t initial, int subtract, acb_srcptr vec1, acb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    acb_dot(res, initial, subtract, vec1, 1, vec2, 1, len, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_vec_dot_rev(acb_t res, const acb_t initial, int subtract, acb_srcptr vec1, acb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    acb_dot(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_poly_mullow(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _acb_poly_mullow(res, poly1, len1, poly2, len2, n, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_ctx_clear(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
    return GR_SUCCESS;
}

int _acb_methods2_initialized = 0;
gr_static_method_table _acb_static_table;
gr_method_tab_t _acb_methods2;

gr_method_tab_input acb_methods2[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_acb_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_acb_ctx_write},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_acb_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_acb_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_acb_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_acb_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_acb_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_acb_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_acb_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_acb_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_acb_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_acb_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_acb_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_acb_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_acb_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_acb_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_acb_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_acb_set_fmpq},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_acb_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_acb_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_acb_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_acb_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_acb_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_acb_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_acb_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_acb_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_acb_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_acb_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_acb_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_acb_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_acb_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_acb_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_acb_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_acb_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_acb_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_acb_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_acb_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_acb_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_acb_is_invertible},
    {GR_METHOD_POW,             (gr_funcptr) _gr_acb_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_acb_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_acb_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_acb_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_acb_pow_fmpq},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_acb_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_acb_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_acb_rsqrt},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_acb_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_acb_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_acb_poly_mullow},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec)
{
    ctx->flags = GR_COMMUTATIVE_RING | GR_FIELD;
    ctx->which_ring = GR_WHICH_RING_CC;
    ctx->sizeof_elem = sizeof(acb_struct);
    ctx->elem_ctx = flint_malloc(sizeof(gr_acb_ctx));
    ctx->size_limit = WORD_MAX;

    if (prec < 2 || prec > WORD_MAX / 4)
        abort();

    ACB_CTX_PREC(ctx) = prec;

    if (!_acb_methods2_initialized)
    {
        gr_method_tab_init_static(&_acb_methods2, _acb_static_table, acb_methods2);
        _acb_methods2_initialized = 1;
    }

    ctx->methods2 = &_acb_methods2;
}
