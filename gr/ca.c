#include "ca.h"
#include "ca_mat.h"
#include "ca_poly.h"
#include "fexpr.h"
#include "gr.h"

#define GR_CA_CTX(ring_ctx) ((ca_ctx_struct *)((ring_ctx)->elem_ctx))

int
_gr_ca_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_WHICH_RING_RR)
        gr_stream_write(out, "Real numbers (ca)");
    else if (ctx->which_ring == GR_WHICH_RING_CC)
        gr_stream_write(out, "Complex numbers (ca)");
    else if (ctx->which_ring == GR_WHICH_RING_RR_ALGEBRAIC)
        gr_stream_write(out, "Real algebraic numbers (ca)");
    else
        gr_stream_write(out, "Complex algebraic numbers (ca)");
    return GR_SUCCESS;
}

void
_gr_ca_init(ca_t x, const gr_ctx_t ctx)
{
    ca_init(x, GR_CA_CTX(ctx));
}

void
_gr_ca_clear(ca_t x, const gr_ctx_t ctx)
{
    ca_clear(x, GR_CA_CTX(ctx));
}

void
_gr_ca_swap(ca_t x, ca_t y, const gr_ctx_t ctx)
{
    ca_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
/* todo: faster real/algebraic constructions */
int
_gr_ca_randtest(ca_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    ca_randtest(res, state, 2, 10, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_WHICH_RING_RR)
    {
        if (ca_check_is_real(res, GR_CA_CTX(ctx)) != T_TRUE)
        {
            ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
        }
    }
    else if (ctx->which_ring == GR_WHICH_RING_RR_ALGEBRAIC)
    {
        if (ca_check_is_real(res, GR_CA_CTX(ctx)) != T_TRUE ||
            ca_check_is_algebraic(res, GR_CA_CTX(ctx)) != T_TRUE)
        {
            ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
        }
    }
    else if (ctx->which_ring == GR_WHICH_RING_CC_ALGEBRAIC)
    {
        if (ca_check_is_algebraic(res, GR_CA_CTX(ctx)) != T_TRUE)
        {
            ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
        }
    }

    return GR_SUCCESS;
}

/* todo */
int
_gr_ca_write(gr_stream_t out, const ca_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, ca_get_str(x, GR_CA_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_ca_zero(ca_t x, const gr_ctx_t ctx)
{
    ca_zero(x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_one(ca_t x, const gr_ctx_t ctx)
{
    ca_one(x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_si(ca_t res, slong v, const gr_ctx_t ctx)
{
    ca_set_si(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_ui(ca_t res, ulong v, const gr_ctx_t ctx)
{
    ca_set_ui(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_fmpz(ca_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    ca_set_fmpz(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_fmpq(ca_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    ca_set_fmpq(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_ca_is_zero(const ca_t x, const gr_ctx_t ctx)
{
    return ca_check_is_zero(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_is_one(const ca_t x, const gr_ctx_t ctx)
{
    return ca_check_is_one(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_is_neg_one(const ca_t x, const gr_ctx_t ctx)
{
    return ca_check_is_neg_one(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_equal(const ca_t x, const ca_t y, const gr_ctx_t ctx)
{
    return ca_check_equal(x, y, GR_CA_CTX(ctx));
}

int
_gr_ca_set(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_set(res, x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_neg(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_neg(res, x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add(ca_t res, const ca_t x, const ca_t y, const gr_ctx_t ctx)
{
    ca_add(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_si(ca_t res, const ca_t x, slong y, const gr_ctx_t ctx)
{
    ca_add_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_ui(ca_t res, const ca_t x, ulong y, const gr_ctx_t ctx)
{
    ca_add_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    ca_add_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    ca_add_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub(ca_t res, const ca_t x, const ca_t y, const gr_ctx_t ctx)
{
    ca_sub(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_si(ca_t res, const ca_t x, slong y, const gr_ctx_t ctx)
{
    ca_sub_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_ui(ca_t res, const ca_t x, ulong y, const gr_ctx_t ctx)
{
    ca_sub_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    ca_sub_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    ca_sub_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul(ca_t res, const ca_t x, const ca_t y, const gr_ctx_t ctx)
{
    ca_mul(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_si(ca_t res, const ca_t x, slong y, const gr_ctx_t ctx)
{
    ca_mul_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_ui(ca_t res, const ca_t x, ulong y, const gr_ctx_t ctx)
{
    ca_mul_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    ca_mul_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    ca_mul_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_inv(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_inv(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_div(ca_t res, const ca_t x, const ca_t y, const gr_ctx_t ctx)
{
    ca_div(res, x, y, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_div_si(ca_t res, const ca_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_div_si(res, x, y, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_ca_div_ui(ca_t res, const ca_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_div_ui(res, x, y, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_ca_div_fmpz(ca_t res, const ca_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_div_fmpz(res, x, y, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    if (fmpq_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_div_fmpq(res, x, y, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
}

truth_t
_gr_ca_is_invertible(const ca_t x, const gr_ctx_t ctx)
{
    return truth_not(ca_check_is_zero(x, GR_CA_CTX(ctx)));
}

int
_gr_ca_pow_ui(ca_t res, const ca_t x, ulong exp, const gr_ctx_t ctx)
{
    ca_pow_ui(res, x, exp, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_pow_si(ca_t res, const ca_t x, slong exp, const gr_ctx_t ctx)
{
    ca_pow_si(res, x, exp, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_pow_fmpz(ca_t res, const ca_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    ca_pow_fmpz(res, x, exp, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_pow_fmpq(ca_t res, const ca_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    ca_pow_fmpq(res, x, exp, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_WHICH_RING_RR || ctx->which_ring == GR_WHICH_RING_RR_ALGEBRAIC)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_pow(ca_t res, const ca_t x, const ca_t exp, const gr_ctx_t ctx)
{
    ca_pow(res, x, exp, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_WHICH_RING_RR || ctx->which_ring == GR_WHICH_RING_RR_ALGEBRAIC)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ctx->which_ring == GR_WHICH_RING_RR_ALGEBRAIC || ctx->which_ring == GR_WHICH_RING_CC_ALGEBRAIC)
    {
        truth_t algebraic;

        algebraic = ca_check_is_algebraic(res, GR_CA_CTX(ctx));

        if (algebraic == T_UNKNOWN)
            return GR_UNABLE;

        if (algebraic == T_FALSE)
            return GR_DOMAIN;
    }

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

truth_t
_gr_ca_is_square(const ca_t x, const gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_WHICH_RING_RR || ctx->which_ring == GR_WHICH_RING_RR_ALGEBRAIC)
    {
        return truth_not(ca_check_is_negative_real(x, GR_CA_CTX(ctx)));
    }
    else
    {
        return T_TRUE;
    }
}

int
_gr_ca_sqrt(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_sqrt(res, x, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_WHICH_RING_RR || ctx->which_ring == GR_WHICH_RING_RR_ALGEBRAIC)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_rsqrt(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_sqrt(res, x, GR_CA_CTX(ctx));
    ca_inv(res, res, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_WHICH_RING_RR || ctx->which_ring == GR_WHICH_RING_RR_ALGEBRAIC)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_abs(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_abs(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_conj(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_conj(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_ca_re(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_re(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_ca_im(ca_t res, const ca_t x, const gr_ctx_t ctx)
{
    ca_im(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_poly_mullow(ca_ptr res,
    ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _ca_poly_mullow(res, poly1, len1, poly2, len2, n, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mat_mul(ca_mat_t res, const ca_mat_t x, const ca_mat_t y, gr_ctx_t ctx)
{
    ca_mat_mul(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_ctx_clear(gr_ctx_t ctx)
{
    ca_ctx_clear(GR_CA_CTX(ctx));
    flint_free(ctx->elem_ctx);
    return GR_SUCCESS;
}

int
_gr_ca_ctx_is_algebraically_closed(gr_ctx_t ctx)
{
    return ctx->which_ring == GR_WHICH_RING_CC || GR_WHICH_RING_CC_ALGEBRAIC;
}

int
_gr_ca_ctx_is_ordered_ring(gr_ctx_t ctx)
{
    return ctx->which_ring == GR_WHICH_RING_RR || GR_WHICH_RING_RR_ALGEBRAIC;
}

int _ca_methods_initialized = 0;

gr_static_method_table _ca_methods;

gr_method_tab_input _ca_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_ca_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_ca_ctx_write},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},

    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) _gr_ca_ctx_is_algebraically_closed},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) _gr_ca_ctx_is_ordered_ring},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},

    {GR_METHOD_INIT,            (gr_funcptr) _gr_ca_init},

    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_ca_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_ca_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_ca_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_ca_write},

    {GR_METHOD_ZERO,            (gr_funcptr) _gr_ca_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_ca_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_ca_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_ca_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_ca_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_ca_equal},

    {GR_METHOD_SET,             (gr_funcptr) _gr_ca_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_ca_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_ca_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_ca_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_ca_set_fmpq},

    {GR_METHOD_NEG,             (gr_funcptr) _gr_ca_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_ca_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_ca_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_ca_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_ca_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_ca_add_fmpq},

    {GR_METHOD_SUB,             (gr_funcptr) _gr_ca_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_ca_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_ca_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_ca_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_ca_sub_fmpq},

    {GR_METHOD_MUL,             (gr_funcptr) _gr_ca_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_ca_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_ca_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_ca_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_ca_mul_fmpq},

    {GR_METHOD_DIV,             (gr_funcptr) _gr_ca_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_ca_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_ca_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_ca_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_ca_div_fmpq},

    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_ca_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_ca_inv},

    {GR_METHOD_POW,             (gr_funcptr) _gr_ca_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_ca_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_ca_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_ca_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_ca_pow_fmpq},

    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_ca_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_ca_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_ca_rsqrt},

    {GR_METHOD_ABS,             (gr_funcptr) _gr_ca_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_ca_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_ca_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_ca_im},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_ca_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_ca_mat_mul},

    {0,                         (gr_funcptr) NULL},
};

void
_gr_ctx_init_ca(gr_ctx_t ctx, int which_ring)
{
    ctx->flags = GR_COMMUTATIVE_RING | GR_FIELD;
    ctx->which_ring = which_ring;
    ctx->sizeof_elem = sizeof(ca_struct);
    ctx->elem_ctx = flint_malloc(sizeof(ca_ctx_struct));
    ctx->size_limit = WORD_MAX;

    ca_ctx_init(GR_CA_CTX(ctx));

    ctx->methods = _ca_methods;

    if (!_ca_methods_initialized)
    {
        gr_method_tab_init(_ca_methods, _ca_methods_input);
        _ca_methods_initialized = 1;
    }
}

void
gr_ctx_init_real_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_WHICH_RING_RR);
}

void
gr_ctx_init_complex_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_WHICH_RING_CC);
}

void
gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_WHICH_RING_RR_ALGEBRAIC);
}

void
gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_WHICH_RING_CC_ALGEBRAIC);
}
