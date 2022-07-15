#include "qqbar.h"
#include "fexpr.h"
#include "gr.h"

typedef struct
{
    /* todo: use which_ring */
    int real_only;      /* field restricted to real algebraics instead of complex? */
    slong deg_limit;    /* todo */
    slong bits_limit;   /* todo */
}
gr_qqbar_ctx;

#define QQBAR_CTX(ring_ctx) ((gr_qqbar_ctx *)((ring_ctx)->elem_ctx))

int
_gr_qqbar_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    if (QQBAR_CTX(ctx)->real_only)
        gr_stream_write(out, "Real algebraic numbers (qqbar)");
    else
        gr_stream_write(out, "Complex algebraic numbers (qqbar)");
    return GR_SUCCESS;
}

void
_gr_qqbar_init(qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_init(x);
}

void
_gr_qqbar_clear(qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_clear(x);
}

void
_gr_qqbar_swap(qqbar_t x, qqbar_t y, const gr_ctx_t ctx)
{
    qqbar_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
int
_gr_qqbar_randtest(qqbar_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    slong deg_limit, bits_limit;
    int rcase;
    
    rcase = n_randint(state, 10);

    if (rcase == 0)
    {
        deg_limit = 4;
        bits_limit = 10;
    }
    else if (rcase <= 3)
    {
        deg_limit = 2;
        bits_limit = 10;
    }
    else
    {
        deg_limit = 1;
        bits_limit = 10;
    }

    if (QQBAR_CTX(ctx)->real_only)
        qqbar_randtest_real(res, state, deg_limit, bits_limit);
    else
        qqbar_randtest(res, state, deg_limit, bits_limit);

    return GR_SUCCESS;
}

/* todo: different styles */


void
qqbar_get_decimal_root_nearest(char ** re_s, char ** im_s, const qqbar_t x, slong default_digits);

int
_gr_qqbar_write(gr_stream_t out, const qqbar_t x, const gr_ctx_t ctx)
{
    char *re_s, *im_s;

    if (qqbar_is_rational(x))
    {
        fmpq_t t;
        fmpq_init(t);
        qqbar_get_fmpq(t, x);
        gr_stream_write_fmpz(out, fmpq_numref(t));
        if (!fmpz_is_one(fmpq_denref(t)))
        {
            gr_stream_write(out, "/");
            gr_stream_write_fmpz(out, fmpq_denref(t));
        }
        fmpq_clear(t);
    }
    else
    {
        qqbar_get_decimal_root_nearest(&re_s, &im_s, x, 6);

        gr_stream_write(out, "Root a = ");

        if (re_s != NULL)
        {
            gr_stream_write_free(out, re_s);
        }

        if (im_s != NULL)
        {
            gr_stream_write(out, " + ");
            gr_stream_write_free(out, im_s);
            gr_stream_write(out, "i");
        }
        gr_stream_write(out, " of ");
        gr_stream_write_free(out, fmpz_poly_get_str_pretty(QQBAR_POLY(x), "a"));
    }

    return GR_SUCCESS;
}

int
_gr_qqbar_zero(qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_zero(x);
    return GR_SUCCESS;
}

int
_gr_qqbar_one(qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_one(x);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_si(qqbar_t res, slong v, const gr_ctx_t ctx)
{
    qqbar_set_si(res, v);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_ui(qqbar_t res, ulong v, const gr_ctx_t ctx)
{
    qqbar_set_ui(res, v);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_fmpz(qqbar_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    qqbar_set_fmpz(res, v);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_fmpq(qqbar_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    qqbar_set_fmpq(res, v);
    return GR_SUCCESS;
}

truth_t
_gr_qqbar_is_zero(const qqbar_t x, const gr_ctx_t ctx)
{
    return qqbar_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_qqbar_is_one(const qqbar_t x, const gr_ctx_t ctx)
{
    return qqbar_is_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_qqbar_is_neg_one(const qqbar_t x, const gr_ctx_t ctx)
{
    return qqbar_is_neg_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_qqbar_equal(const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    return qqbar_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_qqbar_set(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_set(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_set_other(qqbar_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            qqbar_set_fmpz(res, x);
            return GR_SUCCESS;

        case GR_CTX_FMPQ:
            qqbar_set_fmpq(res, x);
            return GR_SUCCESS;

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_QQBAR && !qqbar_is_real(x))
                return GR_DOMAIN;
            qqbar_set(res, x);
            return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int
_gr_qqbar_neg(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_add(qqbar_t res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    qqbar_add(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_add_si(qqbar_t res, const qqbar_t x, slong y, const gr_ctx_t ctx)
{
    qqbar_add_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_add_ui(qqbar_t res, const qqbar_t x, ulong y, const gr_ctx_t ctx)
{
    qqbar_add_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_add_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    qqbar_add_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_add_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    qqbar_add_fmpq(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub(qqbar_t res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    qqbar_sub(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub_si(qqbar_t res, const qqbar_t x, slong y, const gr_ctx_t ctx)
{
    qqbar_sub_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub_ui(qqbar_t res, const qqbar_t x, ulong y, const gr_ctx_t ctx)
{
    qqbar_sub_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    qqbar_sub_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_sub_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    qqbar_sub_fmpq(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul(qqbar_t res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    qqbar_mul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul_si(qqbar_t res, const qqbar_t x, slong y, const gr_ctx_t ctx)
{
    qqbar_mul_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul_ui(qqbar_t res, const qqbar_t x, ulong y, const gr_ctx_t ctx)
{
    qqbar_mul_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    qqbar_mul_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_mul_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    qqbar_mul_fmpq(res, x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_inv(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (qqbar_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_inv(res, x);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div(qqbar_t res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    if (qqbar_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div_si(qqbar_t res, const qqbar_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div_si(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div_ui(qqbar_t res, const qqbar_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div_ui(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div_fmpz(res, x, y);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_div_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    if (fmpq_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_div_fmpq(res, x, y);
        return GR_SUCCESS;
    }
}

truth_t
_gr_qqbar_is_invertible(const qqbar_t x, const gr_ctx_t ctx)
{
    return !qqbar_is_zero(x) ? T_TRUE : T_FALSE;
}

int
_gr_qqbar_pow_ui(qqbar_t res, const qqbar_t x, ulong exp, const gr_ctx_t ctx)
{
    qqbar_pow_ui(res, x, exp);
    return GR_SUCCESS;
}

int
_gr_qqbar_pow_si(qqbar_t res, const qqbar_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0 && qqbar_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_pow_si(res, x, exp);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_pow_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (fmpz_sgn(exp) < 0 && qqbar_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_pow_fmpz(res, x, exp);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_pow_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    if (fmpq_sgn(exp) < 0 && qqbar_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_pow_fmpq(res, x, exp);

        /* todo: don't compute */
        if (QQBAR_CTX(ctx)->real_only && !qqbar_is_real(res))
        {
            qqbar_zero(res);
            return GR_DOMAIN;
        }
        else
        {
            return GR_SUCCESS;
        }
    }
}

int
_gr_qqbar_pow(qqbar_t res, const qqbar_t x, const qqbar_t exp, const gr_ctx_t ctx)
{
    if (qqbar_pow(res, x, exp))
    {
        if (QQBAR_CTX(ctx)->real_only && !qqbar_is_real(res))
        {
            qqbar_zero(res);
            return GR_DOMAIN;
        }
        else
        {
            return GR_SUCCESS;
        }
    }
    else
    {
        return GR_DOMAIN;
    }
}

truth_t
_gr_qqbar_is_square(const qqbar_t x, const gr_ctx_t ctx)
{
    if (QQBAR_CTX(ctx)->real_only)
        return (qqbar_sgn_re(x) >= 0) ? T_TRUE : T_FALSE;
    else
        return T_TRUE;
}

int
_gr_qqbar_sqrt(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (QQBAR_CTX(ctx)->real_only && qqbar_sgn_re(x) < 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_sqrt(res, x);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_rsqrt(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    if (qqbar_is_zero(x) || (QQBAR_CTX(ctx)->real_only && qqbar_sgn_re(x) < 0))
    {
        return GR_DOMAIN;
    }
    else
    {
        qqbar_rsqrt(res, x);
        return GR_SUCCESS;
    }
}

int
_gr_qqbar_abs(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_abs(res, x);
    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_qqbar_conj(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_conj(res, x);
    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_qqbar_re(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_re(res, x);
    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_qqbar_im(qqbar_t res, const qqbar_t x, const gr_ctx_t ctx)
{
    qqbar_im(res, x);
    return GR_SUCCESS;
}

int
_gr_qqbar_cmp(int * res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    if (!qqbar_is_real(x) || !qqbar_is_real(y))
        return GR_DOMAIN;

    *res = qqbar_cmp_re(x, y);
    return GR_SUCCESS;
}

int
_gr_qqbar_cmpabs(int * res, const qqbar_t x, const qqbar_t y, const gr_ctx_t ctx)
{
    *res = qqbar_cmpabs(x, y);
    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_qqbar_ctx_clear(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
    return GR_SUCCESS;
}

truth_t
_gr_qqbar_ctx_is_algebraically_closed(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_QQBAR) ? T_TRUE : T_FALSE;
}

truth_t
_gr_qqbar_ctx_is_ordered_ring(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_QQBAR) ? T_TRUE : T_FALSE;
}

int _qqbar_methods_initialized = 0;

gr_static_method_table _qqbar_methods;

gr_method_tab_input _qqbar_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_qqbar_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_qqbar_ctx_write},
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
                                (gr_funcptr) _gr_qqbar_ctx_is_algebraically_closed},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) _gr_qqbar_ctx_is_ordered_ring},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_qqbar_init},

    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_qqbar_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_qqbar_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_qqbar_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_qqbar_write},

    {GR_METHOD_ZERO,            (gr_funcptr) _gr_qqbar_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_qqbar_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_qqbar_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_qqbar_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_qqbar_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_qqbar_equal},

    {GR_METHOD_SET,             (gr_funcptr) _gr_qqbar_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_qqbar_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_qqbar_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_qqbar_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_qqbar_set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_qqbar_set_other},

    {GR_METHOD_NEG,             (gr_funcptr) _gr_qqbar_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_qqbar_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_qqbar_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_qqbar_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_qqbar_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_qqbar_add_fmpq},

    {GR_METHOD_SUB,             (gr_funcptr) _gr_qqbar_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_qqbar_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_qqbar_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_qqbar_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_qqbar_sub_fmpq},

    {GR_METHOD_MUL,             (gr_funcptr) _gr_qqbar_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_qqbar_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_qqbar_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_qqbar_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_qqbar_mul_fmpq},

    {GR_METHOD_DIV,             (gr_funcptr) _gr_qqbar_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_qqbar_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_qqbar_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_qqbar_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_qqbar_div_fmpq},

    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_qqbar_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_qqbar_inv},

    {GR_METHOD_POW,             (gr_funcptr) _gr_qqbar_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_qqbar_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_qqbar_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_qqbar_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_qqbar_pow_fmpq},

    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_qqbar_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_qqbar_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_qqbar_rsqrt},

    {GR_METHOD_CMP,             (gr_funcptr) _gr_qqbar_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_qqbar_cmpabs},

    {GR_METHOD_ABS,             (gr_funcptr) _gr_qqbar_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_qqbar_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_qqbar_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_qqbar_im},

    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_real_qqbar(gr_ctx_t ctx)
{
    ctx->flags = 0;
    ctx->which_ring = GR_CTX_REAL_ALGEBRAIC_QQBAR;
    ctx->sizeof_elem = sizeof(qqbar_struct);
    ctx->elem_ctx = flint_malloc(sizeof(gr_qqbar_ctx));
    ctx->size_limit = WORD_MAX;

    QQBAR_CTX(ctx)->real_only = 1;
    QQBAR_CTX(ctx)->deg_limit = WORD_MAX;
    QQBAR_CTX(ctx)->bits_limit = WORD_MAX;

    ctx->methods = _qqbar_methods;

    if (!_qqbar_methods_initialized)
    {
        gr_method_tab_init(_qqbar_methods, _qqbar_methods_input);
        _qqbar_methods_initialized = 1;
    }
}

void
gr_ctx_init_complex_qqbar(gr_ctx_t ctx)
{
    ctx->flags = 0;
    ctx->which_ring = GR_CTX_COMPLEX_ALGEBRAIC_QQBAR;
    ctx->sizeof_elem = sizeof(qqbar_struct);
    ctx->elem_ctx = flint_malloc(sizeof(gr_qqbar_ctx));
    ctx->size_limit = WORD_MAX;

    QQBAR_CTX(ctx)->real_only = 0;
    QQBAR_CTX(ctx)->deg_limit = WORD_MAX;
    QQBAR_CTX(ctx)->bits_limit = WORD_MAX;

    ctx->methods = _qqbar_methods;

    if (!_qqbar_methods_initialized)
    {
        gr_method_tab_init(_qqbar_methods, _qqbar_methods_input);
        _qqbar_methods_initialized = 1;
    }
}
