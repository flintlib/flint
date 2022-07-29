#include "qqbar.h"
#include "gr.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_mat.h"

int
_gr_fmpz_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integer ring (fmpz)");
    return GR_SUCCESS;
}

void
_gr_fmpz_init(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_init(x);
}

void
_gr_fmpz_clear(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_clear(x);
}

void
_gr_fmpz_swap(fmpz_t x, fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
int
_gr_fmpz_randtest(fmpz_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    switch (n_randint(state, 4))
    {
        case 0:
            fmpz_randtest(res, state, 100);
            break;
        default:
            fmpz_randtest(res, state, 10);
    }

    return GR_SUCCESS;
}

int
_gr_fmpz_write(gr_stream_t out, const fmpz_t x, const gr_ctx_t ctx)
{
    gr_stream_write_fmpz(out, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_zero(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_zero(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_one(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_one(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_set_si(fmpz_t res, slong v, const gr_ctx_t ctx)
{
    fmpz_set_si(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpz_set_ui(fmpz_t res, ulong v, const gr_ctx_t ctx)
{
    fmpz_set_ui(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpz_set_fmpz(fmpz_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fmpz_set(res, v);
    return GR_SUCCESS;
}

int
_gr_fmpz_set_fmpq(fmpz_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    if (fmpz_is_one(fmpq_denref(v)))
    {
        fmpz_set(res, fmpq_numref(v));
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpz_set_other(fmpz_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            fmpz_set(res, x);
            return GR_SUCCESS;

        case GR_CTX_FMPQ:
            return _gr_fmpz_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            if (qqbar_is_integer(x))
            {
                qqbar_get_fmpz(res, x);
                return GR_SUCCESS;
            }
            return GR_DOMAIN;
    }

    return GR_UNABLE;
}

int
_gr_fmpz_set_str(fmpz_t res, const char * x, const gr_ctx_t ctx)
{
    if (fmpz_set_str(res, x, 10))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_fmpz_get_ui(ulong * res, const fmpz_t x, const gr_ctx_t ctx)
{
    if (fmpz_sgn(x) < 0 || fmpz_cmp_ui(x, UWORD_MAX) > 0)
        return GR_DOMAIN;

    *res = fmpz_get_ui(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_get_si(slong * res, const fmpz_t x, const gr_ctx_t ctx)
{
    if (!fmpz_fits_si(x))
        return GR_DOMAIN;

    *res = fmpz_get_si(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_get_d(double * res, const fmpz_t x, const gr_ctx_t ctx)
{
    *res = fmpz_get_d(x);
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_is_zero(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_is_one(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_is_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_is_neg_one(const fmpz_t x, const gr_ctx_t ctx)
{
    return (*x == -1) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_equal(const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    return fmpz_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_set(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_set(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_neg(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_add(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_add(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_add_si(fmpz_t res, const fmpz_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_add_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_add_ui(fmpz_t res, const fmpz_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_add_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_sub(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_sub(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_sub_si(fmpz_t res, const fmpz_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_sub_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_sub_ui(fmpz_t res, const fmpz_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_sub_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_mul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_mul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_mul_ui(fmpz_t res, const fmpz_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_mul_ui(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_mul_si(fmpz_t res, const fmpz_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_mul_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_addmul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_addmul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_submul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_submul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_mul_two(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_mul_2exp(res, x, 1);
    return GR_SUCCESS;
}

int
_gr_fmpz_sqr(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_mul(res, x, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_inv(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    if (fmpz_is_pm1(x))
    {
        fmpz_set(res, x);
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpz_div(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        if (fmpz_divides(res, x, y))
            return GR_SUCCESS;
        else
            return GR_DOMAIN;
    }
}

truth_t
_gr_fmpz_is_invertible(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_is_pm1(x) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_pow_ui(fmpz_t res, const fmpz_t x, ulong exp, const gr_ctx_t ctx)
{
    if (exp > (ulong) WORD_MAX || exp >= ctx->size_limit)  /* todo: systematic size solution for test code */
    {
        return GR_UNABLE;
    }
    else
    {
        fmpz_pow_ui(res, x, exp);
        return GR_SUCCESS;
    }
}

int
_gr_fmpz_pow_si(fmpz_t res, const fmpz_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0)
    {
        if (exp == -1 && fmpz_is_pm1(x))
        {
            fmpz_set(res, x);
            return GR_SUCCESS;
        }

        return GR_DOMAIN;
    }
    else
    {
        return _gr_fmpz_pow_ui(res, x, exp, ctx);
    }
}

int
_gr_fmpz_pow_fmpz(fmpz_t res, const fmpz_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*exp))
    {
        return _gr_fmpz_pow_si(res, x, *exp, ctx);
    }
    else if (fmpz_is_pm1(x))
    {
        if (fmpz_is_one(x) || fmpz_is_even(exp))
            fmpz_one(res);
        else
            fmpz_set_si(res, -1);
        return GR_SUCCESS;
    }
    else if (fmpz_is_zero(x) && fmpz_sgn(exp) > 0)
    {
        fmpz_zero(res);
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

truth_t
_gr_fmpz_is_square(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_is_square(x) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_sqrt(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    if (fmpz_sgn(x) < 0)
        return GR_DOMAIN;

    if (fmpz_root(res, x, 2))
    {
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpz_rsqrt(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    if (fmpz_is_one(x))
    {
        fmpz_one(res);
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_fmpz_abs(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_abs(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_im(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_zero(res);
    return GR_SUCCESS;
}

int
_gr_fmpz_sgn(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_set_si(res, fmpz_sgn(x));
    return GR_SUCCESS;
}

int
_gr_fmpz_cmp(int * res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    int cmp = fmpz_cmp(x, y);

    if (cmp < 0) cmp = -1;
    if (cmp > 0) cmp = 1;
    *res = cmp;
    return GR_SUCCESS;
}

int
_gr_fmpz_cmpabs(int * res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    int cmp = fmpz_cmpabs(x, y);

    if (cmp < 0) cmp = -1;
    if (cmp > 0) cmp = 1;
    *res = cmp;
    return GR_SUCCESS;
}

int
_gr_fmpz_vec_dot(fmpz_t res, const fmpz_t initial, int subtract, const fmpz * vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)
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

    return GR_SUCCESS;
}

int
_gr_fmpz_vec_dot_rev(fmpz_t res, const fmpz_t initial, int subtract, const fmpz * vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)
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

    return GR_SUCCESS;
}

int
_gr_fmpz_poly_mullow(fmpz * res,
    const fmpz * poly1, slong len1,
    const fmpz * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    if (len1 >= len2)
        _fmpz_poly_mullow(res, poly1, len1, poly2, len2, n);
    else
        _fmpz_poly_mullow(res, poly2, len2, poly1, len1, n);

    return GR_SUCCESS;
}

int
_gr_fmpz_mat_mul(fmpz_mat_t res, const fmpz_mat_t x, const fmpz_mat_t y, gr_ctx_t ctx)
{
    fmpz_mat_mul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fmpz_mat_det(fmpz_t res, const fmpz_mat_t x, const gr_ctx_t ctx)
{
    fmpz_mat_det(res, x);
    return GR_SUCCESS;
}


int _fmpz_methods_initialized = 0;

gr_static_method_table _fmpz_methods;

gr_method_tab_input _fmpz_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fmpz_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fmpz_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fmpz_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fmpz_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fmpz_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fmpz_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fmpz_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fmpz_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fmpz_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fmpz_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_fmpz_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fmpz_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fmpz_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fmpz_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fmpz_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fmpz_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_fmpz_set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_fmpz_set_other},
    {GR_METHOD_SET_STR,         (gr_funcptr) _gr_fmpz_set_str},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_fmpz_set},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_fmpz_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_fmpz_get_si},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_fmpz_get_d},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fmpz_neg},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fmpz_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fmpz_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_fmpz_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_fmpz_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_fmpz_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fmpz_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_fmpz_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_fmpz_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_fmpz_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fmpz_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fmpz_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fmpz_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fmpz_mul},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_fmpz_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_fmpz_submul},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_fmpz_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fmpz_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fmpz_div},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpz_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpz_inv},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_fmpz_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_fmpz_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_fmpz_pow_fmpz},
    {GR_METHOD_POW,             (gr_funcptr) _gr_fmpz_pow_fmpz},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fmpz_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fmpz_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_fmpz_rsqrt},
    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_fmpz_set},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_fmpz_set},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_fmpz_set},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_fmpz_set},
    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_PI,              (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_fmpz_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_fmpz_set},
    {GR_METHOD_RE,              (gr_funcptr) _gr_fmpz_set},
    {GR_METHOD_IM,              (gr_funcptr) _gr_fmpz_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_fmpz_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_fmpz_sgn},
    {GR_METHOD_CMP,             (gr_funcptr) _gr_fmpz_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_fmpz_cmpabs},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_fmpz_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_fmpz_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fmpz_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fmpz_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_fmpz_mat_det},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpz(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_FMPZ;
    ctx->sizeof_elem = sizeof(fmpz);
    ctx->elem_ctx = NULL;
    ctx->size_limit = WORD_MAX;

    ctx->methods = _fmpz_methods;

    if (!_fmpz_methods_initialized)
    {
        gr_method_tab_init(_fmpz_methods, _fmpz_methods_input);
        _fmpz_methods_initialized = 1;
    }
}
