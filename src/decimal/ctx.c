/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "decimal.h"
#include "mag.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_mat.h"

/* Not performance-critical: optimize for size. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/* ------------------------------------------------------------------------- */
/*    Context                                                                */
/* ------------------------------------------------------------------------- */

void
decimal_ctx_clear(gr_ctx_t ctx)
{
    radix_clear(DECIMAL_CTX_RADIX(ctx));
    flint_free(DECIMAL_CTX(ctx));
}

static const char * rnd_names[] = { "down", "up", "floor", "ceil", "near", "near_away", "near_zero" };

int
decimal_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    switch (DECIMAL_CTX_WHICH(ctx))
    {
        case DECIMAL_CTX_FLOAT: status |= gr_stream_write(out, "Decimal floating-point numbers"); break;
        case DECIMAL_CTX_BALL: status |= gr_stream_write(out, "Decimal balls"); break;
        case DECIMAL_CTX_CFLOAT: status |= gr_stream_write(out, "Complex decimal floating-point numbers"); break;
        default: status |= gr_stream_write(out, "Complex decimal balls"); break;
    }

    status |= gr_stream_write(out, " (prec ");
    if (DECIMAL_CTX_IS_EXACT(ctx))
        status |= gr_stream_write(out, "exact");
    else
        status |= gr_stream_write_si(out, DECIMAL_CTX_PREC(ctx));

    if (!DECIMAL_CTX_IS_BALL(ctx))
    {
        status |= gr_stream_write(out, ", rnd ");
        status |= gr_stream_write(out, rnd_names[DECIMAL_CTX_RND(ctx)]);
        if (DECIMAL_CTX_IS_COMPLEX(ctx) && DECIMAL_CTX_RND_IM(ctx) != DECIMAL_CTX_RND(ctx))
        {
            status |= gr_stream_write(out, ", rnd im ");
            status |= gr_stream_write(out, rnd_names[DECIMAL_CTX_RND_IM(ctx)]);
        }
    }
    else
    {
        status |= gr_stream_write(out, ", rad prec ");
        status |= gr_stream_write_si(out, DECIMAL_CTX_RAD_PREC(ctx));
    }

    if (DECIMAL_CTX_E(ctx) != ((FLINT_BITS == 64) ? 19 : 9))
    {
        status |= gr_stream_write(out, ", limb 10^");
        status |= gr_stream_write_si(out, DECIMAL_CTX_E(ctx));
    }

    if (DECIMAL_CTX_HAS_EXP_LIMITS(ctx))
    {
        status |= gr_stream_write(out, ", exp [");
        status |= gr_stream_write_si(out, DECIMAL_CTX_EMIN(ctx));
        status |= gr_stream_write(out, ", ");
        status |= gr_stream_write_si(out, DECIMAL_CTX_EMAX(ctx));
        status |= gr_stream_write(out, "]");
    }

    if (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_ALLOW_INF)
        status |= gr_stream_write(out, ", inf");
    if (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_ALLOW_NAN)
        status |= gr_stream_write(out, ", nan");
    if (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_ALLOW_UNDERFLOW)
        status |= gr_stream_write(out, ", underflow");
    if (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_SLOPPY_RADIUS)
        status |= gr_stream_write(out, ", sloppy radius");

    status |= gr_stream_write(out, ")");
    return status;
}

void
decimal_ctx_set_prec(gr_ctx_t ctx, slong prec)
{
    if (prec != DECIMAL_PREC_EXACT)
    {
        prec = FLINT_MAX(prec, 1);
        prec = FLINT_MIN(prec, DECIMAL_PREC_MAX);
    }
    DECIMAL_CTX_PREC(ctx) = prec;
}

void
decimal_ctx_set_rnd(gr_ctx_t ctx, int rnd)
{
    if (rnd < 0 || rnd >= DECIMAL_RND_NUM)
        flint_throw(FLINT_ERROR, "invalid rounding mode\n");
    DECIMAL_CTX_RND(ctx) = rnd;
    DECIMAL_CTX_RND_IM(ctx) = rnd;
}

void
decimal_ctx_set_rnd_im(gr_ctx_t ctx, int rnd)
{
    if (rnd < 0 || rnd >= DECIMAL_RND_NUM)
        flint_throw(FLINT_ERROR, "invalid rounding mode\n");
    DECIMAL_CTX_RND_IM(ctx) = rnd;
}

void
decimal_ctx_set_rad_prec(gr_ctx_t ctx, slong rad_prec)
{
    slong i;
    ulong p;

    rad_prec = FLINT_MAX(rad_prec, DECMAG_MIN_PREC);
    rad_prec = FLINT_MIN(rad_prec, DECMAG_MAX_PREC);

    DECIMAL_CTX_RAD_PREC(ctx) = rad_prec;
    for (p = 1, i = 0; i < rad_prec; i++)
        p *= 10;
    DECIMAL_CTX_RAD_POW(ctx) = p;
    DECIMAL_CTX_RAD_POW1(ctx) = p / 10;
}

void
decimal_ctx_set_exp_limits(gr_ctx_t ctx, slong emin, slong emax)
{
    DECIMAL_CTX_EMIN(ctx) = emin;
    DECIMAL_CTX_EMAX(ctx) = emax;
}

void
decimal_ctx_set_flags(gr_ctx_t ctx, int flags)
{
    DECIMAL_CTX_FLAGS(ctx) = flags;
}

slong decimal_ctx_get_prec(gr_ctx_t ctx) { return DECIMAL_CTX_PREC(ctx); }
int decimal_ctx_get_rnd(gr_ctx_t ctx) { return DECIMAL_CTX_RND(ctx); }
int decimal_ctx_get_rnd_im(gr_ctx_t ctx) { return DECIMAL_CTX_RND_IM(ctx); }
slong decimal_ctx_get_rad_prec(gr_ctx_t ctx) { return DECIMAL_CTX_RAD_PREC(ctx); }
int decimal_ctx_get_flags(gr_ctx_t ctx) { return DECIMAL_CTX_FLAGS(ctx); }
slong decimal_ctx_get_limb_digits(gr_ctx_t ctx) { return DECIMAL_CTX_E(ctx); }

void
decimal_ctx_get_exp_limits(slong * emin, slong * emax, gr_ctx_t ctx)
{
    *emin = DECIMAL_CTX_EMIN(ctx);
    *emax = DECIMAL_CTX_EMAX(ctx);
}

int
_decimal_ctx_set_real_prec(gr_ctx_t ctx, slong prec_bits)
{
    slong prec;
    prec_bits = FLINT_MAX(prec_bits, 2);
    prec = (slong) ceil(prec_bits * 0.30102999566398119521);
    decimal_ctx_set_prec(ctx, prec);
    return GR_SUCCESS;
}

int
_decimal_ctx_get_real_prec(slong * res, gr_ctx_t ctx)
{
    if (DECIMAL_CTX_IS_EXACT(ctx))
        return GR_UNABLE;
    /* directed rounding to prec digits has relative error up to 10^(1-prec),
       so report slightly fewer bits than log2(10) * prec */
    *res = (slong) (DECIMAL_CTX_PREC(ctx) * 3.3219280948873623479) - 1;
    *res = FLINT_MAX(*res, 2);
    return GR_SUCCESS;
}

static truth_t
_decimal_ctx_is_exact(gr_ctx_t ctx)
{
    return DECIMAL_CTX_IS_EXACT(ctx) ? T_TRUE : T_FALSE;
}

static truth_t
_decimal_ctx_is_ring(gr_ctx_t ctx)
{
    /* the exact ring Z[1/10], provided that no special values are admitted
       and that no results are flushed to zero */
    return (DECIMAL_CTX_IS_EXACT(ctx) && !(DECIMAL_CTX_FLAGS(ctx) & (DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN | DECIMAL_ALLOW_UNDERFLOW))) ? T_TRUE : T_FALSE;
}

static truth_t
_decimal_ctx_is_field(gr_ctx_t ctx)
{
    return T_FALSE;
}

static truth_t
_decimal_ctx_is_approx(gr_ctx_t ctx)
{
    return (_decimal_ctx_is_ring(ctx) == T_TRUE) ? T_FALSE : T_TRUE;
}

/* with rounding, the structure behaves as a real vector space; the exact
   structure only contains decimal fractions */
truth_t
_decimal_ctx_is_vector_space(gr_ctx_t ctx)
{
    return DECIMAL_CTX_IS_EXACT(ctx) ? T_FALSE : T_TRUE;
}

static truth_t
_decimal_ctx_is_canonical(gr_ctx_t ctx)
{
    return (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_ALLOW_NAN) ? T_FALSE : T_TRUE;
}

/* ------------------------------------------------------------------------- */
/*    Method table for decfloat                                              */
/* ------------------------------------------------------------------------- */


int _decfloat_methods_initialized = 0;
gr_static_method_table _decfloat_methods;

gr_method_tab_input _decfloat_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) decimal_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) decimal_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) _decimal_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) _decimal_ctx_is_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _decimal_ctx_is_ring},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _decimal_ctx_is_field},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) _decimal_ctx_is_ring},
    {GR_METHOD_CTX_IS_FINITE,   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING, (gr_funcptr) _decimal_ctx_is_ring},
    {GR_METHOD_CTX_IS_APPROX_COMMUTATIVE_RING, (gr_funcptr) _decimal_ctx_is_approx},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) _decimal_ctx_is_exact},
    {GR_METHOD_CTX_IS_CANONICAL, (gr_funcptr) _decimal_ctx_is_canonical},
    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _decimal_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _decimal_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) decfloat_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) decfloat_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) decfloat_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) decfloat_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) decfloat_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) decfloat_write},
    {GR_METHOD_ZERO,            (gr_funcptr) decfloat_zero},
    {GR_METHOD_ONE,             (gr_funcptr) decfloat_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) decfloat_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) decfloat_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) decfloat_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) decfloat_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) decfloat_equal},
    {GR_METHOD_SET,             (gr_funcptr) decfloat_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) decfloat_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) decfloat_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) decfloat_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) decfloat_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) decfloat_set_d},
    {GR_METHOD_SET_STR,         (gr_funcptr) decfloat_set_str},
    {GR_METHOD_SET_FMPZ_10EXP_FMPZ, (gr_funcptr) decfloat_set_fmpz_10exp_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) decfloat_set_other},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) decfloat_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) decfloat_get_fmpq},
    {GR_METHOD_GET_UI,          (gr_funcptr) decfloat_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) decfloat_get_si},
    {GR_METHOD_GET_D,           (gr_funcptr) decfloat_get_d},

    {GR_METHOD_NEG,             (gr_funcptr) decfloat_neg},
    {GR_METHOD_ADD,             (gr_funcptr) decfloat_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) decfloat_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) decfloat_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) decfloat_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) decfloat_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) decfloat_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) decfloat_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) decfloat_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) decfloat_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) decfloat_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) decfloat_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) decfloat_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) decfloat_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) decfloat_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) decfloat_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) decfloat_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) decfloat_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) decfloat_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) decfloat_inv},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) decfloat_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,   (gr_funcptr) decfloat_mul_2exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ, (gr_funcptr) decfloat_set_fmpz_2exp_fmpz},
    {GR_METHOD_POW_UI,          (gr_funcptr) decfloat_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) decfloat_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) decfloat_pow_fmpz},
    {GR_METHOD_SQRT,            (gr_funcptr) decfloat_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) decfloat_rsqrt},
    {GR_METHOD_POS_INF,         (gr_funcptr) decfloat_pos_inf},
    {GR_METHOD_NEG_INF,         (gr_funcptr) decfloat_neg_inf},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) decfloat_nan},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) decfloat_nan},
    {GR_METHOD_FLOOR,           (gr_funcptr) decfloat_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) decfloat_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) decfloat_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) decfloat_nint},
    {GR_METHOD_ABS,             (gr_funcptr) decfloat_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) decfloat_set},
    {GR_METHOD_RE,              (gr_funcptr) decfloat_set},
    {GR_METHOD_IM,              (gr_funcptr) decfloat_zero},
    {GR_METHOD_SGN,             (gr_funcptr) decfloat_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) decfloat_sgn},
    {GR_METHOD_CMP,             (gr_funcptr) decfloat_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) decfloat_cmpabs},
    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_ATAN2, (gr_funcptr) decfloat_atan2},
    {GR_METHOD_POW,             (gr_funcptr) decfloat_pow},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) decfloat_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) decfloat_vec_dot_rev},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT, (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},
    {0,                         (gr_funcptr) NULL},
};

/* the other tables live in ball.c, cfloat.c and cball.c */
extern int _decball_methods_initialized;
extern gr_static_method_table _decball_methods;
extern gr_method_tab_input _decball_methods_input[];
extern int _deccfloat_methods_initialized;
extern gr_static_method_table _deccfloat_methods;
extern gr_method_tab_input _deccfloat_methods_input[];
extern int _deccball_methods_initialized;
extern gr_static_method_table _deccball_methods;
extern gr_method_tab_input _deccball_methods_input[];

/* the type-specific methods, followed by the special functions shared by
   all types (cfunctions.c) */
extern gr_method_tab_input _decimal_special_methods_input[];
extern gr_method_tab_input _decimal_complex_special_methods_input[];

static void
_decimal_method_tab_init(gr_funcptr * tab, gr_method_tab_input * input, int complex)
{
    gr_method_tab_input * all;
    slong i, n = 0;

    for (i = 0; input[i].function != NULL; i++) n++;
    for (i = 0; _decimal_special_methods_input[i].function != NULL; i++) n++;
    for (i = 0; _decimal_complex_special_methods_input[i].function != NULL; i++) n++;

    all = flint_malloc(sizeof(gr_method_tab_input) * (n + 1));
    n = 0;
    for (i = 0; input[i].function != NULL; i++) all[n++] = input[i];
    for (i = 0; _decimal_special_methods_input[i].function != NULL; i++) all[n++] = _decimal_special_methods_input[i];
    if (complex)
        for (i = 0; _decimal_complex_special_methods_input[i].function != NULL; i++) all[n++] = _decimal_complex_special_methods_input[i];
    all[n].index = 0;
    all[n].function = NULL;

    gr_method_tab_init(tab, all);
    flint_free(all);
}

void
_gr_ctx_init_decimal(gr_ctx_t ctx, decimal_ctx_which which, unsigned int e, slong prec, int rnd, int flags)
{
    decimal_ctx_struct * dctx;

    switch (which)
    {
        case DECIMAL_CTX_FLOAT:
            ctx->which_ring = GR_CTX_DECFLOAT;
            ctx->sizeof_elem = sizeof(decfloat_struct);
            break;
        case DECIMAL_CTX_BALL:
            ctx->which_ring = GR_CTX_DECBALL;
            ctx->sizeof_elem = sizeof(decball_struct);
            break;
        case DECIMAL_CTX_CFLOAT:
            ctx->which_ring = GR_CTX_DECCFLOAT;
            ctx->sizeof_elem = sizeof(deccfloat_struct);
            break;
        default:
            ctx->which_ring = GR_CTX_DECCBALL;
            ctx->sizeof_elem = sizeof(deccball_struct);
            break;
    }
    ctx->size_limit = WORD_MAX;

    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(decimal_ctx_struct));
    dctx = DECIMAL_CTX(ctx);

    radix_init(&dctx->radix, 10, e);
    dctx->which = which;
    dctx->flags = flags;
    dctx->emin = WORD_MIN;
    dctx->emax = WORD_MAX;
    decimal_ctx_set_prec(ctx, prec);
    decimal_ctx_set_rnd(ctx, rnd);
    decimal_ctx_set_rad_prec(ctx, DECMAG_DEFAULT_PREC);

#define INIT_METHODS(tab, input, initialized) \
    do { \
        ctx->methods = tab; \
        if (!initialized) \
        { \
            _decimal_method_tab_init(tab, input, which & 2); \
            initialized = 1; \
        } \
    } while (0)

    switch (which)
    {
        case DECIMAL_CTX_FLOAT: INIT_METHODS(_decfloat_methods, _decfloat_methods_input, _decfloat_methods_initialized); break;
        case DECIMAL_CTX_BALL: INIT_METHODS(_decball_methods, _decball_methods_input, _decball_methods_initialized); break;
        case DECIMAL_CTX_CFLOAT: INIT_METHODS(_deccfloat_methods, _deccfloat_methods_input, _deccfloat_methods_initialized); break;
        default: INIT_METHODS(_deccball_methods, _deccball_methods_input, _deccball_methods_initialized); break;
    }
}

void
gr_ctx_init_decfloat(gr_ctx_t ctx, slong prec, int flags)
{
    _gr_ctx_init_decimal(ctx, DECIMAL_CTX_FLOAT, 0, prec, DECIMAL_RND_NEAR, flags);
}

void
gr_ctx_init_decball(gr_ctx_t ctx, slong prec, int flags)
{
    _gr_ctx_init_decimal(ctx, DECIMAL_CTX_BALL, 0, prec, DECIMAL_RND_NEAR, flags);
}

void
gr_ctx_init_deccfloat(gr_ctx_t ctx, slong prec, int flags)
{
    _gr_ctx_init_decimal(ctx, DECIMAL_CTX_CFLOAT, 0, prec, DECIMAL_RND_NEAR, flags);
}

void
gr_ctx_init_deccball(gr_ctx_t ctx, slong prec, int flags)
{
    _gr_ctx_init_decimal(ctx, DECIMAL_CTX_CBALL, 0, prec, DECIMAL_RND_NEAR, flags);
}

static void
_gr_ctx_init_decimal_randtest(gr_ctx_t ctx, decimal_ctx_which which, flint_rand_t state, slong max_prec)
{
    unsigned int e;
    slong prec;
    int rnd, flags;

    switch (n_randint(state, 4))
    {
        case 0: e = 0; break;
        case 1: e = 1 + n_randint(state, 3); break;
        default: e = 1 + n_randint(state, (FLINT_BITS == 64) ? 19 : 9); break;
    }

    if (n_randint(state, 8) == 0 && !(which & 1))
        prec = DECIMAL_PREC_EXACT;
    else
        prec = 1 + n_randint(state, max_prec);

    rnd = n_randint(state, DECIMAL_RND_NUM);

    flags = 0;
    if (n_randint(state, 2)) flags |= DECIMAL_ALLOW_INF;
    if (n_randint(state, 2)) flags |= DECIMAL_ALLOW_NAN;
    if (n_randint(state, 2)) flags |= DECIMAL_ALLOW_UNDERFLOW;
    if (n_randint(state, 2)) flags |= DECIMAL_SLOPPY_RADIUS;
    if (n_randint(state, 4) == 0) flags |= DECIMAL_WRITE_SCIENTIFIC;

    _gr_ctx_init_decimal(ctx, which, e, prec, rnd, flags);

    if ((which & 2) && n_randint(state, 2))
        decimal_ctx_set_rnd_im(ctx, n_randint(state, DECIMAL_RND_NUM));

    decimal_ctx_set_rad_prec(ctx, DECMAG_MIN_PREC + n_randint(state, DECMAG_MAX_PREC - DECMAG_MIN_PREC + 1));

    if (n_randint(state, 4) == 0)
    {
        slong emax = n_randint(state, 200);
        slong emin = -n_randint(state, 200);
        decimal_ctx_set_exp_limits(ctx, emin, emax);
    }
}

void
gr_ctx_init_decfloat_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec)
{
    _gr_ctx_init_decimal_randtest(ctx, DECIMAL_CTX_FLOAT, state, max_prec);
}

void
gr_ctx_init_decball_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec)
{
    _gr_ctx_init_decimal_randtest(ctx, DECIMAL_CTX_BALL, state, max_prec);
}

void
gr_ctx_init_deccfloat_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec)
{
    _gr_ctx_init_decimal_randtest(ctx, DECIMAL_CTX_CFLOAT, state, max_prec);
}

void
gr_ctx_init_deccball_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec)
{
    _gr_ctx_init_decimal_randtest(ctx, DECIMAL_CTX_CBALL, state, max_prec);
}

POP_OPTIONS
