/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nfloat.h"
#include "gr.h"
#include "gr_mat.h"

int _nfloat_methods_initialized = 0;

gr_static_method_table _nfloat_methods;

gr_method_tab_input _nfloat_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) nfloat_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},

    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _nfloat_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _nfloat_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) nfloat_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) nfloat_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) nfloat_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) nfloat_set},
    {GR_METHOD_RANDTEST,        (gr_funcptr) nfloat_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) nfloat_write},
    {GR_METHOD_ZERO,            (gr_funcptr) nfloat_zero},
    {GR_METHOD_ONE,             (gr_funcptr) nfloat_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) nfloat_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) nfloat_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) nfloat_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) nfloat_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) nfloat_equal},
    {GR_METHOD_SET,             (gr_funcptr) nfloat_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) nfloat_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) nfloat_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) nfloat_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) nfloat_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) nfloat_set_d},
    {GR_METHOD_SET_STR,         (gr_funcptr) nfloat_set_str},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) nfloat_set_other},
/*
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) nfloat_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) nfloat_get_fmpq},
    {GR_METHOD_GET_UI,          (gr_funcptr) nfloat_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) nfloat_get_si},
    {GR_METHOD_GET_D,           (gr_funcptr) nfloat_get_d},
*/

    {GR_METHOD_NEG,             (gr_funcptr) nfloat_neg},
    {GR_METHOD_ADD,             (gr_funcptr) nfloat_add},
/*
    {GR_METHOD_ADD_UI,          (gr_funcptr) nfloat_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) nfloat_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) nfloat_add_fmpz},
*/
    {GR_METHOD_SUB,             (gr_funcptr) nfloat_sub},
/*
    {GR_METHOD_SUB_UI,          (gr_funcptr) nfloat_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) nfloat_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) nfloat_sub_fmpz},
*/
    {GR_METHOD_MUL,             (gr_funcptr) nfloat_mul},
/*
    {GR_METHOD_MUL_UI,          (gr_funcptr) nfloat_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) nfloat_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) nfloat_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) nfloat_mul_two},
*/
    {GR_METHOD_ADDMUL,          (gr_funcptr) nfloat_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) nfloat_submul},
    {GR_METHOD_SQR,             (gr_funcptr) nfloat_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) nfloat_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) nfloat_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) nfloat_div_si},
/*
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) nfloat_div_fmpz},
*/
    {GR_METHOD_INV,             (gr_funcptr) nfloat_inv},

    {GR_METHOD_MUL_2EXP_SI,        (gr_funcptr) nfloat_mul_2exp_si},
/*
    {GR_METHOD_MUL_2EXP_FMPZ,      (gr_funcptr) nfloat_mul_2exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ, (gr_funcptr) nfloat_set_fmpz_2exp_fmpz},
    {GR_METHOD_GET_FMPZ_2EXP_FMPZ, (gr_funcptr) nfloat_get_fmpz_2exp_fmpz},
*/

    {GR_METHOD_POW,             (gr_funcptr) nfloat_pow},

/*
    {GR_METHOD_POW_UI,          (gr_funcptr) nfloat_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) nfloat_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) nfloat_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) nfloat_pow_fmpq},
*/

    {GR_METHOD_SQRT,            (gr_funcptr) nfloat_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) nfloat_rsqrt},

    {GR_METHOD_POS_INF,         (gr_funcptr) nfloat_pos_inf},
    {GR_METHOD_NEG_INF,         (gr_funcptr) nfloat_neg_inf},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) nfloat_nan},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) nfloat_nan},

    {GR_METHOD_FLOOR,           (gr_funcptr) nfloat_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) nfloat_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) nfloat_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) nfloat_nint},

    {GR_METHOD_ABS,             (gr_funcptr) nfloat_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) nfloat_set},
    {GR_METHOD_RE,              (gr_funcptr) nfloat_set},
    {GR_METHOD_IM,              (gr_funcptr) nfloat_im},
    {GR_METHOD_SGN,             (gr_funcptr) nfloat_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) nfloat_sgn},
    {GR_METHOD_CMP,             (gr_funcptr) nfloat_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) nfloat_cmpabs},

    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},

    {GR_METHOD_PI,              (gr_funcptr) nfloat_pi},
    {GR_METHOD_EXP,             (gr_funcptr) nfloat_exp},
    {GR_METHOD_EXPM1,           (gr_funcptr) nfloat_expm1},
    {GR_METHOD_LOG,             (gr_funcptr) nfloat_log},
    {GR_METHOD_LOG1P,           (gr_funcptr) nfloat_log1p},
    {GR_METHOD_SIN,             (gr_funcptr) nfloat_sin},
    {GR_METHOD_COS,             (gr_funcptr) nfloat_cos},
    {GR_METHOD_TAN,             (gr_funcptr) nfloat_tan},
    {GR_METHOD_SINH,            (gr_funcptr) nfloat_sinh},
    {GR_METHOD_COSH,            (gr_funcptr) nfloat_cosh},
    {GR_METHOD_TANH,            (gr_funcptr) nfloat_tanh},
    {GR_METHOD_ATAN,            (gr_funcptr) nfloat_atan},
    {GR_METHOD_GAMMA,            (gr_funcptr) nfloat_gamma},
    {GR_METHOD_ZETA,             (gr_funcptr) nfloat_zeta},

    {GR_METHOD_VEC_INIT,        (gr_funcptr) _nfloat_vec_init},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) _nfloat_vec_clear},
    {GR_METHOD_VEC_SET,         (gr_funcptr) _nfloat_vec_set},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) _nfloat_vec_zero},
    {GR_METHOD_VEC_ADD,                 (gr_funcptr) _nfloat_vec_add},
    {GR_METHOD_VEC_SUB,                 (gr_funcptr) _nfloat_vec_sub},
    {GR_METHOD_VEC_MUL,                 (gr_funcptr) _nfloat_vec_mul},
    {GR_METHOD_VEC_MUL_SCALAR,          (gr_funcptr) _nfloat_vec_mul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR,          (gr_funcptr) _nfloat_vec_addmul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR,          (gr_funcptr) _nfloat_vec_submul_scalar},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _nfloat_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _nfloat_vec_dot_rev},
/*
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) nfloat_poly_mullow},
    {GR_METHOD_POLY_ROOTS_OTHER,(gr_funcptr) nfloat_poly_roots_other},
*/
    {GR_METHOD_MAT_MUL,         (gr_funcptr) nfloat_mat_mul},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIL,  (gr_funcptr) nfloat_mat_nonsingular_solve_tril},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIU,  (gr_funcptr) nfloat_mat_nonsingular_solve_triu},
    {GR_METHOD_MAT_LU,                      (gr_funcptr) nfloat_mat_lu},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT,     (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},

    {0,                         (gr_funcptr) NULL},
};

int
nfloat_ctx_init(gr_ctx_t ctx, slong prec, int flags)
{
    slong nlimbs;

    if (prec <= 0 || prec > NFLOAT_MAX_LIMBS * FLINT_BITS)
        return GR_UNABLE;

    nlimbs = (prec + FLINT_BITS - 1) / FLINT_BITS;

    ctx->which_ring = GR_CTX_NFLOAT;
    ctx->sizeof_elem = sizeof(ulong) * (nlimbs + NFLOAT_HEADER_LIMBS);
    ctx->size_limit = WORD_MAX;

    NFLOAT_CTX_NLIMBS(ctx) = nlimbs;
    NFLOAT_CTX_FLAGS(ctx) = flags;
    NFLOAT_CTX_RND(ctx) = 0;

    ctx->methods = _nfloat_methods;

    if (!_nfloat_methods_initialized)
    {
        gr_method_tab_init(_nfloat_methods, _nfloat_methods_input);
        _nfloat_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

int
nfloat_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_NFLOAT_COMPLEX)
    {
        gr_stream_write(out, "Complex floating-point numbers with prec = ");
        gr_stream_write_si(out, NFLOAT_CTX_PREC(ctx));
        gr_stream_write(out, " (nfloat_complex)");
        return GR_SUCCESS;
    }
    else
    {
        gr_stream_write(out, "Floating-point numbers with prec = ");
        gr_stream_write_si(out, NFLOAT_CTX_PREC(ctx));
        gr_stream_write(out, " (nfloat)");
        return GR_SUCCESS;
    }
}
