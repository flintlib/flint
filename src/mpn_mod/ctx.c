/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "mpn_mod.h"

int _mpn_mod_methods_initialized = 0;

gr_static_method_table _mpn_mod_methods;

#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wcast-function-type"
#endif
gr_method_tab_input _mpn_mod_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) mpn_mod_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) mpn_mod_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) mpn_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) mpn_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) mpn_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_IS_FIELD,(gr_funcptr) mpn_mod_ctx_set_is_field},
    {GR_METHOD_INIT,            (gr_funcptr) mpn_mod_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) mpn_mod_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) mpn_mod_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) mpn_mod_set},
    {GR_METHOD_RANDTEST,        (gr_funcptr) mpn_mod_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) mpn_mod_write},
    {GR_METHOD_ZERO,            (gr_funcptr) mpn_mod_zero},
    {GR_METHOD_ONE,             (gr_funcptr) mpn_mod_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) mpn_mod_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) mpn_mod_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) mpn_mod_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) mpn_mod_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) mpn_mod_equal},
    {GR_METHOD_SET,             (gr_funcptr) mpn_mod_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) mpn_mod_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) mpn_mod_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) mpn_mod_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) mpn_mod_set_other},
    {GR_METHOD_NEG,             (gr_funcptr) mpn_mod_neg},
    {GR_METHOD_ADD,             (gr_funcptr) mpn_mod_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) mpn_mod_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) mpn_mod_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) mpn_mod_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) mpn_mod_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) mpn_mod_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) mpn_mod_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) mpn_mod_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) mpn_mod_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) mpn_mod_mul_si},
    {GR_METHOD_MUL_UI,          (gr_funcptr) mpn_mod_mul_ui},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) mpn_mod_mul_fmpz},
    {GR_METHOD_ADDMUL,          (gr_funcptr) mpn_mod_addmul},
    {GR_METHOD_ADDMUL_SI,       (gr_funcptr) mpn_mod_addmul_si},
    {GR_METHOD_ADDMUL_UI,       (gr_funcptr) mpn_mod_addmul_ui},
    {GR_METHOD_ADDMUL_FMPZ,     (gr_funcptr) mpn_mod_addmul_fmpz},
    {GR_METHOD_SUBMUL,          (gr_funcptr) mpn_mod_submul},
    {GR_METHOD_SUBMUL_SI,       (gr_funcptr) mpn_mod_submul_si},
    {GR_METHOD_SUBMUL_UI,       (gr_funcptr) mpn_mod_submul_ui},
    {GR_METHOD_SUBMUL_FMPZ,     (gr_funcptr) mpn_mod_submul_fmpz},

/*
    {GR_METHOD_MUL_TWO,         (gr_funcptr) mpn_mod_mul_two},
*/
    {GR_METHOD_SQR,             (gr_funcptr) mpn_mod_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) mpn_mod_div},
/*
    {GR_METHOD_DIV_NONUNIQUE,   (gr_funcptr) mpn_mod_div_nonunique},
    {GR_METHOD_DIVIDES,         (gr_funcptr) mpn_mod_divides},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) mpn_mod_is_invertible},
*/


    {GR_METHOD_INV,             (gr_funcptr) mpn_mod_inv},
/*
    {GR_METHOD_POW_SI,          (gr_funcptr) mpn_mod_pow_si},
    {GR_METHOD_POW_UI,          (gr_funcptr) mpn_mod_pow_ui},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) mpn_mod_pow_fmpz},
    {GR_METHOD_SQRT,            (gr_funcptr) mpn_mod_sqrt},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) mpn_mod_is_square},
*/

    {GR_METHOD_VEC_INIT,        (gr_funcptr) _mpn_mod_vec_zero},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) _mpn_mod_vec_clear},
    {GR_METHOD_VEC_SET,         (gr_funcptr) _mpn_mod_vec_set},
    {GR_METHOD_VEC_SWAP,        (gr_funcptr) _mpn_mod_vec_swap},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) _mpn_mod_vec_zero},
    {GR_METHOD_VEC_NEG,         (gr_funcptr) _mpn_mod_vec_neg},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) _mpn_mod_vec_add},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) _mpn_mod_vec_sub},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) _mpn_mod_vec_mul},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) _mpn_mod_vec_mul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR,    (gr_funcptr) _mpn_mod_vec_addmul_scalar},
    {GR_METHOD_SCALAR_MUL_VEC,  (gr_funcptr) _mpn_mod_scalar_mul_vec},

    {GR_METHOD_VEC_DOT,         (gr_funcptr) _mpn_mod_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _mpn_mod_vec_dot_rev},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _mpn_mod_poly_mullow},
    {GR_METHOD_POLY_INV_SERIES, (gr_funcptr) _mpn_mod_poly_inv_series},
    {GR_METHOD_POLY_DIV_SERIES, (gr_funcptr) _mpn_mod_poly_div_series},
    {GR_METHOD_POLY_DIVREM,     (gr_funcptr) _mpn_mod_poly_divrem},
    {GR_METHOD_POLY_DIV,        (gr_funcptr) _mpn_mod_poly_div},
    {GR_METHOD_POLY_GCD,        (gr_funcptr) _mpn_mod_poly_gcd},
    {GR_METHOD_POLY_XGCD,       (gr_funcptr) _mpn_mod_poly_xgcd},
/*
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) mpn_mod_roots_gr_poly},
*/

    {GR_METHOD_MAT_MUL,         (gr_funcptr) mpn_mod_mat_mul},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIL,                 (gr_funcptr) mpn_mod_mat_nonsingular_solve_tril},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIU,                 (gr_funcptr) mpn_mod_mat_nonsingular_solve_triu},
    {GR_METHOD_MAT_LU,          (gr_funcptr) mpn_mod_mat_lu},
    {GR_METHOD_MAT_DET,         (gr_funcptr) mpn_mod_mat_det},
    {0,                         (gr_funcptr) NULL},
};
#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

int
_gr_ctx_init_mpn_mod(gr_ctx_t ctx, mp_srcptr n, mp_size_t nlimbs)
{
    mp_bitcnt_t norm;
    if (nlimbs < MPN_MOD_MIN_LIMBS || nlimbs > MPN_MOD_MAX_LIMBS || n[nlimbs - 1] == 0)
        return GR_UNABLE;

    ctx->which_ring = GR_CTX_MPN_MOD;
    ctx->sizeof_elem = nlimbs * sizeof(mp_limb_t);

    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_mpn_mod_ctx_struct));

    MPN_MOD_CTX_NLIMBS(ctx) = nlimbs;
    flint_mpn_copyi(MPN_MOD_CTX_MODULUS(ctx), n, nlimbs);
    MPN_MOD_CTX_NORM(ctx) = norm = flint_clz(n[nlimbs - 1]);

    if (norm == 0)
        flint_mpn_copyi(MPN_MOD_CTX_MODULUS_NORMED(ctx), n, nlimbs);
    else
        mpn_lshift(MPN_MOD_CTX_MODULUS_NORMED(ctx), n, nlimbs, norm);

    flint_mpn_preinvn(MPN_MOD_CTX_MODULUS_PREINV(ctx), MPN_MOD_CTX_MODULUS_NORMED(ctx), nlimbs);

    MPN_MOD_CTX_IS_PRIME(ctx) = T_UNKNOWN;

    ctx->size_limit = WORD_MAX;

    ctx->methods = _mpn_mod_methods;

    if (!_mpn_mod_methods_initialized)
    {
        gr_method_tab_init(_mpn_mod_methods, _mpn_mod_methods_input);
        _mpn_mod_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

int
gr_ctx_init_mpn_mod(gr_ctx_t ctx, const fmpz_t n)
{
    if (fmpz_sgn(n) <= 0)
        return GR_DOMAIN;

    if (!COEFF_IS_MPZ(*n))
        return GR_UNABLE;

    return _gr_ctx_init_mpn_mod(ctx, COEFF_TO_PTR(*n)->_mp_d, COEFF_TO_PTR(*n)->_mp_size);
}

static const int
randtest_primes[][2] = {
#if FLINT_BITS == 32
    {32, 15},
    {64, -59},
#endif
    {64, 13},
    {128, -159},
    {128, 51},
    {190, 129},
    {192, -237},
    {256, 297},
#if FLINT_BITS == 64
    {1024, -105},
#endif
};

void
gr_ctx_init_mpn_mod_randtest(gr_ctx_t ctx, flint_rand_t state)
{
    fmpz_t n;
    fmpz_init(n);

    if (n_randint(state, 2))
    {
        int i = n_randint(state, sizeof(randtest_primes) / (2 * sizeof(int)));

        fmpz_ui_pow_ui(n, 2, randtest_primes[i][0]);
        fmpz_add_si(n, n, randtest_primes[i][1]);
        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, n));
        GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, n_randint(state, 2) ? T_TRUE : T_UNKNOWN));
    }
    else
    {
        for (;;)
        {
            fmpz_randtest_not_zero(n, state, FLINT_BITS * MPN_MOD_MAX_LIMBS + 10);
            fmpz_abs(n, n);
            if (gr_ctx_init_mpn_mod(ctx, n) == GR_SUCCESS)
                break;
        }
    }

    fmpz_clear(n);
}
