/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_H
#define GR_H

#ifdef GR_INLINES_C
#define GR_INLINE
#else
#define GR_INLINE static inline
#endif

#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

#ifndef CALCIUM_H

typedef enum
{
    T_TRUE,
    T_FALSE,
    T_UNKNOWN
} truth_t;

#endif

GR_INLINE truth_t truth_and(truth_t x, truth_t y)
{
    if (x == T_FALSE || y == T_FALSE)
        return T_FALSE;
    if (x == T_TRUE && y == T_TRUE)
        return T_TRUE;
    return T_UNKNOWN;
}

GR_INLINE truth_t truth_or(truth_t x, truth_t y)
{
    if (x == T_TRUE || y == T_TRUE)
        return T_TRUE;
    if (x == T_FALSE && y == T_FALSE)
        return T_FALSE;
    return T_UNKNOWN;
}

GR_INLINE truth_t truth_not(truth_t x)
{
    if (x == T_TRUE)
        return T_FALSE;
    if (x == T_FALSE)
        return T_TRUE;
    return T_UNKNOWN;
}

GR_INLINE void truth_println(truth_t x)
{
    if (x == T_TRUE) flint_printf("T_TRUE\n");
    if (x == T_FALSE) flint_printf("T_FALSE\n");
    if (x == T_UNKNOWN) flint_printf("T_UNKNOWN\n");
}

typedef int (*gr_funcptr)(void);

/* Copied from Calcium: stream interface allows simple file or string IO. */

typedef struct
{
    FLINT_FILE * fp;
    char * s;
    slong len;
    slong alloc;
}
gr_stream_struct;

typedef gr_stream_struct gr_stream_t[1];

#ifdef FLINT_HAVE_FILE
void gr_stream_init_file(gr_stream_t out, FILE * fp);
#endif

void gr_stream_init_str(gr_stream_t out);
void gr_stream_write(gr_stream_t out, const char * s);
void gr_stream_write_si(gr_stream_t out, slong x);
void gr_stream_write_ui(gr_stream_t out, ulong x);
void gr_stream_write_free(gr_stream_t out, char * s);
void gr_stream_write_fmpz(gr_stream_t out, const fmpz_t x);

#define GR_SUCCESS 0
#define GR_DOMAIN 1
#define GR_UNABLE 2
#define GR_TEST_FAIL 4

#define WARN_UNUSED_RESULT __attribute__((warn_unused_result))

#define GR_MUST_SUCCEED(expr) do { if ((expr) != GR_SUCCESS) { flint_throw(FLINT_ERROR, "GR_MUST_SUCCEED failure: %s", __FILE__); } } while (0)
#define GR_IGNORE(expr) do { int ___unused = (expr); (void) ___unused; } while (0)

typedef void * gr_ptr;
typedef const void * gr_srcptr;
typedef void * gr_ctx_ptr;

#define GR_ENTRY(vec, i, size) ((void *) (((char *) (vec)) + ((i) * (size))))

typedef struct
{
    gr_ptr entries;
    slong alloc;
    slong length;
}
gr_vec_struct;

typedef gr_vec_struct gr_vec_t[1];

GR_INLINE int gr_not_implemented(void) { return GR_UNABLE; }
GR_INLINE int gr_not_in_domain(void) { return GR_DOMAIN; }

typedef enum
{
    GR_METHOD_CTX_WRITE,
    GR_METHOD_CTX_CLEAR,

    /* general ring properties */
    GR_METHOD_CTX_IS_RING,
    GR_METHOD_CTX_IS_COMMUTATIVE_RING,
    GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,
    GR_METHOD_CTX_IS_FIELD,
    GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
    GR_METHOD_CTX_IS_FINITE,
    GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
    GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
    GR_METHOD_CTX_IS_ZERO_RING,

    /* ring properties related to orderings and norms */
    GR_METHOD_CTX_IS_ORDERED_RING,

    /* group properties */
    GR_METHOD_CTX_IS_MULTIPLICATIVE_GROUP,

    /* context properties represented to the representation */
    GR_METHOD_CTX_IS_EXACT,            /* we have no inexact elements */
    GR_METHOD_CTX_IS_CANONICAL,        /* we have no non-canonical representations */

    GR_METHOD_CTX_IS_THREADSAFE,

    GR_METHOD_CTX_HAS_REAL_PREC,
    GR_METHOD_CTX_SET_REAL_PREC,
    GR_METHOD_CTX_GET_REAL_PREC,

    GR_METHOD_CTX_SET_GEN_NAME,
    GR_METHOD_CTX_SET_GEN_NAMES,

    GR_METHOD_INIT,
    GR_METHOD_CLEAR,
    GR_METHOD_SWAP,
    GR_METHOD_SET_SHALLOW,

    GR_METHOD_WRITE,
    GR_METHOD_WRITE_N,

    GR_METHOD_RANDTEST,
    GR_METHOD_RANDTEST_NOT_ZERO,
    GR_METHOD_RANDTEST_SMALL,

    GR_METHOD_ZERO,
    GR_METHOD_ONE,
    GR_METHOD_NEG_ONE,

    GR_METHOD_IS_ZERO,
    GR_METHOD_IS_ONE,
    GR_METHOD_IS_NEG_ONE,

    GR_METHOD_EQUAL,

    GR_METHOD_SET,
    GR_METHOD_SET_UI,
    GR_METHOD_SET_SI,
    GR_METHOD_SET_FMPZ,
    GR_METHOD_SET_FMPQ,
    GR_METHOD_SET_D,
    GR_METHOD_SET_OTHER,
    GR_METHOD_SET_STR,

    GR_METHOD_GET_SI,
    GR_METHOD_GET_UI,
    GR_METHOD_GET_FMPZ,
    GR_METHOD_GET_FMPQ,
    GR_METHOD_GET_D,

    GR_METHOD_GET_FEXPR,
    GR_METHOD_GET_FEXPR_SERIALIZE,
    GR_METHOD_SET_FEXPR,

    GR_METHOD_NEG,

    GR_METHOD_ADD,
    GR_METHOD_ADD_UI,
    GR_METHOD_ADD_SI,
    GR_METHOD_ADD_FMPZ,
    GR_METHOD_ADD_FMPQ,
    GR_METHOD_ADD_OTHER,
    GR_METHOD_OTHER_ADD,

    GR_METHOD_SUB,
    GR_METHOD_SUB_UI,
    GR_METHOD_SUB_SI,
    GR_METHOD_SUB_FMPZ,
    GR_METHOD_SUB_FMPQ,
    GR_METHOD_SUB_OTHER,
    GR_METHOD_OTHER_SUB,

    GR_METHOD_MUL,
    GR_METHOD_MUL_UI,
    GR_METHOD_MUL_SI,
    GR_METHOD_MUL_FMPZ,
    GR_METHOD_MUL_FMPQ,
    GR_METHOD_MUL_OTHER,
    GR_METHOD_OTHER_MUL,

    GR_METHOD_ADDMUL,
    GR_METHOD_ADDMUL_UI,
    GR_METHOD_ADDMUL_SI,
    GR_METHOD_ADDMUL_FMPZ,
    GR_METHOD_ADDMUL_FMPQ,
    GR_METHOD_ADDMUL_OTHER,

    GR_METHOD_SUBMUL,
    GR_METHOD_SUBMUL_UI,
    GR_METHOD_SUBMUL_SI,
    GR_METHOD_SUBMUL_FMPZ,
    GR_METHOD_SUBMUL_FMPQ,
    GR_METHOD_SUBMUL_OTHER,

    GR_METHOD_FMA,
    GR_METHOD_FMS,
    GR_METHOD_FMMA,
    GR_METHOD_FMMS,

    GR_METHOD_MUL_TWO,
    GR_METHOD_SQR,

    GR_METHOD_MUL_2EXP_SI,
    GR_METHOD_MUL_2EXP_FMPZ,
    GR_METHOD_SET_FMPZ_2EXP_FMPZ,
    GR_METHOD_GET_FMPZ_2EXP_FMPZ,

    GR_METHOD_IS_INVERTIBLE,
    GR_METHOD_INV,

    GR_METHOD_DIV,
    GR_METHOD_DIV_UI,
    GR_METHOD_DIV_SI,
    GR_METHOD_DIV_FMPZ,
    GR_METHOD_DIV_FMPQ,
    GR_METHOD_DIV_OTHER,
    GR_METHOD_OTHER_DIV,

    GR_METHOD_DIV_NONUNIQUE,
    GR_METHOD_DIVIDES,

    GR_METHOD_DIVEXACT,
    GR_METHOD_DIVEXACT_UI,
    GR_METHOD_DIVEXACT_SI,
    GR_METHOD_DIVEXACT_FMPZ,
    GR_METHOD_DIVEXACT_FMPQ,
    GR_METHOD_DIVEXACT_OTHER,
    GR_METHOD_OTHER_DIVEXACT,

    GR_METHOD_POW,
    GR_METHOD_POW_UI,
    GR_METHOD_POW_SI,
    GR_METHOD_POW_FMPZ,
    GR_METHOD_POW_FMPQ,
    GR_METHOD_POW_OTHER,
    GR_METHOD_OTHER_POW,

    GR_METHOD_IS_SQUARE,
    GR_METHOD_SQRT,
    GR_METHOD_RSQRT,
    GR_METHOD_HYPOT,

    GR_METHOD_IS_PERFECT_POWER,
    GR_METHOD_FACTOR_PERFECT_POWER,
    GR_METHOD_ROOT_UI,

    GR_METHOD_EUCLIDEAN_DIV,
    GR_METHOD_EUCLIDEAN_REM,
    GR_METHOD_EUCLIDEAN_DIVREM,
    GR_METHOD_GCD,
    GR_METHOD_LCM,

    GR_METHOD_FACTOR,

    GR_METHOD_NUMERATOR,
    GR_METHOD_DENOMINATOR,

    GR_METHOD_FLOOR,
    GR_METHOD_CEIL,
    GR_METHOD_TRUNC,
    GR_METHOD_NINT,

    GR_METHOD_CMP,
    GR_METHOD_CMPABS,
    GR_METHOD_CMP_OTHER,
    GR_METHOD_CMPABS_OTHER,
    GR_METHOD_MIN,
    GR_METHOD_MAX,

    GR_METHOD_ABS,
    GR_METHOD_ABS2,
    GR_METHOD_I,
    GR_METHOD_CONJ,
    GR_METHOD_RE,
    GR_METHOD_IM,
    GR_METHOD_SGN,
    GR_METHOD_CSGN,
    GR_METHOD_ARG,

    GR_METHOD_POS_INF,
    GR_METHOD_NEG_INF,
    GR_METHOD_UINF,
    GR_METHOD_UNDEFINED,
    GR_METHOD_UNKNOWN,

    GR_METHOD_IS_INTEGER,
    GR_METHOD_IS_RATIONAL,
    GR_METHOD_IS_REAL,
    GR_METHOD_IS_POSITIVE_INTEGER,
    GR_METHOD_IS_NONNEGATIVE_INTEGER,
    GR_METHOD_IS_NEGATIVE_INTEGER,
    GR_METHOD_IS_NONPOSITIVE_INTEGER,
    GR_METHOD_IS_POSITIVE_REAL,
    GR_METHOD_IS_NONNEGATIVE_REAL,
    GR_METHOD_IS_NEGATIVE_REAL,
    GR_METHOD_IS_NONPOSITIVE_REAL,

    /* todo: roots of unity */
    GR_METHOD_IS_ROOT_OF_UNITY,
    GR_METHOD_ROOT_OF_UNITY_UI,
    GR_METHOD_ROOT_OF_UNITY_UI_VEC,

    /* Elementary transcendental functions */
    GR_METHOD_PI,
    GR_METHOD_EXP,
    GR_METHOD_EXPM1,
    GR_METHOD_EXP_PI_I,
    GR_METHOD_EXP2,
    GR_METHOD_EXP10,
    GR_METHOD_LOG,
    GR_METHOD_LOG1P,
    GR_METHOD_LOG_PI_I,
    GR_METHOD_LOG2,
    GR_METHOD_LOG10,
    GR_METHOD_SIN,
    GR_METHOD_COS,
    GR_METHOD_SIN_COS,
    GR_METHOD_TAN,
    GR_METHOD_COT,
    GR_METHOD_SEC,
    GR_METHOD_CSC,
    GR_METHOD_SIN_PI,
    GR_METHOD_COS_PI,
    GR_METHOD_SIN_COS_PI,
    GR_METHOD_TAN_PI,
    GR_METHOD_COT_PI,
    GR_METHOD_SEC_PI,
    GR_METHOD_CSC_PI,
    GR_METHOD_SINC,
    GR_METHOD_SINC_PI,
    GR_METHOD_SINH,
    GR_METHOD_COSH,
    GR_METHOD_SINH_COSH,
    GR_METHOD_TANH,
    GR_METHOD_COTH,
    GR_METHOD_SECH,
    GR_METHOD_CSCH,
    GR_METHOD_ASIN,
    GR_METHOD_ACOS,
    GR_METHOD_ATAN,
    GR_METHOD_ATAN2,
    GR_METHOD_ACOT,
    GR_METHOD_ASEC,
    GR_METHOD_ACSC,
    GR_METHOD_ASINH,
    GR_METHOD_ACOSH,
    GR_METHOD_ATANH,
    GR_METHOD_ACOTH,
    GR_METHOD_ASECH,
    GR_METHOD_ACSCH,
    GR_METHOD_ASIN_PI,
    GR_METHOD_ACOS_PI,
    GR_METHOD_ATAN_PI,
    GR_METHOD_ACOT_PI,
    GR_METHOD_ASEC_PI,
    GR_METHOD_ACSC_PI,

    /* Combinatorial and related functions */
    GR_METHOD_FAC,
    GR_METHOD_FAC_UI,
    GR_METHOD_FAC_FMPZ,
    GR_METHOD_FAC_VEC,
    GR_METHOD_RFAC,
    GR_METHOD_RFAC_UI,
    GR_METHOD_RFAC_FMPZ,
    GR_METHOD_RFAC_VEC,
    GR_METHOD_BIN,
    GR_METHOD_BIN_UI,
    GR_METHOD_BIN_UIUI,
    GR_METHOD_BIN_VEC,
    GR_METHOD_BIN_UI_VEC,
    GR_METHOD_RISING_UI,
    GR_METHOD_RISING,
    GR_METHOD_FALLING_UI,
    GR_METHOD_FALLING,
    GR_METHOD_GAMMA,
    GR_METHOD_GAMMA_FMPZ,
    GR_METHOD_GAMMA_FMPQ,
    GR_METHOD_RGAMMA,
    GR_METHOD_LGAMMA,
    GR_METHOD_DIGAMMA,
    GR_METHOD_BETA,
    GR_METHOD_DOUBLEFAC,
    GR_METHOD_DOUBLEFAC_UI,
    GR_METHOD_BARNES_G,
    GR_METHOD_LOG_BARNES_G,
    GR_METHOD_HARMONIC,
    GR_METHOD_HARMONIC_UI,
    GR_METHOD_BERNOULLI_UI,
    GR_METHOD_BERNOULLI_FMPZ,
    GR_METHOD_BERNOULLI_VEC,
    GR_METHOD_FIB_UI,
    GR_METHOD_FIB_FMPZ,
    GR_METHOD_FIB_VEC,
    GR_METHOD_STIRLING_S1U_UIUI,
    GR_METHOD_STIRLING_S1_UIUI,
    GR_METHOD_STIRLING_S2_UIUI,
    GR_METHOD_STIRLING_S1U_UI_VEC,
    GR_METHOD_STIRLING_S1_UI_VEC,
    GR_METHOD_STIRLING_S2_UI_VEC,
    GR_METHOD_EULERNUM_UI,
    GR_METHOD_EULERNUM_FMPZ,
    GR_METHOD_EULERNUM_VEC,
    GR_METHOD_BELLNUM_UI,
    GR_METHOD_BELLNUM_FMPZ,
    GR_METHOD_BELLNUM_VEC,
    GR_METHOD_PARTITIONS_UI,
    GR_METHOD_PARTITIONS_FMPZ,
    GR_METHOD_PARTITIONS_VEC,

    /* Orthogonal polynomials */
    GR_METHOD_CHEBYSHEV_T_FMPZ,
    GR_METHOD_CHEBYSHEV_U_FMPZ,
    GR_METHOD_CHEBYSHEV_T,
    GR_METHOD_CHEBYSHEV_U,
    GR_METHOD_JACOBI_P,
    GR_METHOD_GEGENBAUER_C,
    GR_METHOD_LAGUERRE_L,
    GR_METHOD_HERMITE_H,
    GR_METHOD_LEGENDRE_P,
    GR_METHOD_LEGENDRE_Q,
    GR_METHOD_LEGENDRE_P_ROOT_UI,
    GR_METHOD_SPHERICAL_Y_SI,

    /* Misc real constants */
    GR_METHOD_EULER,
    GR_METHOD_CATALAN,
    GR_METHOD_KHINCHIN,
    GR_METHOD_GLAISHER,

    GR_METHOD_LAMBERTW,
    GR_METHOD_LAMBERTW_FMPZ,

    GR_METHOD_ERF,
    GR_METHOD_ERFC,
    GR_METHOD_ERFCX,
    GR_METHOD_ERFI,
    GR_METHOD_ERFINV,
    GR_METHOD_ERFCINV,
    GR_METHOD_FRESNEL_S,
    GR_METHOD_FRESNEL_C,
    GR_METHOD_FRESNEL,
    GR_METHOD_GAMMA_UPPER,
    GR_METHOD_GAMMA_LOWER,
    GR_METHOD_BETA_LOWER,
    GR_METHOD_EXP_INTEGRAL,
    GR_METHOD_EXP_INTEGRAL_EI,
    GR_METHOD_SIN_INTEGRAL,
    GR_METHOD_COS_INTEGRAL,
    GR_METHOD_SINH_INTEGRAL,
    GR_METHOD_COSH_INTEGRAL,
    GR_METHOD_LOG_INTEGRAL,
    GR_METHOD_DILOG,

    GR_METHOD_BESSEL_J,
    GR_METHOD_BESSEL_Y,
    GR_METHOD_BESSEL_J_Y,
    GR_METHOD_BESSEL_I,
    GR_METHOD_BESSEL_I_SCALED,
    GR_METHOD_BESSEL_K,
    GR_METHOD_BESSEL_K_SCALED,

    GR_METHOD_AIRY,
    GR_METHOD_AIRY_AI,
    GR_METHOD_AIRY_BI,
    GR_METHOD_AIRY_AI_PRIME,
    GR_METHOD_AIRY_BI_PRIME,
    GR_METHOD_AIRY_AI_ZERO,
    GR_METHOD_AIRY_BI_ZERO,
    GR_METHOD_AIRY_AI_PRIME_ZERO,
    GR_METHOD_AIRY_BI_PRIME_ZERO,

    GR_METHOD_COULOMB,
    GR_METHOD_COULOMB_F,
    GR_METHOD_COULOMB_G,
    GR_METHOD_COULOMB_HPOS,
    GR_METHOD_COULOMB_HNEG,

    GR_METHOD_ZETA,
    GR_METHOD_ZETA_UI,
    GR_METHOD_HURWITZ_ZETA,
    GR_METHOD_LERCH_PHI,
    GR_METHOD_STIELTJES,
    GR_METHOD_DIRICHLET_ETA,
    GR_METHOD_DIRICHLET_BETA,
    GR_METHOD_RIEMANN_XI,
    GR_METHOD_ZETA_ZERO,
    GR_METHOD_ZETA_ZERO_VEC,
    GR_METHOD_ZETA_NZEROS,

    GR_METHOD_DIRICHLET_CHI_UI,
    GR_METHOD_DIRICHLET_CHI_FMPZ,
    GR_METHOD_DIRICHLET_L,
    GR_METHOD_DIRICHLET_HARDY_THETA,
    GR_METHOD_DIRICHLET_HARDY_Z,

    GR_METHOD_BERNPOLY_UI,
    GR_METHOD_EULERPOLY_UI,
    GR_METHOD_POLYLOG,
    GR_METHOD_POLYGAMMA,

    GR_METHOD_HYPGEOM_0F1,
    GR_METHOD_HYPGEOM_1F1,
    GR_METHOD_HYPGEOM_2F0,
    GR_METHOD_HYPGEOM_2F1,
    GR_METHOD_HYPGEOM_U,
    GR_METHOD_HYPGEOM_PFQ,

    GR_METHOD_JACOBI_THETA,
    GR_METHOD_JACOBI_THETA_1,
    GR_METHOD_JACOBI_THETA_2,
    GR_METHOD_JACOBI_THETA_3,
    GR_METHOD_JACOBI_THETA_4,
    GR_METHOD_JACOBI_THETA_Q,
    GR_METHOD_JACOBI_THETA_Q_1,
    GR_METHOD_JACOBI_THETA_Q_2,
    GR_METHOD_JACOBI_THETA_Q_3,
    GR_METHOD_JACOBI_THETA_Q_4,

    GR_METHOD_MODULAR_J,
    GR_METHOD_MODULAR_LAMBDA,
    GR_METHOD_MODULAR_DELTA,
    GR_METHOD_HILBERT_CLASS_POLY,

    GR_METHOD_DEDEKIND_ETA,
    GR_METHOD_DEDEKIND_ETA_Q,

    GR_METHOD_EISENSTEIN_E,
    GR_METHOD_EISENSTEIN_G,
    GR_METHOD_EISENSTEIN_G_VEC,

    GR_METHOD_AGM,
    GR_METHOD_AGM1,

    GR_METHOD_ELLIPTIC_K,
    GR_METHOD_ELLIPTIC_E,
    GR_METHOD_ELLIPTIC_PI,
    GR_METHOD_ELLIPTIC_F,
    GR_METHOD_ELLIPTIC_E_INC,
    GR_METHOD_ELLIPTIC_PI_INC,

    GR_METHOD_CARLSON_RF,
    GR_METHOD_CARLSON_RC,
    GR_METHOD_CARLSON_RJ,
    GR_METHOD_CARLSON_RG,
    GR_METHOD_CARLSON_RD,

    GR_METHOD_ELLIPTIC_ROOTS,
    GR_METHOD_ELLIPTIC_INVARIANTS,

    GR_METHOD_WEIERSTRASS_P,
    GR_METHOD_WEIERSTRASS_P_PRIME,
    GR_METHOD_WEIERSTRASS_P_INV,
    GR_METHOD_WEIERSTRASS_ZETA,
    GR_METHOD_WEIERSTRASS_SIGMA,

    GR_METHOD_GEN,
    GR_METHOD_GENS,

    /* Finite field methods */
    GR_METHOD_CTX_FQ_PRIME,
    GR_METHOD_CTX_FQ_DEGREE,
    GR_METHOD_CTX_FQ_ORDER,
    GR_METHOD_FQ_FROBENIUS,
    GR_METHOD_FQ_MULTIPLICATIVE_ORDER,
    GR_METHOD_FQ_NORM,
    GR_METHOD_FQ_TRACE,
    GR_METHOD_FQ_IS_PRIMITIVE,
    GR_METHOD_FQ_PTH_ROOT,

    /* Vector methods */
    GR_METHOD_VEC_INIT,
    GR_METHOD_VEC_CLEAR,
    GR_METHOD_VEC_SWAP,
    GR_METHOD_VEC_SET,
    GR_METHOD_VEC_ZERO,
    GR_METHOD_VEC_EQUAL,
    GR_METHOD_VEC_IS_ZERO,
    GR_METHOD_VEC_NEG,

    GR_METHOD_VEC_NORMALISE,
    GR_METHOD_VEC_NORMALISE_WEAK,

    GR_METHOD_VEC_ADD, GR_METHOD_VEC_SUB, GR_METHOD_VEC_MUL, GR_METHOD_VEC_DIV, GR_METHOD_VEC_DIVEXACT, GR_METHOD_VEC_POW,
    GR_METHOD_VEC_ADD_SCALAR, GR_METHOD_VEC_SUB_SCALAR, GR_METHOD_VEC_MUL_SCALAR, GR_METHOD_VEC_DIV_SCALAR, GR_METHOD_VEC_DIVEXACT_SCALAR, GR_METHOD_VEC_POW_SCALAR,
    GR_METHOD_SCALAR_ADD_VEC, GR_METHOD_SCALAR_SUB_VEC, GR_METHOD_SCALAR_MUL_VEC, GR_METHOD_SCALAR_DIV_VEC, GR_METHOD_SCALAR_DIVEXACT_VEC, GR_METHOD_SCALAR_POW_VEC,
    GR_METHOD_VEC_ADD_OTHER, GR_METHOD_VEC_SUB_OTHER, GR_METHOD_VEC_MUL_OTHER, GR_METHOD_VEC_DIV_OTHER, GR_METHOD_VEC_DIVEXACT_OTHER, GR_METHOD_VEC_POW_OTHER,
    GR_METHOD_OTHER_ADD_VEC, GR_METHOD_OTHER_SUB_VEC, GR_METHOD_OTHER_MUL_VEC, GR_METHOD_OTHER_DIV_VEC, GR_METHOD_OTHER_DIVEXACT_VEC, GR_METHOD_OTHER_POW_VEC,
    GR_METHOD_VEC_ADD_SCALAR_OTHER, GR_METHOD_VEC_SUB_SCALAR_OTHER, GR_METHOD_VEC_MUL_SCALAR_OTHER, GR_METHOD_VEC_DIV_SCALAR_OTHER, GR_METHOD_VEC_DIVEXACT_SCALAR_OTHER, GR_METHOD_VEC_POW_SCALAR_OTHER,
    GR_METHOD_SCALAR_OTHER_ADD_VEC, GR_METHOD_SCALAR_OTHER_SUB_VEC, GR_METHOD_SCALAR_OTHER_MUL_VEC, GR_METHOD_SCALAR_OTHER_DIV_VEC, GR_METHOD_SCALAR_OTHER_DIVEXACT_VEC, GR_METHOD_SCALAR_OTHER_POW_VEC,

    GR_METHOD_VEC_ADD_SCALAR_SI, GR_METHOD_VEC_ADD_SCALAR_UI, GR_METHOD_VEC_ADD_SCALAR_FMPZ, GR_METHOD_VEC_ADD_SCALAR_FMPQ,
    GR_METHOD_VEC_SUB_SCALAR_SI, GR_METHOD_VEC_SUB_SCALAR_UI, GR_METHOD_VEC_SUB_SCALAR_FMPZ, GR_METHOD_VEC_SUB_SCALAR_FMPQ,
    GR_METHOD_VEC_MUL_SCALAR_SI, GR_METHOD_VEC_MUL_SCALAR_UI, GR_METHOD_VEC_MUL_SCALAR_FMPZ, GR_METHOD_VEC_MUL_SCALAR_FMPQ,
    GR_METHOD_VEC_DIV_SCALAR_SI, GR_METHOD_VEC_DIV_SCALAR_UI, GR_METHOD_VEC_DIV_SCALAR_FMPZ, GR_METHOD_VEC_DIV_SCALAR_FMPQ,
    GR_METHOD_VEC_DIVEXACT_SCALAR_SI, GR_METHOD_VEC_DIVEXACT_SCALAR_UI, GR_METHOD_VEC_DIVEXACT_SCALAR_FMPZ, GR_METHOD_VEC_DIVEXACT_SCALAR_FMPQ,
    GR_METHOD_VEC_POW_SCALAR_SI, GR_METHOD_VEC_POW_SCALAR_UI, GR_METHOD_VEC_POW_SCALAR_FMPZ, GR_METHOD_VEC_POW_SCALAR_FMPQ,

    GR_METHOD_VEC_MUL_SCALAR_2EXP_SI,
    GR_METHOD_VEC_MUL_SCALAR_2EXP_FMPZ,

    GR_METHOD_VEC_ADDMUL_SCALAR,
    GR_METHOD_VEC_SUBMUL_SCALAR,
    GR_METHOD_VEC_ADDMUL_SCALAR_SI,
    GR_METHOD_VEC_SUBMUL_SCALAR_SI,

    GR_METHOD_VEC_SUM,
    GR_METHOD_VEC_PRODUCT,

    GR_METHOD_VEC_DOT,
    GR_METHOD_VEC_DOT_REV,
    GR_METHOD_VEC_DOT_UI,
    GR_METHOD_VEC_DOT_SI,
    GR_METHOD_VEC_DOT_FMPZ,

    GR_METHOD_VEC_SET_POWERS,
    GR_METHOD_VEC_RECIPROCALS,

    /* Polynomial methods (todo: rename -> GR_POLY) */
    GR_METHOD_POLY_MULLOW,
    GR_METHOD_POLY_DIV,
    GR_METHOD_POLY_DIVREM,
    GR_METHOD_POLY_DIVEXACT,
    GR_METHOD_POLY_TAYLOR_SHIFT,
    GR_METHOD_POLY_INV_SERIES,
    GR_METHOD_POLY_INV_SERIES_BASECASE,
    GR_METHOD_POLY_DIV_SERIES,
    GR_METHOD_POLY_DIV_SERIES_BASECASE,
    GR_METHOD_POLY_RSQRT_SERIES,
    GR_METHOD_POLY_SQRT_SERIES,
    GR_METHOD_POLY_EXP_SERIES,
    GR_METHOD_POLY_FACTOR,
    GR_METHOD_POLY_ROOTS,
    GR_METHOD_POLY_ROOTS_OTHER,

    /* Matrix methods (todo: rename -> GR_MAT) */
    GR_METHOD_MAT_MUL,
    GR_METHOD_MAT_DET,
    GR_METHOD_MAT_EXP,
    GR_METHOD_MAT_LOG,
    GR_METHOD_MAT_FIND_NONZERO_PIVOT,
    GR_METHOD_MAT_DIAGONALIZATION,

    GR_METHOD_TAB_SIZE
}
gr_method;

typedef gr_funcptr gr_static_method_table[GR_METHOD_TAB_SIZE];

typedef struct
{
    gr_method index;
    gr_funcptr function;
}
gr_method_tab_input;

void gr_method_tab_init(gr_funcptr * methods, gr_method_tab_input * tab);

/* Identify specific rings/fields. */

typedef enum
{
    GR_CTX_FMPZ, GR_CTX_FMPQ, GR_CTX_FMPZI,
    GR_CTX_FMPZ_MOD, GR_CTX_NMOD, GR_CTX_NMOD8, GR_CTX_NMOD32,
    GR_CTX_FQ, GR_CTX_FQ_NMOD, GR_CTX_FQ_ZECH,
    GR_CTX_NF,
    GR_CTX_REAL_ALGEBRAIC_QQBAR, GR_CTX_COMPLEX_ALGEBRAIC_QQBAR,
    GR_CTX_REAL_ALGEBRAIC_CA, GR_CTX_COMPLEX_ALGEBRAIC_CA,
    GR_CTX_RR_CA, GR_CTX_CC_CA,
    GR_CTX_COMPLEX_EXTENDED_CA,
    GR_CTX_RR_ARB, GR_CTX_CC_ACB,
    GR_CTX_REAL_FLOAT_ARF, GR_CTX_COMPLEX_FLOAT_ACF,
    GR_CTX_FMPZ_POLY, GR_CTX_FMPQ_POLY, GR_CTX_GR_POLY,
    GR_CTX_FMPZ_MPOLY, GR_CTX_GR_MPOLY,
    GR_CTX_FMPZ_MPOLY_Q,
    GR_CTX_GR_SERIES, GR_CTX_GR_SERIES_MOD,
    GR_CTX_GR_MAT,
    GR_CTX_GR_VEC,
    GR_CTX_PSL2Z, GR_CTX_DIRICHLET_GROUP, GR_CTX_PERM,
    GR_CTX_FEXPR,
    GR_CTX_UNKNOWN_DOMAIN,
    GR_CTX_WHICH_STRUCTURE_TAB_SIZE
}
gr_which_structure;

/* large enough to hold any context data we want to store inline */
#define GR_CTX_STRUCT_DATA_BYTES (6 * sizeof(ulong))

typedef struct
{
    char data[GR_CTX_STRUCT_DATA_BYTES];
    ulong which_ring;
    slong sizeof_elem;
    gr_funcptr * methods;
    ulong size_limit;
}
gr_ctx_struct;

typedef gr_ctx_struct gr_ctx_t[1];

#define GR_CTX_DATA_AS_PTR(ctx) (*(void **) (&(ctx)->data))

GR_INLINE void * gr_ctx_data_ptr(gr_ctx_t ctx) { return (void *) ctx->data; }
GR_INLINE void * gr_ctx_data_as_ptr(gr_ctx_t ctx) { return (void *) GR_CTX_DATA_AS_PTR(ctx); }

GR_INLINE slong gr_ctx_sizeof_ctx(void)
{
    return sizeof(gr_ctx_struct);
}

GR_INLINE slong gr_ctx_sizeof_elem(gr_ctx_t ctx)
{
    return ctx->sizeof_elem;
}

/* Typedefs for method function pointers. */

typedef void ((*gr_method_init_clear_op)(gr_ptr, gr_ctx_ptr));
typedef void ((*gr_method_swap_op)(gr_ptr, gr_ptr, gr_ctx_ptr));
typedef int ((*gr_method_ctx)(gr_ctx_ptr));
typedef truth_t ((*gr_method_ctx_predicate)(gr_ctx_ptr));
typedef int ((*gr_method_ctx_set_si)(gr_ctx_ptr, slong));
typedef int ((*gr_method_ctx_get_si)(slong *, gr_ctx_ptr));
typedef int ((*gr_method_ctx_stream)(gr_stream_t, gr_ctx_ptr));
typedef int ((*gr_method_ctx_set_str)(gr_ctx_ptr, const char *));
typedef int ((*gr_method_ctx_set_strs)(gr_ctx_ptr, const char **));
typedef int ((*gr_method_stream_in)(gr_stream_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_stream_in_si)(gr_stream_t, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_randtest)(gr_ptr, flint_rand_t state, gr_ctx_ptr));
typedef int ((*gr_method_constant_op)(gr_ptr, gr_ctx_ptr));
typedef int ((*gr_method_constant_op_get_si)(slong *, gr_ctx_ptr));
typedef int ((*gr_method_constant_op_get_fmpz)(fmpz_t, gr_ctx_ptr));
typedef void ((*gr_method_void_unary_op)(gr_ptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op)(gr_ptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_si)(gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_ui)(gr_ptr, ulong, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_fmpz)(gr_ptr, const fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_fmpq)(gr_ptr, const fmpq_t, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_d)(gr_ptr, double, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_other)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_str)(gr_ptr, const char *, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_get_ui)(ulong *, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_get_si)(slong *, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_get_fmpz)(fmpz_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_get_fmpq)(fmpq_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_get_d)(double *, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_get_fmpz_fmpz)(fmpz_t, fmpz_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_with_flag)(gr_ptr, gr_srcptr, int, gr_ctx_ptr));
typedef int ((*gr_method_binary_unary_op)(gr_ptr, gr_ptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_unary_op_with_flag)(gr_ptr, gr_ptr, gr_srcptr, int, gr_ctx_ptr));
typedef int ((*gr_method_quaternary_unary_op)(gr_ptr, gr_ptr, gr_ptr, gr_ptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_si)(gr_ptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_ui)(gr_ptr, gr_srcptr, ulong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_fmpz)(gr_ptr, gr_srcptr, const fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_fmpz_fmpz)(gr_ptr, const fmpz_t, const fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_fmpz_si)(gr_ptr, const fmpz_t, slong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_fmpq)(gr_ptr, gr_srcptr, const fmpq_t, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_other)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr));
typedef int ((*gr_method_other_binary_op)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_si_binary_op)(gr_ptr, slong, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_ui_binary_op)(gr_ptr, ulong, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_fmpz_binary_op)(gr_ptr, const fmpz_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_fmpq_binary_op)(gr_ptr, const fmpq_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_ui_ui)(gr_ptr, ulong, ulong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_ui_si)(gr_ptr, ulong, slong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_get_int)(int *, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_other_get_int)(int *, gr_srcptr, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_binary_op)(gr_ptr, gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_with_flag)(gr_ptr, gr_srcptr, gr_srcptr, int, gr_ctx_ptr));
typedef int ((*gr_method_binary_binary_op_ui_ui)(gr_ptr, gr_ptr, ulong, ulong, gr_ctx_ptr));
typedef int ((*gr_method_ternary_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_ternary_op_with_flag)(gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, int, gr_ctx_ptr));
typedef int ((*gr_method_ternary_unary_op)(gr_ptr, gr_ptr, gr_ptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_quaternary_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_quaternary_op_with_flag)(gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_srcptr, int, gr_ctx_ptr));
typedef int ((*gr_method_quaternary_binary_op)(gr_ptr, gr_ptr, gr_ptr, gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_quaternary_ternary_op)(gr_ptr, gr_ptr, gr_ptr, gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_si_si_quaternary_op)(gr_ptr, slong, slong, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef truth_t ((*gr_method_unary_predicate)(gr_srcptr, gr_ctx_ptr));
typedef truth_t ((*gr_method_binary_predicate)(gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef void ((*gr_method_vec_init_clear_op)(gr_ptr, slong, gr_ctx_ptr));
typedef void ((*gr_method_vec_swap_op)(gr_ptr, gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_constant_op)(gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_op)(gr_ptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_vec_op)(gr_ptr, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_scalar_vec_op)(gr_ptr, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_op_other)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_other_op_vec)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_op_scalar_other)(gr_ptr, gr_srcptr, slong, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr));
typedef int ((*gr_method_scalar_other_op_vec)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op_si)(gr_ptr, gr_srcptr, slong, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op_ui)(gr_ptr, gr_srcptr, slong, ulong, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op_fmpz)(gr_ptr, gr_srcptr, slong, const fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op_fmpq)(gr_ptr, gr_srcptr, slong, const fmpq_t, gr_ctx_ptr));
typedef truth_t ((*gr_method_vec_predicate)(gr_srcptr, slong, gr_ctx_ptr));
typedef truth_t ((*gr_method_vec_vec_predicate)(gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_factor_op)(gr_ptr, gr_vec_t, gr_vec_t, gr_srcptr, int, gr_ctx_ptr));
typedef int ((*gr_method_poly_unary_trunc_op)(gr_ptr, gr_srcptr, slong, slong, gr_ctx_ptr));
typedef int ((*gr_method_poly_binary_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_poly_binary_binary_op)(gr_ptr, gr_ptr, gr_srcptr, slong, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_poly_binary_trunc_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, slong, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_ctx_op)(gr_vec_t, gr_ctx_ptr));

#ifdef FEXPR_H
typedef int ((*gr_method_get_fexpr_op)(fexpr_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_set_fexpr_op)(gr_ptr, fexpr_vec_t, gr_vec_t, const fexpr_t, gr_ctx_ptr));
#endif

/* Macros to retrieve methods (with correct call signature) from context object. */
#define GR_CTX_OP(ctx, NAME) (((gr_method_ctx *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CTX_STREAM(ctx, NAME) (((gr_method_ctx_stream *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CTX_PREDICATE(ctx, NAME) (((gr_method_ctx_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CTX_SET_SI(ctx, NAME) (((gr_method_ctx_set_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CTX_GET_SI(ctx, NAME) (((gr_method_ctx_get_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CTX_SET_STR(ctx, NAME) (((gr_method_ctx_set_str *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CTX_SET_STRS(ctx, NAME) (((gr_method_ctx_set_strs *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_STREAM_IN(ctx, NAME) (((gr_method_stream_in *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_STREAM_IN_SI(ctx, NAME) (((gr_method_stream_in_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_RANDTEST(ctx, NAME) (((gr_method_randtest *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_INIT_CLEAR_OP(ctx, NAME) (((gr_method_init_clear_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SWAP_OP(ctx, NAME) (((gr_method_swap_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CONSTANT_OP(ctx, NAME) (((gr_method_constant_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CONSTANT_OP_GET_SI(ctx, NAME) (((gr_method_constant_op_get_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CONSTANT_OP_GET_FMPZ(ctx, NAME) (((gr_method_constant_op_get_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VOID_UNARY_OP(ctx, NAME) (((gr_method_void_unary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP(ctx, NAME) (((gr_method_unary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_SI(ctx, NAME) (((gr_method_unary_op_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_UI(ctx, NAME) (((gr_method_unary_op_ui *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_FMPZ(ctx, NAME) (((gr_method_unary_op_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_FMPQ(ctx, NAME) (((gr_method_unary_op_fmpq *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_D(ctx, NAME) (((gr_method_unary_op_d *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_OTHER(ctx, NAME) (((gr_method_unary_op_other *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_STR(ctx, NAME) (((gr_method_unary_op_str *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_GET_SI(ctx, NAME) (((gr_method_unary_op_get_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_GET_UI(ctx, NAME) (((gr_method_unary_op_get_ui *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_GET_FMPZ(ctx, NAME) (((gr_method_unary_op_get_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_GET_FMPQ(ctx, NAME) (((gr_method_unary_op_get_fmpq *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_GET_D(ctx, NAME) (((gr_method_unary_op_get_d *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_GET_FMPZ_FMPZ(ctx, NAME) (((gr_method_unary_op_get_fmpz_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_WITH_FLAG(ctx, NAME) (((gr_method_unary_op_with_flag *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_UNARY_OP(ctx, NAME) (((gr_method_binary_unary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_UNARY_OP_WITH_FLAG(ctx, NAME) (((gr_method_binary_unary_op_with_flag *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_TERNARY_UNARY_OP(ctx, NAME) (((gr_method_ternary_unary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_QUATERNARY_UNARY_OP(ctx, NAME) (((gr_method_quaternary_unary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP(ctx, NAME) (((gr_method_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_SI(ctx, NAME) (((gr_method_binary_op_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_UI(ctx, NAME) (((gr_method_binary_op_ui *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_FMPZ(ctx, NAME) (((gr_method_binary_op_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_FMPQ(ctx, NAME) (((gr_method_binary_op_fmpq *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_OTHER(ctx, NAME) (((gr_method_binary_op_other *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_FMPZ_FMPZ(ctx, NAME) (((gr_method_binary_op_fmpz_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_FMPZ_SI(ctx, NAME) (((gr_method_binary_op_fmpz_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_UI_UI(ctx, NAME) (((gr_method_binary_op_ui_ui *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_UI_SI(ctx, NAME) (((gr_method_binary_op_ui_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_OTHER_BINARY_OP(ctx, NAME) (((gr_method_other_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SI_BINARY_OP(ctx, NAME) (((gr_method_si_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UI_BINARY_OP(ctx, NAME) (((gr_method_ui_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_FMPZ_BINARY_OP(ctx, NAME) (((gr_method_fmpz_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_FMPQ_BINARY_OP(ctx, NAME) (((gr_method_fmpq_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_GET_INT(ctx, NAME) (((gr_method_binary_op_get_int *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_OTHER_GET_INT(ctx, NAME) (((gr_method_binary_op_other_get_int *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_BINARY_OP(ctx, NAME) (((gr_method_binary_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_PREDICATE(ctx, NAME) (((gr_method_unary_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_PREDICATE(ctx, NAME) (((gr_method_binary_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_WITH_FLAG(ctx, NAME) (((gr_method_binary_op_with_flag *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_BINARY_OP_UI_UI(ctx, NAME) (((gr_method_binary_binary_op_ui_ui *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_TERNARY_OP(ctx, NAME) (((gr_method_ternary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_TERNARY_OP_WITH_FLAG(ctx, NAME) (((gr_method_ternary_op_with_flag *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_QUATERNARY_OP(ctx, NAME) (((gr_method_quaternary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_QUATERNARY_OP_WITH_FLAG(ctx, NAME) (((gr_method_quaternary_op_with_flag *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_QUATERNARY_BINARY_OP(ctx, NAME) (((gr_method_quaternary_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_QUATERNARY_TERNARY_OP(ctx, NAME) (((gr_method_quaternary_ternary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SI_SI_QUATERNARY_OP(ctx, NAME) (((gr_method_si_si_quaternary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_INIT_CLEAR_OP(ctx, NAME) (((gr_method_vec_init_clear_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SWAP_OP(ctx, NAME) (((gr_method_vec_swap_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_CONSTANT_OP(ctx, NAME) (((gr_method_vec_constant_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_OP(ctx, NAME) (((gr_method_vec_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_VEC_OP(ctx, NAME) (((gr_method_vec_vec_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP(ctx, NAME) (((gr_method_vec_scalar_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SCALAR_VEC_OP(ctx, NAME) (((gr_method_scalar_vec_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_OP_OTHER(ctx, NAME) (((gr_method_vec_op_other *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_OTHER_OP_VEC(ctx, NAME) (((gr_method_other_op_vec *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_OP_SCALAR_OTHER(ctx, NAME) (((gr_method_vec_op_scalar_other *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SCALAR_OTHER_OP_VEC(ctx, NAME) (((gr_method_scalar_other_op_vec *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP_SI(ctx, NAME) (((gr_method_vec_scalar_op_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP_UI(ctx, NAME) (((gr_method_vec_scalar_op_ui *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP_FMPZ(ctx, NAME) (((gr_method_vec_scalar_op_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP_FMPQ(ctx, NAME) (((gr_method_vec_scalar_op_fmpq *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_PREDICATE(ctx, NAME) (((gr_method_vec_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_VEC_PREDICATE(ctx, NAME) (((gr_method_vec_vec_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_FACTOR_OP(ctx, NAME) (((gr_method_factor_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_POLY_BINARY_OP(ctx, NAME) (((gr_method_poly_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_POLY_UNARY_TRUNC_OP(ctx, NAME) (((gr_method_poly_unary_trunc_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_POLY_BINARY_BINARY_OP(ctx, NAME) (((gr_method_poly_binary_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_POLY_BINARY_TRUNC_OP(ctx, NAME) (((gr_method_poly_binary_trunc_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_CTX_OP(ctx, NAME) (((gr_method_vec_ctx_op *) ctx->methods)[GR_METHOD_ ## NAME])
#ifdef FEXPR_H
#define GR_GET_FEXPR_OP(ctx, NAME) (((gr_method_get_fexpr_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SET_FEXPR_OP(ctx, NAME) (((gr_method_set_fexpr_op *) ctx->methods)[GR_METHOD_ ## NAME])
#endif

/* Wrappers to call methods. */

GR_INLINE int gr_ctx_clear(gr_ctx_t ctx) { return GR_CTX_OP(ctx, CTX_CLEAR)(ctx); }
GR_INLINE int gr_ctx_write(gr_stream_t out, gr_ctx_t ctx) { return GR_CTX_STREAM(ctx, CTX_WRITE)(out, ctx); }

GR_INLINE truth_t gr_ctx_is_ring(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_RING)(ctx); }
GR_INLINE truth_t gr_ctx_is_commutative_ring(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_COMMUTATIVE_RING)(ctx); }
GR_INLINE truth_t gr_ctx_is_integral_domain(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_INTEGRAL_DOMAIN)(ctx); }
GR_INLINE truth_t gr_ctx_is_field(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_FIELD)(ctx); }
GR_INLINE truth_t gr_ctx_is_zero_ring(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_ZERO_RING)(ctx); }

GR_INLINE truth_t gr_ctx_is_unique_factorization_domain(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_UNIQUE_FACTORIZATION_DOMAIN)(ctx); }
GR_INLINE truth_t gr_ctx_is_finite(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_FINITE)(ctx); }
GR_INLINE truth_t gr_ctx_is_finite_characteristic(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_FINITE_CHARACTERISTIC)(ctx); }
GR_INLINE truth_t gr_ctx_is_algebraically_closed(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_ALGEBRAICALLY_CLOSED)(ctx); }
GR_INLINE truth_t gr_ctx_is_ordered_ring(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_ORDERED_RING)(ctx); }

GR_INLINE truth_t gr_ctx_is_multiplicative_group(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_MULTIPLICATIVE_GROUP)(ctx); }

GR_INLINE truth_t gr_ctx_is_exact(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_EXACT)(ctx); }
GR_INLINE truth_t gr_ctx_is_canonical(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_CANONICAL)(ctx); }

GR_INLINE truth_t gr_ctx_is_threadsafe(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_THREADSAFE)(ctx); }

GR_INLINE truth_t gr_ctx_has_real_prec(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_HAS_REAL_PREC)(ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_ctx_set_real_prec(gr_ctx_t ctx, slong prec) { return GR_CTX_SET_SI(ctx, CTX_SET_REAL_PREC)(ctx, prec); }
GR_INLINE WARN_UNUSED_RESULT int gr_ctx_get_real_prec(slong * prec, gr_ctx_t ctx) { return GR_CTX_GET_SI(ctx, CTX_GET_REAL_PREC)(prec, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_ctx_set_gen_name(gr_ctx_t ctx, const char * s) { return GR_CTX_SET_STR(ctx, CTX_SET_GEN_NAME)(ctx, s); }
GR_INLINE WARN_UNUSED_RESULT int gr_ctx_set_gen_names(gr_ctx_t ctx, const char ** s) { return GR_CTX_SET_STRS(ctx, CTX_SET_GEN_NAMES)(ctx, s); }

GR_INLINE slong _gr_ctx_get_real_prec(gr_ctx_t ctx)
{
    slong res = 0;
    GR_IGNORE(gr_ctx_get_real_prec(&res, ctx));
    return res;
}

GR_INLINE void gr_init(gr_ptr res, gr_ctx_t ctx) { GR_INIT_CLEAR_OP(ctx, INIT)(res, ctx); }
GR_INLINE void gr_clear(gr_ptr res, gr_ctx_t ctx) { GR_INIT_CLEAR_OP(ctx, CLEAR)(res, ctx); }
GR_INLINE void gr_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx) { GR_SWAP_OP(ctx, SWAP)(x, y, ctx); }
GR_INLINE void gr_set_shallow(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { GR_VOID_UNARY_OP(ctx, SET_SHALLOW)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_randtest(gr_ptr x, flint_rand_t state, gr_ctx_t ctx) { return GR_RANDTEST(ctx, RANDTEST)(x, state, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_randtest_not_zero(gr_ptr x, flint_rand_t state, gr_ctx_t ctx) { return GR_RANDTEST(ctx, RANDTEST_NOT_ZERO)(x, state, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_randtest_small(gr_ptr x, flint_rand_t state, gr_ctx_t ctx) { return GR_RANDTEST(ctx, RANDTEST_SMALL)(x, state, ctx); }
GR_INLINE /* todo: warn? */ int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx) { return GR_STREAM_IN(ctx, WRITE)(out, x, ctx); }
GR_INLINE /* todo: warn? */ int gr_write_n(gr_stream_t out, gr_srcptr x, slong n, gr_ctx_t ctx) { return GR_STREAM_IN_SI(ctx, WRITE_N)(out, x, n, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_zero(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, ZERO)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_one(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, ONE)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_neg_one(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, NEG_ONE)(res, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SET)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_si(gr_ptr res, slong x, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, SET_SI)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, SET_UI)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, SET_FMPZ)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_fmpq(gr_ptr res, const fmpq_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPQ(ctx, SET_FMPQ)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_d(gr_ptr res, double x, gr_ctx_t ctx) { return GR_UNARY_OP_D(ctx, SET_D)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx) { return GR_UNARY_OP_OTHER(ctx, SET_OTHER)(res, x, x_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_str(gr_ptr res, const char * x, gr_ctx_t ctx) { return GR_UNARY_OP_STR(ctx, SET_STR)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_get_si(slong * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_SI(ctx, GET_SI)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_get_ui(ulong * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_UI(ctx, GET_UI)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_get_fmpz(fmpz_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPZ(ctx, GET_FMPZ)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_get_fmpq(fmpq_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPQ(ctx, GET_FMPQ)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_get_d(double * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_D(ctx, GET_D)(res, x, ctx); }

#define GR_FEXPR_SERIALIZE (UWORD(1) << 31)
#define GR_FEXPR_COMPACT (UWORD(1) << 30)

#ifdef FEXPR_H
GR_INLINE WARN_UNUSED_RESULT int gr_get_fexpr(fexpr_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_GET_FEXPR_OP(ctx, GET_FEXPR)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_get_fexpr_serialize(fexpr_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_GET_FEXPR_OP(ctx, GET_FEXPR_SERIALIZE)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_fexpr(gr_ptr res, fexpr_vec_t inputs, gr_vec_t outputs, const fexpr_t x, gr_ctx_t ctx) { return GR_SET_FEXPR_OP(ctx, SET_FEXPR)(res, inputs, outputs, x, ctx); }
#endif

GR_INLINE truth_t gr_is_zero(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_ZERO)(x, ctx); }
GR_INLINE truth_t gr_is_one(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_ONE)(x, ctx); }
GR_INLINE truth_t gr_is_neg_one(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_NEG_ONE)(x, ctx); }

GR_INLINE truth_t gr_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_PREDICATE(ctx, EQUAL)(x, y, ctx); }
GR_INLINE truth_t gr_not_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return truth_not(GR_BINARY_PREDICATE(ctx, EQUAL)(x, y, ctx)); }

GR_INLINE WARN_UNUSED_RESULT int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, NEG)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, ADD)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, ADD_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, ADD_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, ADD_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, ADD_FMPQ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER(ctx, ADD_OTHER)(res, x, y, y_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_other_add(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx) { return GR_OTHER_BINARY_OP(ctx, OTHER_ADD)(res, x, x_ctx, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, SUB)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, SUB_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, SUB_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, SUB_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, SUB_FMPQ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER(ctx, SUB_OTHER)(res, x, y, y_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_other_sub(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx) { return GR_OTHER_BINARY_OP(ctx, OTHER_SUB)(res, x, x_ctx, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, MUL)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, MUL_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, MUL_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, MUL_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, MUL_FMPQ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER(ctx, MUL_OTHER)(res, x, y, y_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_other_mul(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx) { return GR_OTHER_BINARY_OP(ctx, OTHER_MUL)(res, x, x_ctx, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, ADDMUL)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, ADDMUL_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, ADDMUL_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, ADDMUL_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, ADDMUL_FMPQ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER(ctx, ADDMUL_OTHER)(res, x, y, y_ctx, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, SUBMUL)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, SUBMUL_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, SUBMUL_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, SUBMUL_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, SUBMUL_FMPQ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER(ctx, SUBMUL_OTHER)(res, x, y, y_ctx, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, MUL_TWO)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SQR)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_mul_2exp_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, MUL_2EXP_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_2exp_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, MUL_2EXP_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_fmpz_2exp_fmpz(gr_ptr res, const fmpz_t x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ_FMPZ(ctx, SET_FMPZ_2EXP_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPZ_FMPZ(ctx, GET_FMPZ_2EXP_FMPZ)(res1, res2, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, INV)(res, x, ctx); }
GR_INLINE truth_t gr_is_invertible(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_INVERTIBLE)(x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, DIV)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, DIV_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, DIV_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, DIV_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, DIV_FMPQ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER(ctx, DIV_OTHER)(res, x, y, y_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_other_div(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx) { return GR_OTHER_BINARY_OP(ctx, OTHER_DIV)(res, x, x_ctx, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_div_nonunique(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, DIV_NONUNIQUE)(res, x, y, ctx); }
GR_INLINE truth_t gr_divides(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_PREDICATE(ctx, DIVIDES)(x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, DIVEXACT)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_divexact_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, DIVEXACT_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_divexact_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, DIVEXACT_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_divexact_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, DIVEXACT_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_divexact_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, DIVEXACT_FMPQ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_divexact_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER(ctx, DIVEXACT_OTHER)(res, x, y, y_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_other_divexact(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx) { return GR_OTHER_BINARY_OP(ctx, OTHER_DIVEXACT)(res, x, x_ctx, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_euclidean_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, EUCLIDEAN_DIV)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_euclidean_rem(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, EUCLIDEAN_REM)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_euclidean_divrem(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_BINARY_OP(ctx, EUCLIDEAN_DIVREM)(res1, res2, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_gcd(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, GCD)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_lcm(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, LCM)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_numerator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, NUMERATOR)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_denominator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DENOMINATOR)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_factor(gr_ptr c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx) { return GR_FACTOR_OP(ctx, FACTOR)(c, factors, exponents, x, flags, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_pow(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, POW)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, POW_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, POW_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, POW_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, POW_FMPQ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER(ctx, POW_OTHER)(res, x, y, y_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_other_pow(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx) { return GR_OTHER_BINARY_OP(ctx, OTHER_POW)(res, x, x_ctx, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SQRT)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, RSQRT)(res, x, ctx); }
GR_INLINE truth_t gr_is_square(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_SQUARE)(x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_floor(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, FLOOR)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_ceil(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, CEIL)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_trunc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, TRUNC)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_nint(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, NINT)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_abs(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ABS)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_i(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, I)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_conj(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, CONJ)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_re(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, RE)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_im(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, IM)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SGN)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_csgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, CSGN)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_arg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ARG)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_pos_inf(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, POS_INF)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_neg_inf(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, NEG_INF)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_uinf(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, UINF)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_undefined(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, UNDEFINED)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_unknown(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, UNKNOWN)(res, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP_GET_INT(ctx, CMP)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP_GET_INT(ctx, CMPABS)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_cmp_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER_GET_INT(ctx, CMP_OTHER)(res, x, y, y_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_cmpabs_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) { return GR_BINARY_OP_OTHER_GET_INT(ctx, CMPABS_OTHER)(res, x, y, y_ctx, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_gen(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, GEN)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_gens(gr_vec_t res, gr_ctx_t ctx) { return GR_VEC_CTX_OP(ctx, GENS)(res, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_ctx_fq_prime(fmpz_t res, gr_ctx_t ctx) { return GR_CONSTANT_OP_GET_FMPZ(ctx, CTX_FQ_PRIME)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_ctx_fq_degree(slong * res, gr_ctx_t ctx) { return GR_CONSTANT_OP_GET_SI(ctx, CTX_FQ_DEGREE)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_ctx_fq_order(fmpz_t res, gr_ctx_t ctx) { return GR_CONSTANT_OP_GET_FMPZ(ctx, CTX_FQ_ORDER)(res, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_fq_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, FQ_FROBENIUS)(res, x, e, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPZ(ctx, FQ_MULTIPLICATIVE_ORDER)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPZ(ctx, FQ_NORM)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPZ(ctx, FQ_TRACE)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT truth_t gr_fq_is_primitive(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, FQ_IS_PRIMITIVE)(x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, FQ_PTH_ROOT)(res, x, ctx); }

GR_INLINE void _gr_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx) { GR_VEC_INIT_CLEAR_OP(ctx, VEC_INIT)(vec, len, ctx); }
GR_INLINE void _gr_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx) { GR_VEC_INIT_CLEAR_OP(ctx, VEC_CLEAR)(vec, len, ctx); }
GR_INLINE void _gr_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx) { GR_VEC_SWAP_OP(ctx, VEC_SWAP)(vec1, vec2, len, ctx); }

/* todo: warn unused? */
int gr_ctx_print(gr_ctx_t ctx);
int gr_ctx_println(gr_ctx_t ctx);
int gr_print(gr_srcptr x, gr_ctx_t ctx);
int gr_println(gr_srcptr x, gr_ctx_t ctx);
int gr_ctx_get_str(char ** s, gr_ctx_t ctx);
int gr_get_str(char ** s, gr_srcptr x, gr_ctx_t ctx);
int gr_get_str_n(char ** s, gr_srcptr x, slong n, gr_ctx_t ctx);

/* Macros for allocating temporary variables on the stack. */
/* todo: use vector init/clear functions when provided */

#define GR_TMP_VEC_ALLOC_MAX_STACK 1024
#define GR_TMP_ALLOC(size) (((size) <= GR_TMP_VEC_ALLOC_MAX_STACK) ? alloca(size) : flint_malloc(size))
#define GR_TMP_FREE(ptr, size) do { if ((size) > GR_TMP_VEC_ALLOC_MAX_STACK) flint_free(ptr); } while (0)
#define GR_TMP_ALLOC_SMALL(size) alloca(size)

#define GR_TMP_INIT_VEC(vec, len, ctx) \
    do { \
        gr_method_vec_init_clear_op vec_init = GR_VEC_INIT_CLEAR_OP(ctx, VEC_INIT); \
        slong _gr_elem_size = (ctx)->sizeof_elem; \
        (vec) = (gr_ptr) GR_TMP_ALLOC((len) * _gr_elem_size); \
        vec_init((vec), (len), (ctx)); \
    } while (0)

#define GR_TMP_CLEAR_VEC(vec, len, ctx) \
    do { \
        gr_method_vec_init_clear_op vec_clear = GR_VEC_INIT_CLEAR_OP(ctx, VEC_CLEAR); \
        slong _gr_elem_size = (ctx)->sizeof_elem; \
        vec_clear((vec), (len), (ctx)); \
        GR_TMP_FREE(vec, (len) * _gr_elem_size); \
    } while (0)

#define GR_TMP_INIT(x1, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        slong _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(1 * _gr_elem_size); \
        init(x1, (ctx)); \
    } while (0)

#define GR_TMP_INIT2(x1, x2, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        slong _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(2 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
    } while (0)

#define GR_TMP_INIT3(x1, x2, x3, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        slong _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(3 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
        init(x3, (ctx)); \
    } while (0)

#define GR_TMP_INIT4(x1, x2, x3, x4, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        slong _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(4 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        x4 = (gr_ptr) ((char *) x3 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
        init(x3, (ctx)); \
        init(x4, (ctx)); \
    } while (0)

#define GR_TMP_INIT5(x1, x2, x3, x4, x5, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        slong _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(5 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        x4 = (gr_ptr) ((char *) x3 + _gr_elem_size); \
        x5 = (gr_ptr) ((char *) x4 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
        init(x3, (ctx)); \
        init(x4, (ctx)); \
        init(x5, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR(x1, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR2(x1, x2, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR3(x1, x2, x3, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR4(x1, x2, x3, x4, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
        clear(x4, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR5(x1, x2, x3, x4, x5, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
        clear(x4, (ctx)); \
        clear(x5, (ctx)); \
    } while (0)

GR_INLINE gr_ptr gr_heap_init(gr_ctx_t ctx)
{
    gr_ptr ptr;
    ptr = (gr_ptr) flint_malloc(ctx->sizeof_elem);
    gr_init(ptr, ctx);
    return ptr;
}

GR_INLINE void gr_heap_clear(gr_ptr x, gr_ctx_t ctx)
{
    gr_clear(x, ctx);
    flint_free(x);
}

GR_INLINE gr_ptr gr_heap_init_vec(slong len, gr_ctx_t ctx)
{
    gr_ptr ptr;
    ptr = (gr_ptr) flint_malloc(len * ctx->sizeof_elem);
    _gr_vec_init(ptr, len, ctx);
    return ptr;
}

GR_INLINE void gr_heap_clear_vec(gr_ptr x, slong len, gr_ctx_t ctx)
{
    _gr_vec_clear(x, len, ctx);
    flint_free(x);
}

truth_t gr_generic_ctx_predicate(gr_ctx_t ctx);
truth_t gr_generic_ctx_predicate_true(gr_ctx_t ctx);
truth_t gr_generic_ctx_predicate_false(gr_ctx_t ctx);

/* Some base rings */

void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state);

void gr_ctx_init_fmpz(gr_ctx_t ctx);
void gr_ctx_init_fmpq(gr_ctx_t ctx);
void gr_ctx_init_fmpzi(gr_ctx_t ctx);

void gr_ctx_init_fmpz_mod(gr_ctx_t ctx, const fmpz_t n);
void _gr_ctx_init_fmpz_mod_from_ref(gr_ctx_t ctx, const void * fmod_ctx);
void gr_ctx_fmpz_mod_set_primality(gr_ctx_t ctx, truth_t is_prime);

void gr_ctx_init_nmod(gr_ctx_t ctx, ulong n);
void _gr_ctx_init_nmod(gr_ctx_t ctx, void * nmod_t_ref);

void gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n);
void gr_ctx_init_nmod32(gr_ctx_t ctx, unsigned int n);

void gr_ctx_init_real_qqbar(gr_ctx_t ctx);
void gr_ctx_init_complex_qqbar(gr_ctx_t ctx);

void gr_ctx_init_real_arb(gr_ctx_t ctx, slong prec);
void gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec);

void gr_ctx_init_real_float_arf(gr_ctx_t ctx, slong prec);
void gr_ctx_init_complex_float_acf(gr_ctx_t ctx, slong prec);


void gr_ctx_init_real_ca(gr_ctx_t ctx);
void gr_ctx_init_complex_ca(gr_ctx_t ctx);
void gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx);
void gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx);
void gr_ctx_init_complex_extended_ca(gr_ctx_t ctx);
void _gr_ctx_init_ca_from_ref(gr_ctx_t ctx, int which_ring, void * ca_ctx);

void gr_ctx_init_fq(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var);
void gr_ctx_init_fq_nmod(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var);
void gr_ctx_init_fq_zech(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var);

void _gr_ctx_init_fq_from_ref(gr_ctx_t ctx, const void * fq_ctx);
void _gr_ctx_init_fq_nmod_from_ref(gr_ctx_t ctx, const void * fq_nmod_ctx);
void _gr_ctx_init_fq_zech_from_ref(gr_ctx_t ctx, const void * fq_zech_ctx);

void gr_ctx_init_fmpz_poly(gr_ctx_t ctx);
void gr_ctx_init_fmpq_poly(gr_ctx_t ctx);

#ifdef FMPQ_POLY_H
void gr_ctx_init_nf(gr_ctx_t ctx, const fmpq_poly_t poly);
void gr_ctx_init_nf_fmpz_poly(gr_ctx_t ctx, const fmpz_poly_t poly);
void _gr_ctx_init_nf_from_ref(gr_ctx_t ctx, const void * nfctx);
#endif

/* Groups */

void gr_ctx_init_perm(gr_ctx_t ctx, ulong n);
void gr_ctx_init_psl2z(gr_ctx_t ctx);
int gr_ctx_init_dirichlet_group(gr_ctx_t ctx, ulong q);

/* Generic polynomial ring */

typedef struct
{
    gr_ctx_struct * base_ring;
    slong degree_limit;
    char * var;
}
polynomial_ctx_t;

#define POLYNOMIAL_CTX(ring_ctx) ((polynomial_ctx_t *)((ring_ctx)))
#define POLYNOMIAL_ELEM_CTX(ring_ctx) (POLYNOMIAL_CTX(ring_ctx)->base_ring)

void gr_ctx_init_gr_poly(gr_ctx_t ctx, gr_ctx_t base_ring);

/* Multivariate */

#ifdef MPOLY_H
void gr_ctx_init_fmpz_mpoly(gr_ctx_t ctx, slong nvars, const ordering_t ord);
void gr_ctx_init_gr_mpoly(gr_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord);
void gr_ctx_init_fmpz_mpoly_q(gr_ctx_t ctx, slong nvars, const ordering_t ord);
#endif

/* Generic series */

void gr_ctx_init_gr_series(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec);
void gr_ctx_init_gr_series_mod(gr_ctx_t ctx, gr_ctx_t base_ring, slong mod);

/* Generic vectors */

typedef struct
{
    gr_ctx_struct * base_ring;
    int all_sizes;
    slong n;
}
vector_ctx_t;

#define VECTOR_CTX(ring_ctx) ((vector_ctx_t *)((ring_ctx)))

void gr_ctx_init_vector_gr_vec(gr_ctx_t ctx, gr_ctx_t base_ring);
void gr_ctx_init_vector_space_gr_vec(gr_ctx_t ctx, gr_ctx_t base_ring, slong n);

/* Generic matrix ring */

typedef struct
{
    gr_ctx_struct * base_ring;
    int all_sizes;
    slong nrows;
    slong ncols;
}
matrix_ctx_t;

#define MATRIX_CTX(ring_ctx) ((matrix_ctx_t *)((ring_ctx)))

void gr_ctx_init_matrix_domain(gr_ctx_t ctx, gr_ctx_t base_ring);
void gr_ctx_init_matrix_space(gr_ctx_t ctx, gr_ctx_t base_ring, slong nrows, slong ncols);

GR_INLINE void gr_ctx_init_matrix_ring(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)
{
    gr_ctx_init_matrix_space(ctx, base_ring, n, n);
}

/* Expressions */

void gr_ctx_init_fexpr(gr_ctx_t ctx);

/* Coercions */

int gr_ctx_cmp_coercion(gr_ctx_t ctx1, gr_ctx_t ctx2);

/* Testing */

#define GR_TEST_VERBOSE 8
#define GR_TEST_ALWAYS_ABLE 16

/* todo: just have gr_test_structure() */
void gr_test_ring(gr_ctx_t R, slong iters, int test_flags);
void gr_test_multiplicative_group(gr_ctx_t R, slong iters, int test_flags);

#ifdef __cplusplus
}
#endif

#endif
