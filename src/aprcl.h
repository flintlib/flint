/*
    Copyright (C) 2015, Vladimir Glazchev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef APRCL_H
#define APRCL_H

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define SQUARING_SPACE 70

/* Configuration struct */
typedef struct
{
    ulong R;
    fmpz_t s;
    n_factor_t rs;
    fmpz_factor_t qs;
    int * qs_used;
} _aprcl_config;

typedef _aprcl_config aprcl_config[1];

/* Z[unity_root_q, unity_root_p]/(n) struct */
typedef struct
{
    fmpz_mod_poly_t *polys;
    ulong p;
    ulong q;
    fmpz_mod_ctx_t ctx;
} _unity_zpq;

typedef _unity_zpq unity_zpq[1];

/* Z[unity_root]/(n) struct */
typedef struct
{
    fmpz_mod_poly_t poly;
    ulong p;
    ulong exp;
    fmpz_mod_ctx_t ctx;
} _unity_zp;

typedef _unity_zp unity_zp[1];

/* Primality test status */
typedef enum
{
    UNKNOWN,
    PRIME,
    COMPOSITE,
    PROBABPRIME
} primality_test_status;

/* Useful functions */
FLINT_DLL int _aprcl_p_ind(const aprcl_config conf, ulong p);

FLINT_DLL ulong aprcl_p_power_in_q(ulong q, ulong p);

FLINT_DLL int aprcl_is_mul_coprime_ui_ui(ulong x, ulong y, const fmpz_t n);

FLINT_DLL int aprcl_is_mul_coprime_ui_fmpz(ulong x, const fmpz_t y, const fmpz_t n);

/* 
                            Primality tests
--------------------------------------------------------------------------------
*/

FLINT_DLL int aprcl_is_prime(const fmpz_t n);

/* Gauss test configuration */
FLINT_DLL void aprcl_config_gauss_init(aprcl_config conf, const fmpz_t n);

FLINT_DLL void aprcl_config_gauss_init_min_R(aprcl_config conf,
        const fmpz_t n, ulong R);

FLINT_DLL void aprcl_config_gauss_clear(aprcl_config conf);

/* Jacobi test configuration */
FLINT_DLL ulong aprcl_R_value(const fmpz_t n);

FLINT_DLL void aprcl_config_jacobi_init(aprcl_config conf, const fmpz_t n);

FLINT_DLL void aprcl_config_jacobi_clear(aprcl_config conf);

/*  Gauss sums primality test */
FLINT_DLL int aprcl_is_prime_gauss(const fmpz_t n);

FLINT_DLL int aprcl_is_prime_gauss_min_R(const fmpz_t n, ulong R);

FLINT_DLL primality_test_status _aprcl_is_prime_gauss(const fmpz_t n,
        const aprcl_config config);

FLINT_DLL int _aprcl_is_gausspower_2q_equal_first(ulong q, const fmpz_t n);

FLINT_DLL int _aprcl_is_gausspower_2q_equal_second(ulong q, const fmpz_t n);

FLINT_DLL slong _aprcl_is_gausspower_from_unity_p(ulong q, ulong r, const fmpz_t n);

/* Jacobi sums primality test */
FLINT_DLL int aprcl_is_prime_jacobi(const fmpz_t n);

FLINT_DLL primality_test_status _aprcl_is_prime_jacobi(const fmpz_t n,
        const aprcl_config config);

FLINT_DLL slong _aprcl_is_prime_jacobi_check_pk(const unity_zp j,
        const fmpz_t u, ulong v);

FLINT_DLL slong _aprcl_is_prime_jacobi_check_21(ulong q,
        const fmpz_t n);

FLINT_DLL slong _aprcl_is_prime_jacobi_check_22(const unity_zp j,
        const fmpz_t u, ulong v, ulong q);

FLINT_DLL slong _aprcl_is_prime_jacobi_check_2k(const unity_zp j, const unity_zp j2_1,
        const unity_zp j2_2, const fmpz_t u, ulong v);

FLINT_DLL int _aprcl_is_prime_jacobi_additional_test(const fmpz_t n, ulong p);

/* Final division functions */
FLINT_DLL int aprcl_is_prime_divisors_in_residue(const fmpz_t n,
        const fmpz_t s, ulong r);

FLINT_DLL int aprcl_is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r);

/*
                        Z[unity_root]/(n) operations
--------------------------------------------------------------------------------
*/

/* Memory management */
FLINT_DLL void unity_zp_init(unity_zp f, ulong p, ulong exp, const fmpz_t n);

FLINT_DLL void unity_zp_clear(unity_zp f);

FLINT_DLL void unity_zp_copy(unity_zp f, const unity_zp g);

FLINT_DLL void unity_zp_swap(unity_zp f, unity_zp g);

FLINT_DLL void unity_zp_set_zero(unity_zp f);

/* Comparison */
FLINT_DLL slong unity_zp_is_unity(unity_zp f);

FLINT_DLL int unity_zp_equal(unity_zp f, unity_zp g);

/* Output */
FLINT_DLL void unity_zp_print(const unity_zp f);

/* Coefficient management */
FLINT_DLL void unity_zp_coeff_set_fmpz(unity_zp f, ulong ind, const fmpz_t x);

FLINT_DLL void unity_zp_coeff_set_ui(unity_zp f, ulong ind, ulong x);

FLINT_DLL void unity_zp_coeff_add_fmpz(unity_zp f, ulong ind, const fmpz_t x);

FLINT_DLL void unity_zp_coeff_add_ui(unity_zp f, ulong ind, ulong x);

FLINT_DLL void unity_zp_coeff_inc(unity_zp f, ulong ind);

FLINT_DLL void unity_zp_coeff_dec(unity_zp f, ulong ind);

/* Scalar multiplication */
FLINT_DLL void unity_zp_mul_scalar_fmpz(unity_zp f,
        const unity_zp g, const fmpz_t s);

FLINT_DLL void unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s);

/* Addition */
FLINT_DLL void unity_zp_add(unity_zp f, const unity_zp g, const unity_zp h);

/* General multiplication and squaring */
FLINT_DLL void unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h);

FLINT_DLL void unity_zp_sqr(unity_zp f, const unity_zp g);

/* Special multiplication and squaring */
FLINT_DLL void unity_zp_mul_inplace(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_sqr_inplace(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_ar1(fmpz_t * t);

FLINT_DLL void unity_zp_ar2(fmpz_t * t);

FLINT_DLL void unity_zp_ar3(fmpz_t * t);

FLINT_DLL void unity_zp_ar4(fmpz_t * t);

FLINT_DLL void unity_zp_mul3(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_mul9(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_mul4(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_mul8(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_mul16(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_mul5(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_mul7(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_mul11(unity_zp f,
        const unity_zp g, const unity_zp h, fmpz_t * t);

FLINT_DLL void unity_zp_sqr3(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr9(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr4(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr8(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr16(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr5(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr7(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr11(unity_zp f, const unity_zp g, fmpz_t * t);

/* Powering functions */
FLINT_DLL void unity_zp_pow_fmpz(unity_zp f,
        const unity_zp g, const fmpz_t pow);

FLINT_DLL void unity_zp_pow_ui(unity_zp f, const unity_zp g, ulong pow);

FLINT_DLL ulong _unity_zp_pow_select_k(const fmpz_t n);

FLINT_DLL void unity_zp_pow_2k_fmpz(unity_zp f,
        const unity_zp g, const fmpz_t pow);

FLINT_DLL void unity_zp_pow_2k_ui(unity_zp f, const unity_zp g, ulong pow);

FLINT_DLL void unity_zp_pow_sliding_fmpz(unity_zp f,
        unity_zp g, const fmpz_t pow);

/* Cyclotomic reduction */
FLINT_DLL void _unity_zp_reduce_cyclotomic_divmod(unity_zp f);

FLINT_DLL void _unity_zp_reduce_cyclotomic(unity_zp f);

FLINT_DLL void unity_zp_reduce_cyclotomic(unity_zp f, const unity_zp g);

/* Automorphism and inverse computation */
FLINT_DLL void unity_zp_aut(unity_zp f, const unity_zp g, ulong x);

FLINT_DLL void unity_zp_aut_inv(unity_zp f, const unity_zp g, ulong x);

/* Jacobi sum computation. */
FLINT_DLL mp_ptr aprcl_f_table(const ulong q);

FLINT_DLL void _unity_zp_jacobi_sum_pq_general(unity_zp f,
        const mp_ptr table, ulong p, ulong q, ulong k, ulong a, ulong b);

FLINT_DLL void unity_zp_jacobi_sum_pq(unity_zp f, ulong q, ulong p);

FLINT_DLL void unity_zp_jacobi_sum_2q_one(unity_zp f, ulong q);

FLINT_DLL void unity_zp_jacobi_sum_2q_two(unity_zp f, ulong q);


/* 
                Z[unity_root_q, unity_root_p]/(n) operations 
--------------------------------------------------------------------------------
*/

/* Memory management */
FLINT_DLL void unity_zpq_init(unity_zpq f, ulong q, ulong p, const fmpz_t n);

FLINT_DLL void unity_zpq_clear(unity_zpq f);

FLINT_DLL void unity_zpq_copy(unity_zpq f, const unity_zpq g);

FLINT_DLL void unity_zpq_swap(unity_zpq f, unity_zpq g);

/* Comparison */
FLINT_DLL int unity_zpq_equal(const unity_zpq f, const unity_zpq g);

FLINT_DLL slong unity_zpq_p_unity(const unity_zpq f);

FLINT_DLL int unity_zpq_is_p_unity(const unity_zpq f);

FLINT_DLL int unity_zpq_is_p_unity_generator(const unity_zpq f);

/* Coefficient management */
FLINT_DLL void unity_zpq_coeff_set_fmpz(unity_zpq f,
        slong i, slong j, const fmpz_t x);

FLINT_DLL void unity_zpq_coeff_set_ui(unity_zpq f,
        slong i, slong j, ulong x);

FLINT_DLL void unity_zpq_coeff_add(unity_zpq f,
        slong i, slong j, const fmpz_t x);

FLINT_DLL void unity_zpq_coeff_add_ui(unity_zpq f,
        slong i, slong j, ulong x);

/* Scalar multiplication */
FLINT_DLL void unity_zpq_scalar_mul_ui(unity_zpq f,
        const unity_zpq g, ulong s);

/* Addition and multiplication */
FLINT_DLL void unity_zpq_add(unity_zpq f,
        const unity_zpq g, const unity_zpq h);

FLINT_DLL void unity_zpq_mul(unity_zpq f,
        const unity_zpq g, const unity_zpq h);

FLINT_DLL void _unity_zpq_mul_unity_p(unity_zpq f);

FLINT_DLL void unity_zpq_mul_unity_p_pow(
        unity_zpq f, const unity_zpq g, slong p);

/* Powering */
FLINT_DLL void unity_zpq_pow(unity_zpq f,
        const unity_zpq g, const fmpz_t p);

FLINT_DLL void unity_zpq_pow_ui(unity_zpq f,
        const unity_zpq g, ulong pow);

/* Gauss sum computation */
FLINT_DLL void unity_zpq_gauss_sum(unity_zpq f,
        ulong q, ulong p);

FLINT_DLL void unity_zpq_gauss_sum_character_pow(unity_zpq f,
        ulong q, ulong p, ulong pow);

FLINT_DLL void unity_zpq_gauss_sum_sigma_pow(
        unity_zpq f, ulong q, ulong p);

#ifdef __cplusplus
}
#endif

#endif

