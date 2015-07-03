/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015, Vladimir Glazchev

******************************************************************************/

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

typedef struct
{
    ulong R;
    fmpz_t s;
    n_factor_t rs;
    fmpz_factor_t qs;
} _aprcl_config;

typedef _aprcl_config aprcl_config[1];

typedef struct
{
    fmpz_mod_poly_t *polys;
    ulong p;
    ulong q;
    fmpz_t n;
} _unity_zpq;

typedef _unity_zpq unity_zpq[1];

typedef struct
{
    fmpz_mod_poly_t poly;
    ulong p;
    ulong exp;
    fmpz_t n;
} _unity_zp;

typedef _unity_zp unity_zp[1];

typedef enum
{
    UNKNOWN,
    PRIME,
    COMPOSITE,
    PROBABPRIME
} primality_test_status;

int _p_ind(const aprcl_config conf, ulong p);
int _is_gausspower_2q_equal_first(ulong q, const fmpz_t n);
int _is_gausspower_2q_equal_second(ulong q, const fmpz_t n);
slong _is_gausspower_from_unity_p(ulong q, ulong r, const fmpz_t n);
primality_test_status _is_prime_gauss(const fmpz_t n, const aprcl_config config);
int is_prime_gauss_min_R(const fmpz_t n, ulong R);
int is_prime_gauss(const fmpz_t n);


int _is_prime_jacobi_additional_test(const fmpz_t n, ulong p);
slong _is_prime_jacobi_check_pk(const unity_zp j, const fmpz_t u, ulong v);
slong _is_prime_jacobi_check_21(ulong q, const fmpz_t n);
slong _is_prime_jacobi_check_22(const unity_zp j
    , const fmpz_t u, ulong v, ulong q);
slong _is_prime_jacobi_check_2k(const unity_zp j, const unity_zp j2_1
    , const unity_zp j2_2, const fmpz_t u, ulong v);
primality_test_status _is_prime_jacobi(const fmpz_t n, const aprcl_config config);
int is_prime_jacobi(const fmpz_t n);

int is_prime_aprcl(const fmpz_t n);

int is_mul_coprime_ui_ui(ulong x, ulong y, const fmpz_t n);
int is_mul_coprime_ui_fmpz(ulong x, const fmpz_t y, const fmpz_t n);

int is_prime_divisors_in_residue(const fmpz_t n, const fmpz_t s, ulong r);
int is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r);

/* Z[unity_root] operation v2. */

void unity_zp_init(unity_zp f, ulong p, ulong exp, const fmpz_t n);
void unity_zp_clear(unity_zp f);
void unity_zp_copy(unity_zp f, const unity_zp g);
void unity_zp_swap(unity_zp f, unity_zp g);

void unity_zp_set_zero(unity_zp f);

slong unity_zp_is_unity(unity_zp f);
int unity_zp_equal(unity_zp f, unity_zp g);

void unity_zp_print(const unity_zp f);

void unity_zp_coeff_inc(unity_zp f, ulong ind);
void unity_zp_coeff_dec(unity_zp f, ulong ind);
void unity_zp_coeff_set_fmpz(unity_zp value, ulong ind, const fmpz_t x);
void unity_zp_coeff_set_ui(unity_zp value, ulong ind, ulong x);
void unity_zp_coeff_add_fmpz(unity_zp f, ulong ind, const fmpz_t x);
void unity_zp_coeff_add_ui(unity_zp f, ulong ind, ulong x);

void unity_zp_mul_scalar_fmpz(unity_zp f, const unity_zp g, const fmpz_t s);
void unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s);

void unity_zp_add(unity_zp f, const unity_zp g, const unity_zp h);
void unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h);
void unity_zp_sqr(unity_zp f, const unity_zp g);
void unity_zp_sqr_inplace(unity_zp f, const unity_zp g, fmpz_t * t);

void unity_zp_sqr5(unity_zp f, const unity_zp g, fmpz_t * t);
void unity_zp_sqr7(unity_zp f, const unity_zp g, fmpz_t * t);

ulong _unity_zp_pow_2k_find_k(const fmpz_t n);
void unity_zp_pow_2k_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow);
void unity_zp_pow_2k_ui(unity_zp f, const unity_zp g, ulong pow);

void unity_zp_pow_sliding_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow);

void unity_zp_pow_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow);
void unity_zp_pow_ui(unity_zp f, const unity_zp g, ulong pow);

void _unity_zp_reduce_cyclotomic(unity_zp f);
void unity_zp_reduce_cyclotomic(unity_zp f, const unity_zp g);

void unity_zp_aut(unity_zp f, const unity_zp g, ulong x);
void unity_zp_aut_inv(unity_zp f, const unity_zp g, ulong x);

/* Jacobi sum computation. */
void _jacobi_pq_general(unity_zp f, const mp_ptr table, ulong p, ulong q, ulong k, ulong a, ulong b);
void jacobi_pq(unity_zp f, ulong q, ulong p);
void jacobi_2q_one(unity_zp f, ulong q);
void jacobi_2q_two(unity_zp f, ulong q);

/* Z[unity_root_q, unity_root_p] operations. */
void unity_zpq_init(unity_zpq value, ulong q, ulong p, const fmpz_t n);
void unity_zpq_clear(unity_zpq value);

void unity_zpq_copy(unity_zpq f, const unity_zpq g);
void unity_zpq_swap(unity_zpq f, unity_zpq g);

void unity_zpq_coeff_set_fmpz(unity_zpq value, ulong i, ulong j, const fmpz_t x);
void unity_zpq_coeff_set_ui(unity_zpq value, ulong i, ulong j, ulong x);

void unity_zpq_coeff_add(unity_zpq value, ulong i, ulong j, const fmpz_t x);
void unity_zpq_coeff_add_ui(unity_zpq value, ulong i, ulong j, ulong x);

int unity_zpq_equal(const unity_zpq f, const unity_zpq g);

void unity_zpq_add(unity_zpq result, const unity_zpq left, const unity_zpq right);
void unity_zpq_mul(unity_zpq result, const unity_zpq left, const unity_zpq right);
void unity_zpq_pow(unity_zpq f, const unity_zpq g, const fmpz_t p);
void unity_zpq_pow_ui(unity_zpq f, const unity_zpq g, ulong pow);

void _unity_zpq_mul_unity_p(unity_zpq f);
void unity_zpq_mul_unity_p_pow(unity_zpq f, const unity_zpq g, ulong p);

void unity_zpq_scalar_mul_ui(unity_zpq f, const unity_zpq g, ulong s);

void unity_zpq_gauss_sum(unity_zpq value, ulong q, ulong p);
void unity_zpq_gauss_sum_character_pow(unity_zpq value, ulong q, ulong p, ulong pow);
void unity_zpq_gauss_sum_sigma_pow(unity_zpq value, ulong q, ulong p);

ulong unity_zpq_p_unity(const unity_zpq f);
int unity_zpq_is_p_unity(const unity_zpq f);
int unity_zpq_is_p_unity_generator(const unity_zpq f);


/* Initial step functions. */
void aprcl_config_init(aprcl_config conf, const fmpz_t n);
void aprcl_config_init_min_R(aprcl_config conf, const fmpz_t n, ulong R);
void aprcl_config_clear(aprcl_config conf);
void _aprcl_config_update(aprcl_config conf);

ulong _R_value(const fmpz_t n);
void _jacobi_config_update(aprcl_config conf);
void jacobi_config_init(aprcl_config conf, const fmpz_t n);
void jacobi_config_clear(aprcl_config conf);

mp_ptr f_table(const ulong q);
ulong p_power_in_q(ulong q, ulong p);

#ifdef __cplusplus
}
#endif

#endif
