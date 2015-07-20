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

#define SQUARING_SPACE 50

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

    ulong r;
    fmpz_t ninv;
} _unity_zp;

typedef _unity_zp unity_zp[1];

typedef struct
{
    fmpz_poly_t poly;
    ulong p;
    ulong exp;
    fmpz_t n;
    fmpz_t nr;

    fmpz_t ninv;
    ulong r;
} _unity_zp_mont;

typedef _unity_zp_mont unity_zp_mont[1];

typedef enum
{
    UNKNOWN,
    PRIME,
    COMPOSITE,
    PROBABPRIME
} primality_test_status;

/* Primality testing useful functions */
int _p_ind(const aprcl_config conf, ulong p);

int is_mul_coprime_ui_ui(ulong x, ulong y, const fmpz_t n);
int is_mul_coprime_ui_fmpz(ulong x, const fmpz_t y, const fmpz_t n);

int is_prime_divisors_in_residue(const fmpz_t n, const fmpz_t s, ulong r);
int is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r);

/* Gauss sums primality test */

int is_prime_gauss(const fmpz_t n);
int is_prime_gauss_min_R(const fmpz_t n, ulong R);

primality_test_status _is_prime_gauss(const fmpz_t n,
        const aprcl_config config);

int _is_gausspower_2q_equal_first(ulong q, const fmpz_t n);
int _is_gausspower_2q_equal_second(ulong q, const fmpz_t n);

slong _is_gausspower_from_unity_p(ulong q, ulong r, const fmpz_t n);


/* Jacobi sums primality test */

int is_prime_jacobi(const fmpz_t n);
primality_test_status _is_prime_jacobi(const fmpz_t n,
        const aprcl_config config);

slong _is_prime_jacobi_check_pk(const unity_zp j, const fmpz_t u, ulong v);
slong _is_prime_jacobi_check_21(ulong q, const fmpz_t n);
slong _is_prime_jacobi_check_22(const unity_zp j
    , const fmpz_t u, ulong v, ulong q);
slong _is_prime_jacobi_check_2k(const unity_zp j, const unity_zp j2_1
    , const unity_zp j2_2, const fmpz_t u, ulong v);

int _is_prime_jacobi_additional_test(const fmpz_t n, ulong p);


int is_prime_aprcl(const fmpz_t n); /* now just returns is_prime_jacobi(n) */

/* Z[unity_root]/(n) operation */

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
FLINT_DLL void unity_zp_sqr(unity_zp f, const unity_zp g);
FLINT_DLL void unity_zp_sqr_inplace(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_ar1(fmpz_t * t);
FLINT_DLL void unity_zp_ar2(fmpz_t * t);
FLINT_DLL void unity_zp_ar3(fmpz_t * t);
FLINT_DLL void unity_zp_ar4(fmpz_t * t);

FLINT_DLL void unity_zp_sqr3(unity_zp f, const unity_zp g, fmpz_t * t);
FLINT_DLL void unity_zp_sqr9(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr4(unity_zp f, const unity_zp g, fmpz_t * t);
FLINT_DLL void unity_zp_sqr8(unity_zp f, const unity_zp g, fmpz_t * t);
FLINT_DLL void unity_zp_sqr16(unity_zp f, const unity_zp g, fmpz_t * t);

FLINT_DLL void unity_zp_sqr5(unity_zp f, const unity_zp g, fmpz_t * t);
FLINT_DLL void unity_zp_sqr7(unity_zp f, const unity_zp g, fmpz_t * t);

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
void _jacobi_pq_general(unity_zp f, const mp_ptr table, ulong p,
        ulong q, ulong k, ulong a, ulong b);
void jacobi_pq(unity_zp f, ulong q, ulong p);
void jacobi_2q_one(unity_zp f, ulong q);
void jacobi_2q_two(unity_zp f, ulong q);

/* Z[unity_root_q, unity_root_p] operations. */
void unity_zpq_init(unity_zpq value, ulong q, ulong p, const fmpz_t n);
void unity_zpq_clear(unity_zpq value);

void unity_zpq_copy(unity_zpq f, const unity_zpq g);
void unity_zpq_swap(unity_zpq f, unity_zpq g);

void unity_zpq_coeff_set_fmpz(unity_zpq value,
        ulong i, ulong j, const fmpz_t x);
void unity_zpq_coeff_set_ui(unity_zpq value, ulong i, ulong j, ulong x);

void unity_zpq_coeff_add(unity_zpq value, ulong i, ulong j, const fmpz_t x);
void unity_zpq_coeff_add_ui(unity_zpq value, ulong i, ulong j, ulong x);

int unity_zpq_equal(const unity_zpq f, const unity_zpq g);

void unity_zpq_add(unity_zpq result, const unity_zpq left,
        const unity_zpq right);
void unity_zpq_mul(unity_zpq result, const unity_zpq left,
        const unity_zpq right);
void unity_zpq_pow(unity_zpq f, const unity_zpq g, const fmpz_t p);
void unity_zpq_pow_ui(unity_zpq f, const unity_zpq g, ulong pow);

void _unity_zpq_mul_unity_p(unity_zpq f);
void unity_zpq_mul_unity_p_pow(unity_zpq f, const unity_zpq g, ulong p);

void unity_zpq_scalar_mul_ui(unity_zpq f, const unity_zpq g, ulong s);

void unity_zpq_gauss_sum(unity_zpq value, ulong q, ulong p);
void unity_zpq_gauss_sum_character_pow(unity_zpq value,
        ulong q, ulong p, ulong pow);
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


/* Z[unity_root] functions for Montgomery form */

/* Memory allocation and free */
void unity_zp_mont_init(unity_zp_mont f, ulong p, ulong exp,
        const fmpz_t n, const fmpz_t ninv);
void unity_zp_mont_clear(unity_zp_mont f);

void unity_zp_mont_set_zero(unity_zp_mont f);

void unity_zp_mont_swap(unity_zp_mont f, unity_zp_mont g);
void unity_zp_mont_copy(unity_zp_mont f, const unity_zp_mont g);

/* Conversion from simple to Montgomery form and back */
void unity_zp_to_mont(unity_zp_mont f, const unity_zp g);
void unity_zp_from_mont(unity_zp f, unity_zp_mont g);

/* Montgomery reduction */
void unity_zp_mont_ninv(fmpz_t ninv, const unity_zp_mont f);
void unity_zp_mont_reduction(unity_zp_mont f);


/* Modular functions for fmpz_t */

static __inline__
void
mod_redc_mont(fmpz_t f, const fmpz_t g, const fmpz_t n,
        const fmpz_t ninv, fmpz_t m, fmpz_t t, ulong r)
{
    fmpz_fdiv_r_2exp(m, g, r);
    fmpz_mul(t, m, ninv);
    fmpz_fdiv_r_2exp(m, t, r);
    fmpz_mul(t, n, m);
    fmpz_add(f, g, t);
    fmpz_fdiv_q_2exp(f, f, r);
    if (fmpz_cmp(n, f) <= 0)
        fmpz_sub(f, f, n);
}

static __inline__
void
mod_mul(fmpz_t f, const fmpz_t g, const fmpz_t h, const fmpz_t n,
        const fmpz_t ninv, fmpz_t m, fmpz_t t, ulong r)
{
    fmpz_mul(f, g, h);
    mod_redc_mont(f, f, n, ninv, m, t, r);
}

static __inline__
void
mod_add(fmpz_t f, const fmpz_t g, const fmpz_t h, const fmpz_t n)
{
    fmpz_add(f, g, h);
    if (fmpz_cmp(f, n) >= 0)
        fmpz_sub(f, f, n);
}

static __inline__
void
mod_sub(fmpz_t f, const fmpz_t g, const fmpz_t h, const fmpz_t n)
{
    fmpz_sub(f, g, h);
    if (fmpz_cmp_ui(f, 0) < 0)
        fmpz_add(f, f, n);
}

/* Cyclotomic reduction */
void _unity_zp_mont_reduce_cyclotomic(unity_zp_mont f);

/* Multiplication and squaring */
void unity_zp_mont_sqr(unity_zp_mont f, const unity_zp_mont g);
void unity_zp_mont_mul(unity_zp_mont f, const unity_zp_mont g,
        const unity_zp_mont h);

void
unity_zp_mont_sqr7(unity_zp_mont f, const unity_zp_mont g, fmpz_t * t);

void unity_zp_mont_sqr_inplace(unity_zp_mont f,
        const unity_zp_mont g, fmpz_t * t);

/* Powering functions */
int unity_zp_is_inplace(ulong p, ulong exp);
void _unity_zp_pow_mont_fmpz(unity_zp_mont f, const unity_zp_mont g,
        const fmpz_t pow, const fmpz_t r);
void unity_zp_pow_mont_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow);

#ifdef __cplusplus
}
#endif

#endif
