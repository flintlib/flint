/*
    Copyright (C) 2007, David Howden
    Copyright (C) 2010, 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_H
#define NMOD_POLY_H

#ifdef NMOD_POLY_INLINES_C
#define NMOD_POLY_INLINE
#else
#define NMOD_POLY_INLINE static inline
#endif

#include "nmod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NMOD_POLY_HGCD_CUTOFF  100      /* HGCD: Basecase -> Recursion      */
#define NMOD_POLY_GCD_CUTOFF  340       /* GCD:  Euclidean -> HGCD          */
#define NMOD_POLY_SMALL_GCD_CUTOFF 200  /* GCD (small n): Euclidean -> HGCD */

typedef struct
{
   ulong res;
   ulong lc;
   slong len0;
   slong len1;
   slong off;
} nmod_poly_res_struct;

typedef nmod_poly_res_struct nmod_poly_res_t[1];

typedef struct
{
    nmod_mat_struct * A;
    nmod_poly_struct * poly1;
    nmod_poly_struct * poly2;
    nmod_poly_struct * poly2inv;
}
nmod_poly_matrix_precompute_arg_t;

typedef struct
{
    nmod_mat_struct * A;
    nmod_poly_struct * res;
    nmod_poly_struct * poly1;
    nmod_poly_struct * poly3;
    nmod_poly_struct * poly3inv;
}
nmod_poly_compose_mod_precomp_preinv_arg_t;

/* Memory management  ********************************************************/

void nmod_poly_init(nmod_poly_t poly, ulong n);
void nmod_poly_init_preinv(nmod_poly_t poly, ulong n, ulong ninv);

void nmod_poly_init2(nmod_poly_t poly, ulong n, slong alloc);
void nmod_poly_init2_preinv(nmod_poly_t poly, ulong n, ulong ninv, slong alloc);

void nmod_poly_realloc(nmod_poly_t poly, slong alloc);

void nmod_poly_clear(nmod_poly_t poly);

void nmod_poly_fit_length(nmod_poly_t poly, slong alloc);

NMOD_POLY_INLINE
void nmod_poly_init_mod(nmod_poly_t poly, const nmod_t mod)
{
    poly->coeffs = NULL;
    poly->alloc = 0;
    poly->length = 0;
    poly->mod = mod;
}

NMOD_POLY_INLINE
void nmod_poly_set_mod(nmod_poly_t poly, const nmod_t mod)
{
    poly->mod = mod;
}

NMOD_POLY_INLINE
void _nmod_poly_set_length(nmod_poly_t poly, slong len)
{
    poly->length = len;
}

NMOD_POLY_INLINE
void _nmod_poly_normalise(nmod_poly_t poly)
{
    while (poly->length && (poly->coeffs[poly->length - 1] == WORD(0)))
        poly->length--;
}

/* Polynomial parameters  ****************************************************/

NMOD_POLY_INLINE
slong nmod_poly_length(const nmod_poly_t poly)
{
    return poly->length;
}

NMOD_POLY_INLINE
slong nmod_poly_degree(const nmod_poly_t poly)
{
    return poly->length - 1;
}

NMOD_POLY_INLINE
ulong nmod_poly_modulus(const nmod_poly_t poly)
{
    return poly->mod.n;
}

flint_bitcnt_t nmod_poly_max_bits(const nmod_poly_t poly);

NMOD_POLY_INLINE
nn_ptr nmod_poly_lead(const nmod_poly_t poly)
{
    if (poly->length)
        return poly->coeffs + (poly->length - 1);
    else
        return NULL;
}

/* Assignment and basic manipulation  ****************************************/

void nmod_poly_set(nmod_poly_t a, const nmod_poly_t b);

NMOD_POLY_INLINE
void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)
{
    FLINT_SWAP(nmod_poly_struct, *poly1, *poly2);
}

NMOD_POLY_INLINE
void nmod_poly_zero(nmod_poly_t res)
{
    res->length = 0;
}

NMOD_POLY_INLINE
void nmod_poly_one(nmod_poly_t res)
{
    nmod_poly_fit_length(res, 1);
    res->length = (res->mod.n != UWORD(1));
    res->coeffs[0] = 1;
}

void nmod_poly_set_trunc(nmod_poly_t res, const nmod_poly_t poly, slong len);

NMOD_POLY_INLINE
void nmod_poly_truncate(nmod_poly_t poly, slong len)
{
    nmod_poly_set_trunc(poly, poly, len);
}

void _nmod_poly_reverse(nn_ptr output, nn_srcptr input, slong len, slong m);
void nmod_poly_reverse(nmod_poly_t output, const nmod_poly_t input, slong m);

/* Comparison  ***************************************************************/

int nmod_poly_equal(const nmod_poly_t a, const nmod_poly_t b);

int nmod_poly_equal_nmod(const nmod_poly_t poly, ulong cst);
int nmod_poly_equal_ui(const nmod_poly_t poly, ulong cst);

int nmod_poly_equal_trunc(const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

NMOD_POLY_INLINE
int nmod_poly_is_zero(const nmod_poly_t poly)
{
    return (poly->length == 0);
}

NMOD_POLY_INLINE
int nmod_poly_is_one(const nmod_poly_t poly)
{
    return (poly->mod.n == UWORD(1)) || (poly->length == 1 && poly->coeffs[0] == 1);
}

/* bogus for non-prime modulus */
NMOD_POLY_INLINE
int nmod_poly_is_unit(const nmod_poly_t poly)
{
    return (poly->mod.n == UWORD(1)) || ((poly->length == 1) && poly->coeffs[0] != 0);
}

NMOD_POLY_INLINE
int nmod_poly_is_gen(const nmod_poly_t poly)
{
    return (poly->mod.n == UWORD(1)) ||
           (poly->length == 2 && poly->coeffs[0] == 0 && poly->coeffs[1] == 1);
}

NMOD_POLY_INLINE
int nmod_poly_is_monic(const nmod_poly_t poly)
{
    return (poly->length && poly->coeffs[(poly->length - 1)] == 1);
}

/* Randomisation  ************************************************************/

void nmod_poly_randtest(nmod_poly_t poly, flint_rand_t state, slong len);

NMOD_POLY_INLINE
void nmod_poly_randtest_not_zero(nmod_poly_t poly, flint_rand_t state, slong len)
{
    do {
        nmod_poly_randtest(poly, state, len);
    } while (nmod_poly_is_zero(poly));
}

void nmod_poly_randtest_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);

void nmod_poly_randtest_monic(nmod_poly_t poly, flint_rand_t state, slong len);
void nmod_poly_randtest_monic_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);
void nmod_poly_randtest_monic_primitive(nmod_poly_t poly, flint_rand_t state, slong len);

void nmod_poly_randtest_trinomial(nmod_poly_t poly, flint_rand_t state, slong len);
int nmod_poly_randtest_trinomial_irreducible(nmod_poly_t poly, flint_rand_t state, slong len, slong max_attempts);

void nmod_poly_randtest_pentomial(nmod_poly_t poly, flint_rand_t state, slong len);
int nmod_poly_randtest_pentomial_irreducible(nmod_poly_t poly, flint_rand_t state, slong len, slong max_attempts);

void nmod_poly_randtest_sparse_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);

/* Getting and setting coefficients  *****************************************/

NMOD_POLY_INLINE
ulong nmod_poly_get_coeff_ui(const nmod_poly_t poly, slong j)
{
    return (j >= poly->length) ? 0 : poly->coeffs[j];
}

void nmod_poly_set_coeff_ui(nmod_poly_t poly, slong j, ulong c);

/* Input and output  *********************************************************/

char * nmod_poly_get_str(const nmod_poly_t poly);

char * nmod_poly_get_str_pretty(const nmod_poly_t poly, const char * x);

int nmod_poly_set_str(nmod_poly_t poly, const char * s);

#ifdef FLINT_HAVE_FILE
int nmod_poly_fprint(FILE * f, const nmod_poly_t poly);
int nmod_poly_fprint_pretty(FILE * f, const nmod_poly_t a, const char * x);

int nmod_poly_fread(FILE * f, nmod_poly_t poly);
#endif

int nmod_poly_print(const nmod_poly_t a);
int nmod_poly_print_pretty(const nmod_poly_t a, const char * x);

int nmod_poly_read(nmod_poly_t poly);

/* Shifting  *****************************************************************/

void _nmod_poly_shift_left(nn_ptr res, nn_srcptr poly, slong len, slong k);
void nmod_poly_shift_left(nmod_poly_t res, const nmod_poly_t poly, slong k);

void _nmod_poly_shift_right(nn_ptr res, nn_srcptr poly, slong len, slong k);
void nmod_poly_shift_right(nmod_poly_t res, const nmod_poly_t poly, slong k);

/* Addition and subtraction  *************************************************/

void _nmod_poly_add(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod);
void nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_sub(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod);
void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

void nmod_poly_add_series(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);
void nmod_poly_sub_series(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

void nmod_poly_add_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c);
void nmod_poly_sub_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c);

void nmod_poly_neg(nmod_poly_t res, const nmod_poly_t poly1);

/* Scalar multiplication and division  ***************************************/

void nmod_poly_scalar_mul_nmod(nmod_poly_t res, const nmod_poly_t poly, ulong c);

void nmod_poly_scalar_addmul_nmod(nmod_poly_t res, const nmod_poly_t poly, ulong c);

void _nmod_poly_make_monic(nn_ptr output, nn_srcptr input, slong len, nmod_t mod);
void nmod_poly_make_monic(nmod_poly_t output, const nmod_poly_t input);

/* Bit packing and unpacking aand reduction  **********************************/

void _nmod_poly_KS2_pack1(nn_ptr res, nn_srcptr op, slong n, slong s, ulong b, ulong k, slong r);
void _nmod_poly_KS2_pack(nn_ptr res, nn_srcptr op, slong n, slong s, ulong b, ulong k, slong r);

void _nmod_poly_KS2_unpack1(nn_ptr res, nn_srcptr op, slong n, ulong b, ulong k);
void _nmod_poly_KS2_unpack2(nn_ptr res, nn_srcptr op, slong n, ulong b, ulong k);
void _nmod_poly_KS2_unpack3(nn_ptr res, nn_srcptr op, slong n, ulong b, ulong k);
void _nmod_poly_KS2_unpack(nn_ptr res, nn_srcptr op, slong n, ulong b, ulong k);

void _nmod_poly_KS2_reduce(nn_ptr res, slong s, nn_srcptr op, slong n, ulong w, nmod_t mod);

void _nmod_poly_KS2_recover_reduce1(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod);
void _nmod_poly_KS2_recover_reduce2(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod);
void _nmod_poly_KS2_recover_reduce2b(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong FLINT_UNUSED(b), nmod_t mod);
void _nmod_poly_KS2_recover_reduce3(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod);
void _nmod_poly_KS2_recover_reduce(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod);

void _nmod_poly_bit_pack(nn_ptr res, nn_srcptr poly, slong len, flint_bitcnt_t bits);
void nmod_poly_bit_pack(fmpz_t f, const nmod_poly_t poly, flint_bitcnt_t bit_size);

void _nmod_poly_bit_unpack(nn_ptr res, slong len, nn_srcptr mpn, flint_bitcnt_t bits, nmod_t mod);
void nmod_poly_bit_unpack(nmod_poly_t poly, const fmpz_t f, flint_bitcnt_t bit_size);

/* Multiplication  ***********************************************************/

void _nmod_poly_mul_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod);
void nmod_poly_mul_classical(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong trunc, nmod_t mod);
void nmod_poly_mullow_classical(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong trunc);

void _nmod_poly_mulhigh_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong start, nmod_t mod);
void nmod_poly_mulhigh_classical(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong start);

void _nmod_poly_mul_KS(nn_ptr out, nn_srcptr in1, slong len1, nn_srcptr in2, slong len2, flint_bitcnt_t bits, nmod_t mod);
void nmod_poly_mul_KS(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, flint_bitcnt_t bits);

void _nmod_poly_mul_KS2(nn_ptr res, nn_srcptr op1, slong n1, nn_srcptr op2, slong n2, nmod_t mod);
void nmod_poly_mul_KS2(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mul_KS4(nn_ptr res, nn_srcptr op1, slong n1, nn_srcptr op2, slong n2, nmod_t mod);
void nmod_poly_mul_KS4(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow_KS(nn_ptr out, nn_srcptr in1, slong len1, nn_srcptr in2, slong len2, flint_bitcnt_t bits, slong n, nmod_t mod);
void nmod_poly_mullow_KS(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, flint_bitcnt_t bits, slong n);

void _nmod_poly_mul(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod);
void nmod_poly_mul(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong trunc, nmod_t mod);
void nmod_poly_mullow(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong trunc);

void _nmod_poly_mulhigh(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong n, nmod_t mod);
void nmod_poly_mulhigh(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

void _nmod_poly_mulmod(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nn_srcptr f, slong lenf, nmod_t mod);
void nmod_poly_mulmod(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t f);

void _nmod_poly_mulmod_preinv(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod);
void nmod_poly_mulmod_preinv(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t f, const nmod_poly_t finv);

int _nmod_poly_invmod(ulong *A, const ulong *B, slong lenB, const ulong *P, slong lenP, const nmod_t mod);
int nmod_poly_invmod(nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t P);

/* Powering  *****************************************************************/

void _nmod_poly_pow_binexp(nn_ptr res, nn_srcptr poly, slong len, ulong e, nmod_t mod);
void nmod_poly_pow_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e);

void _nmod_poly_pow(nn_ptr res, nn_srcptr poly, slong len, ulong e, nmod_t mod);
void nmod_poly_pow(nmod_poly_t res, const nmod_poly_t poly, ulong e);

void _nmod_poly_pow_trunc_binexp(nn_ptr res, nn_srcptr poly, ulong e, slong trunc, nmod_t mod);
void nmod_poly_pow_trunc_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e, slong trunc);

void _nmod_poly_pow_trunc(nn_ptr res, nn_srcptr poly, ulong e, slong trunc, nmod_t mod);
void nmod_poly_pow_trunc(nmod_poly_t res, const nmod_poly_t poly, ulong e, slong trunc);

void _nmod_poly_powmod_ui_binexp(nn_ptr res, nn_srcptr poly, ulong e, nn_srcptr f, slong lenf, nmod_t mod);
void nmod_poly_powmod_ui_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e, const nmod_poly_t f);

void _nmod_poly_powmod_fmpz_binexp(nn_ptr res, nn_srcptr poly, fmpz_t e, nn_srcptr f, slong lenf, nmod_t mod);
void nmod_poly_powmod_fmpz_binexp(nmod_poly_t res, const nmod_poly_t poly, fmpz_t e, const nmod_poly_t f);

void _nmod_poly_powmod_ui_binexp_preinv (nn_ptr res, nn_srcptr poly, ulong e, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod);
void nmod_poly_powmod_ui_binexp_preinv(nmod_poly_t res, const nmod_poly_t poly, ulong e, const nmod_poly_t f, const nmod_poly_t finv);

void _nmod_poly_powmod_fmpz_binexp_preinv (nn_ptr res, nn_srcptr poly, fmpz_t e, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod);
void nmod_poly_powmod_fmpz_binexp_preinv(nmod_poly_t res, const nmod_poly_t poly, fmpz_t e, const nmod_poly_t f, const nmod_poly_t finv);

void _nmod_poly_powmod_x_ui_preinv (nn_ptr res, ulong e, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod);
void nmod_poly_powmod_x_ui_preinv(nmod_poly_t res, ulong e, const nmod_poly_t f, const nmod_poly_t finv);

void _nmod_poly_powmod_x_fmpz_preinv (nn_ptr res, fmpz_t e, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod);
void nmod_poly_powmod_x_fmpz_preinv(nmod_poly_t res, fmpz_t e, const nmod_poly_t f,                             const nmod_poly_t finv);

void _nmod_poly_powers_mod_preinv_naive(nn_ptr * res, nn_srcptr f,
		 slong flen, slong n, nn_srcptr g, slong glen, nn_srcptr ginv,
		                              slong ginvlen, const nmod_t mod);

void nmod_poly_powers_mod_naive(nmod_poly_struct * res,
                            const nmod_poly_t f, slong n, const nmod_poly_t g);

void _nmod_poly_powers_mod_preinv_threaded_pool(nn_ptr * res,
	       nn_srcptr f, slong flen, slong n, nn_srcptr g, slong glen,
			    nn_srcptr ginv, slong ginvlen, const nmod_t mod,
                              thread_pool_handle * threads, slong num_threads);

void
_nmod_poly_powers_mod_preinv_threaded(nn_ptr * res, nn_srcptr f,
		                 slong flen, slong n, nn_srcptr g, slong glen,
                              nn_srcptr ginv, slong ginvlen, const nmod_t mod);

void nmod_poly_powers_mod_bsgs(nmod_poly_struct * res,
                            const nmod_poly_t f, slong n, const nmod_poly_t g);

/* Division  *****************************************************************/

void _nmod_poly_divrem_basecase_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong A_len, nn_srcptr B, slong B_len, ulong invB, nmod_t mod);

void _nmod_poly_divrem_basecase(nn_ptr Q, nn_ptr R, nn_srcptr A, slong A_len, nn_srcptr B, slong B_len, nmod_t mod);
void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_divrem(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_div(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_div(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_rem(nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_rem(nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_divexact(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_divexact(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_inv_series_basecase(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, nmod_t mod);
void nmod_poly_inv_series_basecase(nmod_poly_t Qinv, const nmod_poly_t Q, slong n);

void _nmod_poly_inv_series_newton(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, nmod_t mod);
void nmod_poly_inv_series_newton(nmod_poly_t Qinv, const nmod_poly_t Q, slong n);

void _nmod_poly_inv_series(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, nmod_t mod);
void nmod_poly_inv_series(nmod_poly_t Qinv, const nmod_poly_t Q, slong n);

void _nmod_poly_div_series_basecase(nn_ptr Q, nn_srcptr A, slong Alen, nn_srcptr B, slong Blen, slong n, nmod_t mod);
void nmod_poly_div_series_basecase(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B, slong n);

void _nmod_poly_div_series(nn_ptr Q, nn_srcptr A, slong Alen, nn_srcptr B, slong Blen, slong n, nmod_t mod);
void nmod_poly_div_series(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B, slong n);

void _nmod_poly_div_newton_n_preinv(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr FLINT_UNUSED(B), slong lenB, nn_srcptr Binv, slong lenBinv, nmod_t mod);
void nmod_poly_div_newton_n_preinv (nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t Binv);

void _nmod_poly_divrem_newton_n_preinv (nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr Binv, slong lenBinv, nmod_t mod);
void nmod_poly_divrem_newton_n_preinv(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t Binv);

ulong _nmod_poly_div_root(nn_ptr Q, nn_srcptr A, slong len, ulong c, nmod_t mod);

ulong nmod_poly_div_root(nmod_poly_t Q, const nmod_poly_t A, ulong c);

/* Divisibility testing  *****************************************************/

int _nmod_poly_divides_classical(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
int nmod_poly_divides_classical(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B);

int _nmod_poly_divides(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
int nmod_poly_divides(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B);

ulong nmod_poly_remove(nmod_poly_t f, const nmod_poly_t p);

/* Derivative  ***************************************************************/

void _nmod_poly_derivative(nn_ptr x_prime, nn_srcptr x, slong len, nmod_t mod);
void nmod_poly_derivative(nmod_poly_t x_prime, const nmod_poly_t x);

void _nmod_poly_integral(nn_ptr x_int, nn_srcptr x, slong len, nmod_t mod);
void nmod_poly_integral(nmod_poly_t x_int, const nmod_poly_t x);

/* Evaluation  ***************************************************************/

ulong _nmod_poly_evaluate_nmod(nn_srcptr poly, slong len, ulong c, nmod_t mod);
ulong _nmod_poly_evaluate_nmod_precomp(nn_srcptr poly, slong len, ulong c, ulong c_precomp, nmod_t mod);
ulong _nmod_poly_evaluate_nmod_precomp_lazy(nn_srcptr poly, slong len, ulong c, ulong c_precomp, nmod_t mod);
ulong nmod_poly_evaluate_nmod(const nmod_poly_t poly, ulong c);

void _nmod_poly_evaluate_nmod_vec(nn_ptr ys, nn_srcptr coeffs, slong len, nn_srcptr xs, slong n, nmod_t mod);
void nmod_poly_evaluate_nmod_vec(nn_ptr ys, const nmod_poly_t poly, nn_srcptr xs, slong n);

void _nmod_poly_evaluate_nmod_vec_iter(nn_ptr ys, nn_srcptr coeffs, slong len, nn_srcptr xs, slong n, nmod_t mod);
void nmod_poly_evaluate_nmod_vec_iter(nn_ptr ys, const nmod_poly_t poly, nn_srcptr xs, slong n);

void _nmod_poly_evaluate_nmod_vec_fast_precomp(nn_ptr vs, nn_srcptr poly, slong plen, const nn_ptr * tree, slong len, nmod_t mod);

void _nmod_poly_evaluate_nmod_vec_fast(nn_ptr ys, nn_srcptr coeffs, slong len, nn_srcptr xs, slong n, nmod_t mod);
void nmod_poly_evaluate_nmod_vec_fast(nn_ptr ys, const nmod_poly_t poly, nn_srcptr xs, slong n);

void nmod_mat_one_addmul(nmod_mat_t dest, const nmod_mat_t mat, ulong c);

void nmod_poly_evaluate_mat_horner(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c);

void nmod_poly_evaluate_mat_paterson_stockmeyer(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c);

NMOD_POLY_INLINE
void nmod_poly_evaluate_mat(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c)
{
    if (poly->length < 5 || c->r * poly->length < 425)
    {
        nmod_poly_evaluate_mat_horner(dest, poly, c);
    }
    else
    {
        nmod_poly_evaluate_mat_paterson_stockmeyer(dest, poly, c);
    }
}

/* Subproduct tree  **********************************************************/

nn_ptr * _nmod_poly_tree_alloc(slong len);

void _nmod_poly_tree_free(nn_ptr * tree, slong len);

void _nmod_poly_tree_build(nn_ptr * tree, nn_srcptr roots, slong len, nmod_t mod);

/* Interpolation  ************************************************************/

void _nmod_poly_interpolate_nmod_vec_newton(nn_ptr poly, nn_srcptr xs, nn_srcptr ys, slong n, nmod_t mod);
void nmod_poly_interpolate_nmod_vec_newton(nmod_poly_t poly, nn_srcptr xs, nn_srcptr ys, slong n);

void _nmod_poly_interpolate_nmod_vec_barycentric(nn_ptr poly, nn_srcptr xs, nn_srcptr ys, slong n, nmod_t mod);
void nmod_poly_interpolate_nmod_vec_barycentric(nmod_poly_t poly, nn_srcptr xs, nn_srcptr ys, slong n);

void _nmod_poly_interpolate_nmod_vec(nn_ptr poly, nn_srcptr xs, nn_srcptr ys, slong n, nmod_t mod);
void nmod_poly_interpolate_nmod_vec(nmod_poly_t poly, nn_srcptr xs, nn_srcptr ys, slong n);

void _nmod_poly_interpolate_nmod_vec_fast(nn_ptr poly, nn_srcptr xs, nn_srcptr ys, slong len, nmod_t mod);
void nmod_poly_interpolate_nmod_vec_fast(nmod_poly_t poly, nn_srcptr xs, nn_srcptr ys, slong n);

void _nmod_poly_interpolate_nmod_vec_fast_precomp(nn_ptr poly, nn_srcptr ys,
    const nn_ptr * tree, nn_srcptr weights, slong len, nmod_t mod);

void _nmod_poly_interpolation_weights(nn_ptr w, const nn_ptr * tree, slong len, nmod_t mod);

/* Composition  **************************************************************/

void _nmod_poly_compose_horner(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod);
void nmod_poly_compose_horner(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_compose(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod);
void nmod_poly_compose(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2);

/* obsolete implementation */
#define _nmod_poly_compose_divconquer _nmod_poly_compose
#define nmod_poly_compose_divconquer nmod_poly_compose

/* Taylor shift  *************************************************************/

void _nmod_poly_taylor_shift_horner(nn_ptr poly, ulong c, slong len, nmod_t mod);
void nmod_poly_taylor_shift_horner(nmod_poly_t g, const nmod_poly_t f, ulong c);

void _nmod_poly_taylor_shift_convolution(nn_ptr poly, ulong c, slong len, nmod_t mod);
void nmod_poly_taylor_shift_convolution(nmod_poly_t g, const nmod_poly_t f, ulong c);

void _nmod_poly_taylor_shift(nn_ptr poly, ulong c, slong len, nmod_t mod);
void nmod_poly_taylor_shift(nmod_poly_t g, const nmod_poly_t f, ulong c);

/* Modular composition  ******************************************************/

void _nmod_poly_compose_mod_brent_kung(nn_ptr res, nn_srcptr f, slong lenf, nn_srcptr g, nn_srcptr h, slong lenh, nmod_t mod);
void nmod_poly_compose_mod_brent_kung(nmod_poly_t res, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t h);

void _nmod_poly_reduce_matrix_mod_poly(nmod_mat_t A, const nmod_mat_t B, const nmod_poly_t f);

void _nmod_poly_precompute_matrix(nmod_mat_t A, nn_srcptr poly1, nn_srcptr poly2,
               slong len2, nn_srcptr poly2inv, slong len2inv, nmod_t mod);

void _nmod_poly_precompute_matrix_worker(void * arg_ptr);

void nmod_poly_precompute_matrix(nmod_mat_t A, const nmod_poly_t poly1,
                          const nmod_poly_t poly2, const nmod_poly_t poly2inv);

void _nmod_poly_compose_mod_brent_kung_precomp_preinv(nn_ptr res, nn_srcptr poly1,
                            slong len1, const nmod_mat_t A, nn_srcptr poly3,
                            slong len3, nn_srcptr poly3inv, slong len3inv,
                            nmod_t mod);

void _nmod_poly_compose_mod_brent_kung_precomp_preinv_worker(void * arg_ptr);

void nmod_poly_compose_mod_brent_kung_precomp_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_mat_t A,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv);

void _nmod_poly_compose_mod_brent_kung_preinv(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, nn_srcptr poly3, slong len3,
                            nn_srcptr poly3inv, slong len3inv, nmod_t mod);

void nmod_poly_compose_mod_brent_kung_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv);

void
_nmod_poly_compose_mod_brent_kung_vec_preinv(nmod_poly_struct * res,
        const nmod_poly_struct * polys, slong FLINT_UNUSED(lenpolys), slong l,
        nn_srcptr g, slong glen, nn_srcptr poly, slong len,
        nn_srcptr polyinv, slong leninv, nmod_t mod);

void nmod_poly_compose_mod_brent_kung_vec_preinv(nmod_poly_struct * res,
                    const nmod_poly_struct * polys, slong len1, slong n,
                    const nmod_poly_t g, const nmod_poly_t poly,
		    const nmod_poly_t polyinv);

void _nmod_poly_compose_mod_brent_kung_vec_preinv_worker(void * arg_ptr);

void
nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(nmod_poly_struct * res,
           const nmod_poly_struct * polys, slong len1, slong n,
                          const nmod_poly_t g, const nmod_poly_t poly,
                     const nmod_poly_t polyinv, thread_pool_handle * threads,
                                                                slong num_threads);

void _nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(
        nmod_poly_struct * res,
        const nmod_poly_struct * polys,
        slong FLINT_UNUSED(lenpolys), slong l,
        nn_srcptr g, slong glen,
        nn_srcptr poly, slong len,
        nn_srcptr polyinv, slong leninv,
        nmod_t mod,
        thread_pool_handle * threads,
        slong num_threads);

void nmod_poly_compose_mod_brent_kung_vec_preinv_threaded(
        nmod_poly_struct * res,
        const nmod_poly_struct * polys,
        slong len1, slong n,
        const nmod_poly_t g, const nmod_poly_t poly,
        const nmod_poly_t polyinv);

void _nmod_poly_compose_mod_horner(nn_ptr res, nn_srcptr f, slong lenf, nn_srcptr g, nn_srcptr h, slong lenh, nmod_t mod);
void nmod_poly_compose_mod_horner(nmod_poly_t res, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t h);

void _nmod_poly_compose_mod(nn_ptr res, nn_srcptr f, slong lenf, nn_srcptr g, nn_srcptr h, slong lenh, nmod_t mod);
void nmod_poly_compose_mod(nmod_poly_t res, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t h);

/* Power series composition and reversion ************************************/

void _nmod_poly_compose_series(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong n, nmod_t mod);
void nmod_poly_compose_series(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

void _nmod_poly_revert_series(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, nmod_t mod);
void nmod_poly_revert_series(nmod_poly_t Qinv, const nmod_poly_t Q, slong n);

/* norms *********************************************************************/

NMOD_POLY_INLINE slong _nmod_poly_hamming_weight(nn_srcptr a, slong len)
{
    slong i, sum = 0;
    for (i = 0; i < len; i++)
        sum += !(a[i] == 0);
    return sum;
}

NMOD_POLY_INLINE slong nmod_poly_hamming_weight(const nmod_poly_t A)
{
    return _nmod_poly_hamming_weight(A->coeffs, A->length);
}

/* Greatest common divisor  **************************************************/

slong _nmod_poly_gcd_euclidean(nn_ptr G, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_gcd_euclidean(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_hgcd_recursive(nn_ptr *M, slong *lenM, nn_ptr A, slong *lenA, nn_ptr B, slong *lenB, nn_srcptr a, slong lena, nn_srcptr b, slong lenb, nn_ptr P, nmod_t mod, int flag, nmod_poly_res_t res);

slong nmod_poly_hgcd_ref(nmod_poly_t m11, nmod_poly_t m12, nmod_poly_t m21, nmod_poly_t m22, nmod_poly_t A, nmod_poly_t B, const nmod_poly_t a, const nmod_poly_t b);

slong _nmod_poly_hgcd(nn_ptr *M, slong *lenM, nn_ptr A, slong *lenA, nn_ptr B, slong *lenB, nn_srcptr a, slong lena, nn_srcptr b, slong lenb, nmod_t mod);
slong nmod_poly_hgcd(nmod_poly_t m11, nmod_poly_t m12, nmod_poly_t m21, nmod_poly_t m22, nmod_poly_t A, nmod_poly_t B, const nmod_poly_t a, const nmod_poly_t b);

slong _nmod_poly_gcd_hgcd(nn_ptr G, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_gcd_hgcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_gcd(nn_ptr G, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_gcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_xgcd_euclidean(nn_ptr res, nn_ptr s, nn_ptr t, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod);
void nmod_poly_xgcd_euclidean(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T, const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_xgcd_hgcd(nn_ptr G, nn_ptr S, nn_ptr T, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_xgcd_hgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T, const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_xgcd(nn_ptr G, nn_ptr S, nn_ptr T, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
void nmod_poly_xgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T, const nmod_poly_t A, const nmod_poly_t B);

ulong _nmod_poly_resultant_euclidean(nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod);
ulong nmod_poly_resultant_euclidean(const nmod_poly_t f, const nmod_poly_t g);

ulong _nmod_poly_resultant_hgcd(nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
ulong nmod_poly_resultant_hgcd(const nmod_poly_t A, const nmod_poly_t B);

ulong _nmod_poly_resultant(nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod);
ulong nmod_poly_resultant(const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_gcdinv(ulong *G, ulong *S, const ulong *A, slong lenA, const ulong *B, slong lenB, const nmod_t mod);
void nmod_poly_gcdinv(nmod_poly_t G, nmod_poly_t S, const nmod_poly_t A, const nmod_poly_t B);

/* Discriminant **************************************************************/

ulong _nmod_poly_discriminant(nn_srcptr poly, slong len, nmod_t mod);
ulong nmod_poly_discriminant(const nmod_poly_t f);

/* Square roots **************************************************************/

void _nmod_poly_invsqrt_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void nmod_poly_invsqrt_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_sqrt_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void nmod_poly_sqrt_series(nmod_poly_t g, const nmod_poly_t h, slong n);

int _nmod_poly_sqrt(nn_ptr s, nn_srcptr p, slong len, nmod_t mod);
int nmod_poly_sqrt(nmod_poly_t b, const nmod_poly_t a);

/* Power sums ****************************************************************/

void _nmod_poly_power_sums_naive(nn_ptr res, nn_srcptr poly, slong len, slong n, nmod_t mod);
void nmod_poly_power_sums_naive(nmod_poly_t res, const nmod_poly_t poly, slong n);

void _nmod_poly_power_sums_schoenhage(nn_ptr res, nn_srcptr poly, slong len, slong n, nmod_t mod);
void nmod_poly_power_sums_schoenhage(nmod_poly_t res, const nmod_poly_t poly, slong n);

void _nmod_poly_power_sums(nn_ptr res, nn_srcptr poly, slong len, slong n, nmod_t mod);
void nmod_poly_power_sums(nmod_poly_t res, const nmod_poly_t poly, slong n);

void _nmod_poly_power_sums_to_poly_naive(nn_ptr res, nn_srcptr poly, slong len, nmod_t mod);
void nmod_poly_power_sums_to_poly_naive(nmod_poly_t res, const nmod_poly_t Q);

void _nmod_poly_power_sums_to_poly_schoenhage(nn_ptr res, nn_srcptr poly, slong len, nmod_t mod);
void nmod_poly_power_sums_to_poly_schoenhage(nmod_poly_t res, const nmod_poly_t Q);

void _nmod_poly_power_sums_to_poly(nn_ptr res, nn_srcptr poly, slong len, nmod_t mod);
void nmod_poly_power_sums_to_poly(nmod_poly_t res, const nmod_poly_t Q);

/* Transcendental functions **************************************************/

void _nmod_poly_atan_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void nmod_poly_atan_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_tan_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void nmod_poly_tan_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_asin_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void nmod_poly_asin_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_sin_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod);
void nmod_poly_sin_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_cos_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod);
void nmod_poly_cos_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_asinh_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void nmod_poly_asinh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_atanh_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void nmod_poly_atanh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_sinh_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod);
void nmod_poly_sinh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_cosh_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod);
void nmod_poly_cosh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_tanh_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod);
void nmod_poly_tanh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_log_series(nn_ptr res, nn_srcptr f, slong flen, slong n, nmod_t mod);
void nmod_poly_log_series(nmod_poly_t res, const nmod_poly_t f, slong n);

void  _nmod_poly_exp_expinv_series(nn_ptr f, nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void _nmod_poly_exp_series(nn_ptr f, nn_srcptr h, slong hlen, slong n, nmod_t mod);
void nmod_poly_exp_series(nmod_poly_t f, const nmod_poly_t h, slong n);

/* Special polynomials *******************************************************/

int _nmod_poly_conway(nn_ptr op, ulong prime, slong deg);

ulong _nmod_poly_conway_rand(slong * degree, flint_rand_t state, int type);

/* Products  *****************************************************************/

void nmod_poly_product_roots_nmod_vec(nmod_poly_t poly, nn_srcptr xs, slong n);

void _nmod_poly_product_roots_nmod_vec(nn_ptr poly, nn_srcptr xs, slong n, nmod_t mod);

void _nmod_poly_split_rabin(nmod_poly_t a, nmod_poly_t b, const nmod_poly_t f, nmod_poly_t t, nmod_poly_t t2, flint_rand_t randstate);

int nmod_poly_find_distinct_nonzero_roots(ulong * roots, const nmod_poly_t P);

/* CRT ***********************************************************************/

/* instructions do A = B + I*(C - B) mod M */
typedef struct
{
    slong a_idx; /* index of A */
    slong b_idx; /* index of B */
    slong c_idx; /* index of C */
    nmod_poly_t idem;     /* I */
    nmod_poly_t modulus;  /* M */
} _nmod_poly_multi_crt_prog_instr;

typedef struct
{
    _nmod_poly_multi_crt_prog_instr * prog; /* straight line program */
    slong length; /* length of prog */
    slong alloc;  /* alloc of prog */
    slong localsize; /* length of outputs required in nmod_poly_multi_crt_run */
    slong temp1loc; /* index of temporary used in run */
    slong temp2loc; /* index of another temporary used in run */
    int good;   /* the moduli are good for CRT, essentially relatively prime */
} nmod_poly_multi_crt_struct;

typedef nmod_poly_multi_crt_struct nmod_poly_multi_crt_t[1];

void nmod_poly_multi_crt_init(nmod_poly_multi_crt_t CRT);
void nmod_poly_multi_crt_clear(nmod_poly_multi_crt_t CRT);

int nmod_poly_multi_crt_precompute(nmod_poly_multi_crt_t CRT, const nmod_poly_struct * moduli, slong len);
int nmod_poly_multi_crt_precompute_p(nmod_poly_multi_crt_t CRT, const nmod_poly_struct * const * moduli, slong len);

void nmod_poly_multi_crt_precomp(nmod_poly_t output, const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * values);
void nmod_poly_multi_crt_precomp_p(nmod_poly_t output, const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * const * values);

int nmod_poly_multi_crt(nmod_poly_t output, const nmod_poly_struct * moduli, const nmod_poly_struct * values, slong len);

NMOD_POLY_INLINE
slong _nmod_poly_multi_crt_local_size(const nmod_poly_multi_crt_t CRT)
{
    return CRT->localsize;
}

void _nmod_poly_multi_crt_run(nmod_poly_struct * outputs, const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * inputs);
void _nmod_poly_multi_crt_run_p(nmod_poly_struct * outputs, const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * const * inputs);

/* Inflation and deflation ***************************************************/

slong nmod_poly_deflation(const nmod_poly_t input);

void nmod_poly_deflate(nmod_poly_t result, const nmod_poly_t input, slong deflation);
void nmod_poly_inflate(nmod_poly_t result, const nmod_poly_t input, slong inflation);

/* Characteristic polynomial and minimal polynomial */
/* FIXME: These should be moved to nmod_mat.h. */

void _nmod_mat_charpoly_berkowitz(nn_ptr p, const nmod_mat_t M, nmod_t mod);
void nmod_mat_charpoly_berkowitz(nmod_poly_t p, const nmod_mat_t M);
void nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t M);
void nmod_mat_charpoly(nmod_poly_t p, const nmod_mat_t M);

void nmod_mat_minpoly_with_gens(nmod_poly_t p, const nmod_mat_t X, ulong * P);

void nmod_mat_minpoly(nmod_poly_t p, const nmod_mat_t M);

/* Berlekamp-Massey Algorithm - see nmod_poly/berlekamp_massey.c for more info ************/

typedef struct
{
    slong npoints;
    nmod_poly_t R0, R1;
    nmod_poly_t V0, V1;
    nmod_poly_t qt, rt;
    nmod_poly_t points;
} nmod_berlekamp_massey_struct;

typedef nmod_berlekamp_massey_struct nmod_berlekamp_massey_t[1];

void nmod_berlekamp_massey_init(nmod_berlekamp_massey_t B, ulong p);
void nmod_berlekamp_massey_clear(nmod_berlekamp_massey_t B);

void nmod_berlekamp_massey_start_over(nmod_berlekamp_massey_t B);

void nmod_berlekamp_massey_set_prime(nmod_berlekamp_massey_t B, ulong p);

void nmod_berlekamp_massey_print(const nmod_berlekamp_massey_t B);

void nmod_berlekamp_massey_add_point(nmod_berlekamp_massey_t B, ulong a);
void nmod_berlekamp_massey_add_points(nmod_berlekamp_massey_t B, const ulong * a, slong count);
void nmod_berlekamp_massey_add_zeros(nmod_berlekamp_massey_t B, slong count);

int nmod_berlekamp_massey_reduce(nmod_berlekamp_massey_t B);

NMOD_POLY_INLINE const ulong * nmod_berlekamp_massey_points(const nmod_berlekamp_massey_t B)
{
    return B->points->coeffs;
}

NMOD_POLY_INLINE slong nmod_berlekamp_massey_point_count(const nmod_berlekamp_massey_t B)
{
    return B->points->length;
}

NMOD_POLY_INLINE const nmod_poly_struct * nmod_berlekamp_massey_V_poly(const nmod_berlekamp_massey_t B)
{
    return B->V1;
}

NMOD_POLY_INLINE const nmod_poly_struct * nmod_berlekamp_massey_R_poly(const nmod_berlekamp_massey_t B)
{
    return B->R1;
}

#ifdef __cplusplus
}
#endif

#endif
