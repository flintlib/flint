/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2010, 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_H
#define NMOD_POLY_H

#include "nmod_poly_mini.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NMOD_DIVREM_DIVCONQUER_CUTOFF  300
#define NMOD_DIV_DIVCONQUER_CUTOFF     300 /* Must be <= NMOD_DIVREM_DIVCONQUER_CUTOFF */

#define NMOD_POLY_HGCD_CUTOFF  100      /* HGCD: Basecase -> Recursion      */
#define NMOD_POLY_GCD_CUTOFF  340       /* GCD:  Euclidean -> HGCD          */
#define NMOD_POLY_SMALL_GCD_CUTOFF 200  /* GCD (small n): Euclidean -> HGCD */

NMOD_POLY_INLINE
slong NMOD_DIVREM_BC_ITCH(slong lenA, slong lenB, nmod_t mod)
{
    const flint_bitcnt_t bits = 
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);
    
    if (bits <= FLINT_BITS)
        return lenA;
    else if (bits <= 2 * FLINT_BITS)
        return 2*(lenA + lenB - 1);
    else
        return 3*(lenA + lenB - 1);
}

NMOD_POLY_INLINE
slong NMOD_DIV_BC_ITCH(slong lenA, slong lenB, nmod_t mod)
{
    const flint_bitcnt_t bits = 
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);
    
    if (bits <= FLINT_BITS)
        return lenA - lenB + 1;
    else if (bits <= 2 * FLINT_BITS)
        return 2*lenA;
    else
        return 3*lenA;
}

NMOD_POLY_INLINE
slong NMOD_DIVREM_DC_ITCH(slong lenB, nmod_t mod)
{
    slong i = 0;
    
    while (lenB > NMOD_DIVREM_DIVCONQUER_CUTOFF + i)
    {
        lenB = (lenB + 1)/2;
        i++;
    }
    if (lenB > NMOD_DIVREM_DIVCONQUER_CUTOFF)
        lenB = NMOD_DIVREM_DIVCONQUER_CUTOFF;

    return NMOD_DIVREM_BC_ITCH(2*lenB - 1, lenB, mod) + 2*lenB - 1;
}

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

/* zn_poly helper functions  ************************************************

Copyright (C) 2007, 2008 David Harvey

*/

#ifdef __GMP_H__
NMOD_POLY_INLINE
int signed_mpn_sub_n(ulong_ptr res, ulong_srcptr op1, ulong_srcptr op2, slong n)
{
   if (mpn_cmp(op1, op2, n) >= 0)
   {
      mpn_sub_n(res, op1, op2, n);
      return 0;
   }
   else
   {
      mpn_sub_n(res, op2, op1, n);
      return 1;
   }
}
#endif

/* Polynomial parameters  ****************************************************/

NMOD_POLY_INLINE
flint_bitcnt_t nmod_poly_max_bits(const nmod_poly_t poly)
{
    return _nmod_vec_max_bits(poly->coeffs, poly->length);
}

NMOD_POLY_INLINE
ulong_ptr nmod_poly_lead(const nmod_poly_t poly)
{
    if (poly->length)
        return poly->coeffs + (poly->length - 1);
    else
        return NULL;
}

/* Assignment and basic manipulation  ****************************************/

FLINT_DLL void nmod_poly_set_trunc(nmod_poly_t res,
                                              const nmod_poly_t poly, slong n);

FLINT_DLL void _nmod_poly_reverse(ulong_ptr output,
                                          ulong_srcptr input, slong len, slong m);

FLINT_DLL void nmod_poly_reverse(nmod_poly_t output, 
                                             const nmod_poly_t input, slong m);

/* Comparison  ***************************************************************/

FLINT_DLL int nmod_poly_equal_trunc(const nmod_poly_t poly1, 
                                             const nmod_poly_t poly2, slong n);

/* bogus for non-prime modulus */
NMOD_POLY_INLINE
int nmod_poly_is_unit(const nmod_poly_t poly)
{
    return (poly->length == 1) && poly->coeffs[0] != 0;
}

NMOD_POLY_INLINE
int nmod_poly_is_gen(const nmod_poly_t poly)
{
    return (poly->mod.n == 0) ||
           (poly->length == 2 && poly->coeffs[0] == 0 && poly->coeffs[1] == 1);
}

/* Randomisation  ************************************************************/

FLINT_DLL void nmod_poly_randtest(nmod_poly_t poly, flint_rand_t state, slong len);

NMOD_POLY_INLINE
void nmod_poly_randtest_not_zero(nmod_poly_t poly, flint_rand_t state, slong len)
{
    do {
        nmod_poly_randtest(poly, state, len);
    } while (nmod_poly_is_zero(poly));
}
      
FLINT_DLL void nmod_poly_randtest_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);

FLINT_DLL void nmod_poly_randtest_monic(nmod_poly_t poly, flint_rand_t state, slong len);

FLINT_DLL void nmod_poly_randtest_monic_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);

FLINT_DLL void nmod_poly_randtest_monic_primitive(nmod_poly_t poly, flint_rand_t state, slong len);

FLINT_DLL void nmod_poly_randtest_trinomial(nmod_poly_t poly, flint_rand_t state, slong len);

FLINT_DLL int nmod_poly_randtest_trinomial_irreducible(nmod_poly_t poly, flint_rand_t state,
                                         slong len, slong max_attempts);

FLINT_DLL void nmod_poly_randtest_pentomial(nmod_poly_t poly, flint_rand_t state, slong len);

FLINT_DLL int nmod_poly_randtest_pentomial_irreducible(nmod_poly_t poly, flint_rand_t state,
                                         slong len, slong max_attempts);

FLINT_DLL void nmod_poly_randtest_sparse_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);

/* Getting and setting coefficients  *****************************************/

NMOD_POLY_INLINE
ulong nmod_poly_get_coeff_ui(const nmod_poly_t poly, slong j)
{
    return (j >= poly->length) ? 0 : poly->coeffs[j];
}

/* Input and output  *********************************************************/

FLINT_DLL char * nmod_poly_get_str(const nmod_poly_t poly);

FLINT_DLL char * nmod_poly_get_str_pretty(const nmod_poly_t poly, const char * x);

FLINT_DLL int nmod_poly_set_str(nmod_poly_t poly, const char * s);

#if defined (H_STDIO)               \
  || defined (_H_STDIO)             \
  || defined (_STDIO_H)             \
  || defined (_STDIO_H_)            \
  || defined (__STDIO_H)            \
  || defined (__STDIO_H__)          \
  || defined (_STDIO_INCLUDED)      \
  || defined (__dj_include_stdio_h_)\
  || defined (_FILE_DEFINED)        \
  || defined (__STDIO__)            \
  || defined (_MSL_STDIO_H)         \
  || defined (_STDIO_H_INCLUDED)    \
  || defined (_ISO_STDIO_ISO_H)     \
  || defined (__STDIO_LOADED)       \
  || defined (_STDIO)

#include "flint-impl.h"

FLINT_DLL int nmod_poly_fread(FILE * f, nmod_poly_t poly);

NMOD_POLY_INLINE
int nmod_poly_fprint(FILE * f, const nmod_poly_t poly)
{
    char *s;
    int r;

    s = nmod_poly_get_str(poly);
    r = fputs(s, f);
    flint_free(s);

    return (r < 0) ? r : 1;
}

FLINT_DLL int nmod_poly_fprint_pretty(FILE * f, const nmod_poly_t a, const char * x);

NMOD_POLY_INLINE
int nmod_poly_print(const nmod_poly_t a)
{
    size_t r;
    slong i;

    r = printf(WORD_FMT "d " WORD_FMT "u", a->length, a->mod.n);

    if (a->length == 0)
        return r;
    else
        if (r > 0)
            r = printf(" ");

    for (i = 0; (r > 0) && (i < a->length); i++)
        r = printf(" " WORD_FMT "u", a->coeffs[i]);

    return (int) r;
}

NMOD_POLY_INLINE
int nmod_poly_print_pretty(const nmod_poly_t a, const char * x)
{
    return nmod_poly_fprint_pretty(stdout, a, x);
}

NMOD_POLY_INLINE
int nmod_poly_read(nmod_poly_t poly)
{
    return nmod_poly_fread(stdin, poly);
}
#endif

/* Addition and subtraction  *************************************************/

FLINT_DLL void nmod_poly_add_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c);

FLINT_DLL void nmod_poly_add_series(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

FLINT_DLL void nmod_poly_sub_series(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

FLINT_DLL void nmod_poly_sub_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c);

FLINT_DLL void nmod_poly_neg(nmod_poly_t res, const nmod_poly_t poly1);

/* Scalar multiplication and division  ***************************************/

FLINT_DLL void nmod_poly_scalar_mul_nmod(nmod_poly_t res, 
                                         const nmod_poly_t poly1, ulong c);

FLINT_DLL void nmod_poly_scalar_addmul_nmod(nmod_poly_t A, const nmod_poly_t B,
                                                                      ulong x);

/* Bit packing and unpacking aand reduction  **********************************/

FLINT_DLL void _nmod_poly_KS2_pack1(ulong_ptr res, ulong_srcptr op, slong n, slong s,
                                                    ulong b, ulong k, slong r);

FLINT_DLL void _nmod_poly_KS2_pack(ulong_ptr res, ulong_srcptr op, slong n, slong s,
                                                    ulong b, ulong k, slong r);

FLINT_DLL void _nmod_poly_KS2_unpack1(ulong_ptr res, ulong_srcptr op, slong n, ulong b,
                                                                      ulong k);

FLINT_DLL void _nmod_poly_KS2_unpack2(ulong_ptr res, ulong_srcptr op, slong n, ulong b,
                                                                      ulong k);

FLINT_DLL void _nmod_poly_KS2_unpack3(ulong_ptr res, ulong_srcptr op, slong n, ulong b,
                                                                      ulong k);

FLINT_DLL void _nmod_poly_KS2_unpack(ulong_ptr res, ulong_srcptr op, slong n, ulong b,
                                                                      ulong k);

FLINT_DLL void _nmod_poly_KS2_reduce(ulong_ptr res, slong s, ulong_srcptr op, 
                                                 slong n, ulong w, nmod_t mod);

FLINT_DLL void _nmod_poly_KS2_recover_reduce1(ulong_ptr res, slong s, ulong_srcptr op1,
                                  ulong_srcptr op2, slong n, ulong b, nmod_t mod);

FLINT_DLL void _nmod_poly_KS2_recover_reduce2(ulong_ptr res, slong s, ulong_srcptr op1,
                                  ulong_srcptr op2, slong n, ulong b, nmod_t mod);

FLINT_DLL void _nmod_poly_KS2_recover_reduce2b(ulong_ptr res, slong s, ulong_srcptr op1,
                                  ulong_srcptr op2, slong n, ulong b, nmod_t mod);

FLINT_DLL void _nmod_poly_KS2_recover_reduce3(ulong_ptr res, slong s, ulong_srcptr op1,
                                  ulong_srcptr op2, slong n, ulong b, nmod_t mod);

FLINT_DLL void _nmod_poly_KS2_recover_reduce(ulong_ptr res, slong s, ulong_srcptr op1,
                                  ulong_srcptr op2, slong n, ulong b, nmod_t mod);

FLINT_DLL void _nmod_poly_bit_pack(ulong_ptr res, ulong_srcptr poly, 
                                                  slong len, flint_bitcnt_t bits);

FLINT_DLL void _nmod_poly_bit_unpack(ulong_ptr res, slong len, 
                                  ulong_srcptr mpn, flint_bitcnt_t bits, nmod_t mod);

FLINT_DLL void nmod_poly_bit_pack(fmpz_t f, const nmod_poly_t poly,
                   flint_bitcnt_t bit_size);

FLINT_DLL void nmod_poly_bit_unpack(nmod_poly_t poly, const fmpz_t f, flint_bitcnt_t bit_size);

/* Multiplication  ***********************************************************/

FLINT_DLL void _nmod_poly_mul_classical(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                                       ulong_srcptr poly2, slong len2, nmod_t mod);

FLINT_DLL void nmod_poly_mul_classical(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

FLINT_DLL void _nmod_poly_mullow_classical(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                           ulong_srcptr poly2, slong len2, slong trunc, nmod_t mod);

FLINT_DLL void nmod_poly_mullow_classical(nmod_poly_t res, 
                  const nmod_poly_t poly1, const nmod_poly_t poly2, slong trunc);

FLINT_DLL void _nmod_poly_mulhigh_classical(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                           ulong_srcptr poly2, slong len2, slong start, nmod_t mod);

FLINT_DLL void nmod_poly_mulhigh_classical(nmod_poly_t res, 
                  const nmod_poly_t poly1, const nmod_poly_t poly2, slong start);

FLINT_DLL void _nmod_poly_mul_KS(ulong_ptr out, ulong_srcptr in1, slong len1, 
                        ulong_srcptr in2, slong len2, flint_bitcnt_t bits, nmod_t mod);

FLINT_DLL void nmod_poly_mul_KS(nmod_poly_t res, 
             const nmod_poly_t poly1, const nmod_poly_t poly2, flint_bitcnt_t bits);

FLINT_DLL void _nmod_poly_mul_KS2(ulong_ptr res, ulong_srcptr op1, slong n1,
                                            ulong_srcptr op2, slong n2, nmod_t mod);

FLINT_DLL void nmod_poly_mul_KS2(nmod_poly_t res,
                               const nmod_poly_t poly1, const nmod_poly_t poly2);

FLINT_DLL void _nmod_poly_mul_KS4(ulong_ptr res, ulong_srcptr op1, slong n1,
                                            ulong_srcptr op2, slong n2, nmod_t mod);

FLINT_DLL void nmod_poly_mul_KS4(nmod_poly_t res,
                               const nmod_poly_t poly1, const nmod_poly_t poly2);

FLINT_DLL void _nmod_poly_mullow_KS(ulong_ptr out, ulong_srcptr in1, slong len1,
               ulong_srcptr in2, slong len2, flint_bitcnt_t bits, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_mullow_KS(nmod_poly_t res, const nmod_poly_t poly1, 
                             const nmod_poly_t poly2, flint_bitcnt_t bits, slong n);

FLINT_DLL void _nmod_poly_mullow(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                           ulong_srcptr poly2, slong len2, slong trunc, nmod_t mod);

FLINT_DLL void nmod_poly_mullow(nmod_poly_t res, const nmod_poly_t poly1, 
                                          const nmod_poly_t poly2, slong trunc);

FLINT_DLL void _nmod_poly_mulhigh(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                               ulong_srcptr poly2, slong len2, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_mulhigh(nmod_poly_t res, const nmod_poly_t poly1, 
                                              const nmod_poly_t poly2, slong n);

FLINT_DLL void _nmod_poly_mulmod(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                             ulong_srcptr poly2, slong len2, ulong_srcptr f,
                            slong lenf, nmod_t mod);

FLINT_DLL void nmod_poly_mulmod(nmod_poly_t res,
    const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t f);

FLINT_DLL void _nmod_poly_mulmod_preinv(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                              ulong_srcptr poly2, slong len2, ulong_srcptr f,
                              slong lenf, ulong_srcptr finv, slong lenfinv, nmod_t mod);

FLINT_DLL void nmod_poly_mulmod_preinv(nmod_poly_t res, const nmod_poly_t poly1,
                        const nmod_poly_t poly2, const nmod_poly_t f,
                        const nmod_poly_t finv);

FLINT_DLL int _nmod_poly_invmod(ulong_ptr A,
                      ulong_srcptr B, slong lenB,
                      ulong_srcptr P, slong lenP, const nmod_t mod);

FLINT_DLL int nmod_poly_invmod(nmod_poly_t A, 
                     const nmod_poly_t B, const nmod_poly_t P);

/* Powering  *****************************************************************/

FLINT_DLL void _nmod_poly_pow_binexp(ulong_ptr res, 
                              ulong_srcptr poly, slong len, ulong e, nmod_t mod);

FLINT_DLL void nmod_poly_pow_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e);

FLINT_DLL void _nmod_poly_pow(ulong_ptr res, ulong_srcptr poly, slong len, ulong e, nmod_t mod);

FLINT_DLL void nmod_poly_pow(nmod_poly_t res, const nmod_poly_t poly, ulong e);

FLINT_DLL void _nmod_poly_pow_trunc_binexp(ulong_ptr res, ulong_srcptr poly, 
                                              ulong e, slong trunc, nmod_t mod);

FLINT_DLL void nmod_poly_pow_trunc_binexp(nmod_poly_t res, 
                                  const nmod_poly_t poly, ulong e, slong trunc);

FLINT_DLL void _nmod_poly_pow_trunc(ulong_ptr res, ulong_srcptr poly, 
                                              ulong e, slong trunc, nmod_t mod);

FLINT_DLL void nmod_poly_pow_trunc(nmod_poly_t res, 
                                  const nmod_poly_t poly, ulong e, slong trunc);

FLINT_DLL void nmod_poly_powmod_ui_binexp(nmod_poly_t res, 
                           const nmod_poly_t poly, ulong e,
                           const nmod_poly_t f);

FLINT_DLL void _nmod_poly_powmod_ui_binexp(ulong_ptr res, ulong_srcptr poly, 
                                ulong e, ulong_srcptr f, slong lenf, nmod_t mod);

FLINT_DLL void nmod_poly_powmod_fmpz_binexp(nmod_poly_t res,
                           const nmod_poly_t poly, fmpz_t e,
                           const nmod_poly_t f);

FLINT_DLL void _nmod_poly_powmod_fmpz_binexp(ulong_ptr res, ulong_srcptr poly,
                                fmpz_t e, ulong_srcptr f, slong lenf, nmod_t mod);

#ifdef __GMP_H__
FLINT_DLL void _nmod_poly_powmod_mpz_binexp(ulong_ptr res, ulong_srcptr poly, 
                                mpz_srcptr e, ulong_srcptr f,
                                slong lenf, nmod_t mod);

FLINT_DLL void nmod_poly_powmod_mpz_binexp(nmod_poly_t res, 
                           const nmod_poly_t poly, mpz_srcptr e,
                           const nmod_poly_t f);
#endif

FLINT_DLL void _nmod_poly_powmod_ui_binexp_preinv (ulong_ptr res, ulong_srcptr poly,
                                    ulong e, ulong_srcptr f, slong lenf,
                                    ulong_srcptr finv, slong lenfinv, nmod_t mod);

FLINT_DLL void nmod_poly_powmod_ui_binexp_preinv(nmod_poly_t res, 
                           const nmod_poly_t poly, ulong e,
                           const nmod_poly_t f, const nmod_poly_t finv);

FLINT_DLL void _nmod_poly_powmod_fmpz_binexp_preinv (ulong_ptr res, ulong_srcptr poly,
                                    fmpz_t e, ulong_srcptr f, slong lenf,
                                    ulong_srcptr finv, slong lenfinv, nmod_t mod);

FLINT_DLL void nmod_poly_powmod_fmpz_binexp_preinv(nmod_poly_t res,
                           const nmod_poly_t poly, fmpz_t e,
                           const nmod_poly_t f, const nmod_poly_t finv);

FLINT_DLL void _nmod_poly_powmod_x_ui_preinv (ulong_ptr res, ulong e, ulong_srcptr f, slong lenf,
                               ulong_srcptr finv, slong lenfinv, nmod_t mod);

FLINT_DLL void nmod_poly_powmod_x_ui_preinv(nmod_poly_t res, ulong e, const nmod_poly_t f,
                             const nmod_poly_t finv);

#ifdef __GMP_H__
FLINT_DLL void _nmod_poly_powmod_mpz_binexp_preinv (ulong_ptr res, ulong_srcptr poly,
                                    mpz_srcptr e, ulong_srcptr f, slong lenf,
                                    ulong_srcptr finv, slong lenfinv, nmod_t mod);

FLINT_DLL void nmod_poly_powmod_mpz_binexp_preinv(nmod_poly_t res,
                           const nmod_poly_t poly, mpz_srcptr e,
                           const nmod_poly_t f, const nmod_poly_t finv);
#endif

FLINT_DLL void _nmod_poly_powmod_x_fmpz_preinv (ulong_ptr res, fmpz_t e, ulong_srcptr f, slong lenf,
                               ulong_srcptr finv, slong lenfinv, nmod_t mod);

FLINT_DLL void nmod_poly_powmod_x_fmpz_preinv(nmod_poly_t res, fmpz_t e, const nmod_poly_t f,                             const nmod_poly_t finv);

FLINT_DLL void _nmod_poly_powers_mod_preinv_naive(ulong_ptr * res, ulong_srcptr f,
		 slong flen, slong n, ulong_srcptr g, slong glen, ulong_srcptr ginv,
		                              slong ginvlen, const nmod_t mod);

FLINT_DLL void nmod_poly_powers_mod_naive(nmod_poly_struct * res,
                            const nmod_poly_t f, slong n, const nmod_poly_t g);

#ifdef THREAD_POOL_H
FLINT_DLL void _nmod_poly_powers_mod_preinv_threaded_pool(ulong_ptr * res,
	       ulong_srcptr f, slong flen, slong n, ulong_srcptr g, slong glen,
			    ulong_srcptr ginv, slong ginvlen, const nmod_t mod,
                              thread_pool_handle * threads, slong num_threads);
#endif

FLINT_DLL void
_nmod_poly_powers_mod_preinv_threaded(ulong_ptr * res, ulong_srcptr f,
		                 slong flen, slong n, ulong_srcptr g, slong glen,
                              ulong_srcptr ginv, slong ginvlen, const nmod_t mod);

FLINT_DLL void nmod_poly_powers_mod_bsgs(nmod_poly_struct * res,
                            const nmod_poly_t f, slong n, const nmod_poly_t g);

/* Division  *****************************************************************/

FLINT_DLL void _nmod_poly_divrem_basecase(ulong_ptr Q, ulong_ptr R, ulong_ptr W,
               ulong_srcptr A, slong A_len, ulong_srcptr B, slong B_len, nmod_t mod);

FLINT_DLL void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, 
                                     const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL void _nmod_poly_divrem_divconquer_recursive(ulong_ptr Q, ulong_ptr BQ, 
         ulong_ptr W, ulong_ptr V, ulong_srcptr A, ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void _nmod_poly_divrem_divconquer(ulong_ptr Q, ulong_ptr R, 
                 ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void nmod_poly_divrem_divconquer(nmod_poly_t Q, nmod_poly_t R,
                                     const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL void _nmod_poly_divrem_q0(ulong_ptr Q, ulong_ptr R, 
                             ulong_srcptr A, ulong_srcptr B, slong lenA, nmod_t mod);

FLINT_DLL void _nmod_poly_divrem_q1(ulong_ptr Q, ulong_ptr R, 
                 ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void _nmod_poly_div_basecase(ulong_ptr Q, ulong_ptr W,
               ulong_srcptr A, slong A_len, ulong_srcptr B, slong B_len, nmod_t mod);

FLINT_DLL void nmod_poly_div_basecase(nmod_poly_t Q, const nmod_poly_t A,
                                                          const nmod_poly_t B);

FLINT_DLL void _nmod_poly_div_divconquer_recursive(ulong_ptr Q, 
         ulong_ptr W, ulong_ptr V, ulong_srcptr A, ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void _nmod_poly_div_divconquer(ulong_ptr Q, ulong_srcptr A, slong lenA, 
                                          ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void nmod_poly_div_divconquer(nmod_poly_t Q,
                                     const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL void _nmod_poly_rem_basecase(ulong_ptr R, 
                               ulong_ptr W, ulong_srcptr A, slong lenA, 
                                          ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void nmod_poly_rem_basecase(nmod_poly_t R,
                                     const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL void _nmod_poly_rem_q1(ulong_ptr R, 
                 ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void _nmod_poly_inv_series_basecase(ulong_ptr Qinv, 
                                 ulong_srcptr Q, slong Qlen, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_inv_series_basecase(nmod_poly_t Qinv, 
                                                 const nmod_poly_t Q, slong n);

FLINT_DLL void _nmod_poly_inv_series_newton(ulong_ptr Qinv, 
                                 ulong_srcptr Q, slong Qlen, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_inv_series_newton(nmod_poly_t Qinv, 
                                                 const nmod_poly_t Q, slong n);

NMOD_POLY_INLINE
void _nmod_poly_inv_series(ulong_ptr Qinv, ulong_srcptr Q, 
                                               slong Qlen, slong n, nmod_t mod)
{
    _nmod_poly_inv_series_newton(Qinv, Q, Qlen, n, mod);
}

NMOD_POLY_INLINE
void nmod_poly_inv_series(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)
{
    nmod_poly_inv_series_newton(Qinv, Q, n);
}

FLINT_DLL void _nmod_poly_div_series_basecase(ulong_ptr Q, ulong_srcptr A, 
                     slong Alen, ulong_srcptr B, slong Blen, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_div_series_basecase(nmod_poly_t Q,
                            const nmod_poly_t A, const nmod_poly_t B, slong n);

FLINT_DLL void _nmod_poly_div_series(ulong_ptr Q, ulong_srcptr A, slong Alen,
                                 ulong_srcptr B, slong Blen, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_div_series(nmod_poly_t Q, const nmod_poly_t A, 
                                                 const nmod_poly_t B, slong n);

FLINT_DLL void _nmod_poly_div_newton(ulong_ptr Q, ulong_srcptr A, slong Alen, 
                                          ulong_srcptr B, slong Blen, nmod_t mod);

FLINT_DLL void nmod_poly_div_newton(nmod_poly_t Q, const nmod_poly_t A,
                                                          const nmod_poly_t B);

FLINT_DLL void _nmod_poly_divrem_newton(ulong_ptr Q, ulong_ptr R, 
                 ulong_srcptr A, slong Alen, ulong_srcptr B, slong Blen, nmod_t mod);

FLINT_DLL void nmod_poly_divrem_newton(nmod_poly_t Q, nmod_poly_t R, 
                                     const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL void _nmod_poly_div_newton_n_preinv (ulong_ptr Q, 
                   ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, 
                                    ulong_srcptr Binv, slong lenBinv, nmod_t mod);

FLINT_DLL void nmod_poly_div_newton_n_preinv (nmod_poly_t Q,
             const nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t Binv);

FLINT_DLL void _nmod_poly_divrem_newton_n_preinv (ulong_ptr Q, 
         ulong_ptr R, ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB,
                                    ulong_srcptr Binv, slong lenBinv, nmod_t mod);

FLINT_DLL void nmod_poly_divrem_newton_n_preinv(nmod_poly_t Q, nmod_poly_t R,
             const nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t Binv);

FLINT_DLL ulong _nmod_poly_div_root(ulong_ptr Q, 
                              ulong_srcptr A, slong len, ulong c, nmod_t mod);

FLINT_DLL ulong nmod_poly_div_root(nmod_poly_t Q,
                                             const nmod_poly_t A, ulong c);

/* Divisibility testing  *****************************************************/

FLINT_DLL int _nmod_poly_divides_classical(ulong_ptr Q, ulong_srcptr A, slong lenA,
                                          ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL int nmod_poly_divides_classical(nmod_poly_t Q,
		                     const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL int _nmod_poly_divides(ulong_ptr Q, ulong_srcptr A, slong lenA,
                                          ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL int nmod_poly_divides(nmod_poly_t Q,
		                     const nmod_poly_t A, const nmod_poly_t B);

/* Derivative  ***************************************************************/

FLINT_DLL void _nmod_poly_derivative(ulong_ptr x_prime, 
                                           ulong_srcptr x, slong len, nmod_t mod);

FLINT_DLL void nmod_poly_derivative(nmod_poly_t x_prime, const nmod_poly_t x);

FLINT_DLL void _nmod_poly_integral(ulong_ptr 
                                    x_int, ulong_srcptr x, slong len, nmod_t mod);

FLINT_DLL void nmod_poly_integral(nmod_poly_t x_int, const nmod_poly_t x);

/* Evaluation  ***************************************************************/

FLINT_DLL void _nmod_poly_evaluate_fmpz(fmpz_t rop,
                             ulong_srcptr poly, const slong len, const fmpz_t c);

FLINT_DLL void nmod_poly_evaluate_fmpz(fmpz_t rop,
                                       const nmod_poly_t poly, const fmpz_t c);

FLINT_DLL void _nmod_poly_evaluate_nmod_vec(ulong_ptr ys,
               ulong_srcptr coeffs, slong len, ulong_srcptr xs, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_evaluate_nmod_vec(ulong_ptr ys, 
                                const nmod_poly_t poly, ulong_srcptr xs, slong n);

FLINT_DLL void _nmod_poly_evaluate_nmod_vec_iter(ulong_ptr ys,
               ulong_srcptr coeffs, slong len, ulong_srcptr xs, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_evaluate_nmod_vec_iter(ulong_ptr ys,
                                const nmod_poly_t poly, ulong_srcptr xs, slong n);

FLINT_DLL void _nmod_poly_evaluate_nmod_vec_fast_precomp(ulong_ptr vs, 
       ulong_srcptr poly, slong plen, const ulong_ptr * tree, slong len, nmod_t mod);

FLINT_DLL void _nmod_poly_evaluate_nmod_vec_fast(ulong_ptr ys, 
               ulong_srcptr coeffs, slong len, ulong_srcptr xs, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_evaluate_nmod_vec_fast(ulong_ptr ys,
                                const nmod_poly_t poly, ulong_srcptr xs, slong n);

FLINT_DLL void nmod_mat_one_addmul(nmod_mat_t dest, 
                                            const nmod_mat_t mat, ulong c);

FLINT_DLL void nmod_poly_evaluate_mat_horner(nmod_mat_t dest,
                                   const nmod_poly_t poly, const nmod_mat_t c);

FLINT_DLL void nmod_poly_evaluate_mat_paterson_stockmeyer(nmod_mat_t dest,
                                   const nmod_poly_t poly, const nmod_mat_t c);

NMOD_POLY_INLINE
void nmod_poly_evaluate_mat(nmod_mat_t dest,
	const nmod_poly_t poly, const nmod_mat_t c)
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

FLINT_DLL ulong_ptr * _nmod_poly_tree_alloc(slong len);

FLINT_DLL void _nmod_poly_tree_free(ulong_ptr * tree, slong len);

FLINT_DLL void _nmod_poly_tree_build(ulong_ptr * tree, ulong_srcptr roots,
    slong len, nmod_t mod);

/* Interpolation  ************************************************************/

FLINT_DLL void _nmod_poly_interpolate_nmod_vec_newton(ulong_ptr poly, ulong_srcptr xs,
                        ulong_srcptr ys, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_interpolate_nmod_vec_newton(nmod_poly_t poly,
                        ulong_srcptr xs, ulong_srcptr ys, slong n);

FLINT_DLL void _nmod_poly_interpolate_nmod_vec_barycentric(ulong_ptr poly, ulong_srcptr xs,
                        ulong_srcptr ys, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_interpolate_nmod_vec_barycentric(nmod_poly_t poly,
                        ulong_srcptr xs, ulong_srcptr ys, slong n);

FLINT_DLL void _nmod_poly_interpolate_nmod_vec(ulong_ptr poly, ulong_srcptr xs,
                        ulong_srcptr ys, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_interpolate_nmod_vec(nmod_poly_t poly,
                        ulong_srcptr xs, ulong_srcptr ys, slong n);

FLINT_DLL void nmod_poly_interpolate_nmod_vec_fast(nmod_poly_t poly,
                                    ulong_srcptr xs, ulong_srcptr ys, slong n);

FLINT_DLL void _nmod_poly_interpolate_nmod_vec_fast(ulong_ptr poly,
                            ulong_srcptr xs, ulong_srcptr ys, slong len, nmod_t mod);

FLINT_DLL void _nmod_poly_interpolate_nmod_vec_fast_precomp(ulong_ptr poly, ulong_srcptr ys,
    const ulong_ptr * tree, ulong_srcptr weights, slong len, nmod_t mod);

FLINT_DLL void _nmod_poly_interpolation_weights(ulong_ptr w, const ulong_ptr * tree,
    slong len, nmod_t mod);

/* Composition  **************************************************************/

FLINT_DLL void _nmod_poly_compose_horner(ulong_ptr res, ulong_srcptr poly1, 
                            slong len1, ulong_srcptr poly2, slong len2, nmod_t mod);

FLINT_DLL void nmod_poly_compose_horner(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

FLINT_DLL void _nmod_poly_compose_divconquer(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                                       ulong_srcptr poly2, slong len2, nmod_t mod);
FLINT_DLL void nmod_poly_compose_divconquer(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

FLINT_DLL void _nmod_poly_compose(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                                       ulong_srcptr poly2, slong len2, nmod_t mod);

FLINT_DLL void nmod_poly_compose(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

/* Taylor shift  *************************************************************/

FLINT_DLL void _nmod_poly_taylor_shift_horner(ulong_ptr poly, ulong c,
    slong len, nmod_t mod);

FLINT_DLL void nmod_poly_taylor_shift_horner(nmod_poly_t g,
    const nmod_poly_t f, ulong c);

FLINT_DLL void _nmod_poly_taylor_shift_convolution(ulong_ptr poly, ulong c,
    slong len, nmod_t mod);

FLINT_DLL void nmod_poly_taylor_shift_convolution(nmod_poly_t g,
    const nmod_poly_t f, ulong c);

/* Modular composition  ******************************************************/

FLINT_DLL void _nmod_poly_compose_mod_brent_kung(ulong_ptr res, ulong_srcptr f, slong lenf,
                            ulong_srcptr g, ulong_srcptr h, slong lenh, nmod_t mod);

FLINT_DLL void nmod_poly_compose_mod_brent_kung(nmod_poly_t res, 
                    const nmod_poly_t f, const nmod_poly_t g,
                    const nmod_poly_t h);

FLINT_DLL void _nmod_poly_reduce_matrix_mod_poly(nmod_mat_t A, const nmod_mat_t B,
                          const nmod_poly_t f);

FLINT_DLL void _nmod_poly_precompute_matrix(nmod_mat_t A, ulong_srcptr poly1, ulong_srcptr poly2,
               slong len2, ulong_srcptr poly2inv, slong len2inv, nmod_t mod);

FLINT_DLL void _nmod_poly_precompute_matrix_worker(void * arg_ptr);

FLINT_DLL void nmod_poly_precompute_matrix(nmod_mat_t A, const nmod_poly_t poly1,
                          const nmod_poly_t poly2, const nmod_poly_t poly2inv);

FLINT_DLL void _nmod_poly_compose_mod_brent_kung_precomp_preinv(ulong_ptr res, ulong_srcptr poly1,
                            slong len1, const nmod_mat_t A, ulong_srcptr poly3,
                            slong len3, ulong_srcptr poly3inv, slong len3inv,
                            nmod_t mod);

FLINT_DLL void _nmod_poly_compose_mod_brent_kung_precomp_preinv_worker(void * arg_ptr);

FLINT_DLL void nmod_poly_compose_mod_brent_kung_precomp_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_mat_t A,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv);

FLINT_DLL void _nmod_poly_compose_mod_brent_kung_preinv(ulong_ptr res, ulong_srcptr poly1, slong len1,
                            ulong_srcptr poly2, ulong_srcptr poly3, slong len3,
                            ulong_srcptr poly3inv, slong len3inv, nmod_t mod);

FLINT_DLL void nmod_poly_compose_mod_brent_kung_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv);

FLINT_DLL void _nmod_poly_compose_mod_brent_kung_vec_preinv(nmod_poly_struct * res,
                 const nmod_poly_struct * polys, slong len1, slong l,
                 ulong_srcptr g, slong glen, ulong_srcptr poly, slong len,
		 ulong_srcptr polyinv,slong leninv, nmod_t mod);

FLINT_DLL void nmod_poly_compose_mod_brent_kung_vec_preinv(nmod_poly_struct * res,
                    const nmod_poly_struct * polys, slong len1, slong n,
                    const nmod_poly_t g, const nmod_poly_t poly,
		    const nmod_poly_t polyinv);

FLINT_DLL void _nmod_poly_compose_mod_brent_kung_vec_preinv_worker(void * arg_ptr);

#ifdef THREAD_POOL_H
FLINT_DLL void
nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(nmod_poly_struct * res,
           const nmod_poly_struct * polys, slong len1, slong n,
                          const nmod_poly_t g, const nmod_poly_t poly,
                     const nmod_poly_t polyinv, thread_pool_handle * threads,
                                                                slong num_threads);

FLINT_DLL void _nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(
                 nmod_poly_struct * res, const nmod_poly_struct * polys,
                 slong lenpolys, slong l, ulong_srcptr g, slong glen,
                 ulong_srcptr poly, slong len, ulong_srcptr polyinv, slong leninv,
                 nmod_t mod, thread_pool_handle * threads, slong num_threads);
#endif

FLINT_DLL void nmod_poly_compose_mod_brent_kung_vec_preinv_threaded(nmod_poly_struct * res,
                                            const nmod_poly_struct * polys,
                                            slong len1, slong n,
                                            const nmod_poly_t g,
					    const nmod_poly_t poly,
                                            const nmod_poly_t polyinv);

FLINT_DLL void _nmod_poly_compose_mod_horner(ulong_ptr res,
    ulong_srcptr f, slong lenf, ulong_srcptr g, ulong_srcptr h, slong lenh, nmod_t mod);

FLINT_DLL void nmod_poly_compose_mod_horner(nmod_poly_t res, 
                    const nmod_poly_t f, const nmod_poly_t g,
                    const nmod_poly_t h);

FLINT_DLL void _nmod_poly_compose_mod(ulong_ptr res, ulong_srcptr f, slong lenf, 
                            ulong_srcptr g,
                            ulong_srcptr h, slong lenh, nmod_t mod);

FLINT_DLL void nmod_poly_compose_mod(nmod_poly_t res, 
                    const nmod_poly_t f, const nmod_poly_t g,
                    const nmod_poly_t h);

/* Power series composition and reversion ************************************/

FLINT_DLL void _nmod_poly_compose_series_horner(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                            ulong_srcptr poly2, slong len2, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_compose_series_horner(nmod_poly_t res, 
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

FLINT_DLL void _nmod_poly_compose_series_brent_kung(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                            ulong_srcptr poly2, slong len2, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_compose_series_brent_kung(nmod_poly_t res, 
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

FLINT_DLL void _nmod_poly_compose_series(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                            ulong_srcptr poly2, slong len2, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_compose_series(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

FLINT_DLL void _nmod_poly_revert_series_lagrange(ulong_ptr Qinv, ulong_srcptr Q, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_revert_series_lagrange(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n);

FLINT_DLL void _nmod_poly_revert_series_lagrange_fast(ulong_ptr Qinv, ulong_srcptr Q,
    slong n, nmod_t mod);

FLINT_DLL void nmod_poly_revert_series_lagrange_fast(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n);

FLINT_DLL void _nmod_poly_revert_series_newton(ulong_ptr Qinv, ulong_srcptr Q, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_revert_series_newton(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n);

FLINT_DLL void _nmod_poly_revert_series(ulong_ptr Qinv, ulong_srcptr Q, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_revert_series(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n);

FLINT_DLL void _nmod_poly_compose_series_divconquer(ulong_ptr res, ulong_srcptr poly1, slong len1, 
                                                 ulong_srcptr poly2, slong len2, 
                                                 slong N, nmod_t mod);

FLINT_DLL void nmod_poly_compose_series_divconquer(nmod_poly_t res, 
    const nmod_poly_t poly1, const nmod_poly_t poly2, slong N);

/* norms *********************************************************************/

NMOD_POLY_INLINE slong _nmod_poly_hamming_weight(ulong_srcptr a, slong len)
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

FLINT_DLL slong _nmod_poly_gcd_euclidean(ulong_ptr G, 
                   ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void nmod_poly_gcd_euclidean(nmod_poly_t G, 
                                     const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL slong _nmod_poly_hgcd_recursive(ulong_ptr *M, slong *lenM, 
    ulong_ptr A, slong *lenA, ulong_ptr B, slong *lenB, 
    ulong_srcptr a, slong lena, ulong_srcptr b, slong lenb, 
    ulong_ptr P, nmod_t mod, int flag, nmod_poly_res_t res);

FLINT_DLL slong _nmod_poly_hgcd(ulong_ptr *M, slong *lenM, 
                     ulong_ptr A, slong *lenA, ulong_ptr B, slong *lenB, 
                     ulong_srcptr a, slong lena, ulong_srcptr b, slong lenb, 
                     nmod_t mod);

FLINT_DLL slong nmod_poly_hgcd_ref(
        nmod_poly_t m11, nmod_poly_t m12, nmod_poly_t m21, nmod_poly_t m22,
        nmod_poly_t A, nmod_poly_t B, const nmod_poly_t a, const nmod_poly_t b);

FLINT_DLL slong nmod_poly_hgcd(
        nmod_poly_t m11, nmod_poly_t m12, nmod_poly_t m21, nmod_poly_t m22,
        nmod_poly_t A, nmod_poly_t B, const nmod_poly_t a, const nmod_poly_t b);

FLINT_DLL slong _nmod_poly_gcd_hgcd(ulong_ptr G, ulong_srcptr A, slong lenA, 
                                   ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void nmod_poly_gcd_hgcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL slong _nmod_poly_gcd(ulong_ptr G, ulong_srcptr A, slong lenA, 
                              ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void nmod_poly_gcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL slong _nmod_poly_xgcd_euclidean(ulong_ptr res, ulong_ptr s, ulong_ptr t, 
           ulong_srcptr poly1, slong len1, ulong_srcptr poly2, slong len2, nmod_t mod);

FLINT_DLL void nmod_poly_xgcd_euclidean(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                                     const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL slong _nmod_poly_xgcd_hgcd(ulong_ptr G, ulong_ptr S, ulong_ptr T, 
                          ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, 
                          nmod_t mod);

FLINT_DLL void nmod_poly_xgcd_hgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                         const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL slong _nmod_poly_xgcd(ulong_ptr G, ulong_ptr S, ulong_ptr T, 
            ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL void nmod_poly_xgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                                   const nmod_poly_t A, const nmod_poly_t B);

FLINT_DLL ulong _nmod_poly_resultant_euclidean(ulong_srcptr poly1, slong len1, 
                               ulong_srcptr poly2, slong len2, nmod_t mod);

FLINT_DLL ulong nmod_poly_resultant_euclidean(const nmod_poly_t f, const nmod_poly_t g);

FLINT_DLL ulong _nmod_poly_resultant_hgcd(ulong_srcptr A, slong lenA, 
                         ulong_srcptr B, slong lenB, nmod_t mod);

FLINT_DLL ulong nmod_poly_resultant_hgcd(const nmod_poly_t A, const nmod_poly_t B);

NMOD_POLY_INLINE
ulong _nmod_poly_resultant(ulong_srcptr poly1, slong len1, 
                     ulong_srcptr poly2, slong len2, nmod_t mod)
{
    const slong cutoff = FLINT_BIT_COUNT(mod.n) <= 8 ? 
                        NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;

    if (len1 < cutoff)
        return _nmod_poly_resultant_euclidean(poly1, len1, poly2, len2, mod);
    else
        return _nmod_poly_resultant_hgcd(poly1, len1, poly2, len2, mod);
}

NMOD_POLY_INLINE
ulong nmod_poly_resultant(const nmod_poly_t f, const nmod_poly_t g)
{
    const slong cutoff = FLINT_BIT_COUNT(f->mod.n) <= 8 ? 
        NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;

    if (FLINT_MAX(f->length, g->length) < cutoff)
        return nmod_poly_resultant_euclidean(f, g);
    else
        return nmod_poly_resultant_hgcd(f, g);
}

FLINT_DLL slong _nmod_poly_gcdinv(ulong_ptr G, ulong_ptr S,
                        ulong_srcptr A, slong lenA,
                        ulong_srcptr B, slong lenB, 
                        const nmod_t mod);

FLINT_DLL void nmod_poly_gcdinv(nmod_poly_t G, nmod_poly_t S, 
                      const nmod_poly_t A, const nmod_poly_t B);

/* Discriminant **************************************************************/

FLINT_DLL ulong _nmod_poly_discriminant(ulong_srcptr poly, slong len, nmod_t mod);

FLINT_DLL ulong nmod_poly_discriminant(const nmod_poly_t f);

/* Square roots **************************************************************/

FLINT_DLL void _nmod_poly_invsqrt_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_invsqrt_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_sqrt_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_sqrt_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL int _nmod_poly_sqrt(ulong_ptr s, ulong_srcptr p, slong len, nmod_t mod);

FLINT_DLL int nmod_poly_sqrt(nmod_poly_t b, const nmod_poly_t a);

/* Power sums ****************************************************************/

FLINT_DLL void _nmod_poly_power_sums_naive(ulong_ptr res, ulong_srcptr poly, slong len, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_power_sums_naive(nmod_poly_t res, const nmod_poly_t poly, slong n);

FLINT_DLL void _nmod_poly_power_sums_schoenhage(ulong_ptr res, ulong_srcptr poly, slong len, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_power_sums_schoenhage(nmod_poly_t res, const nmod_poly_t poly, slong n);

FLINT_DLL void _nmod_poly_power_sums(ulong_ptr res, ulong_srcptr poly, slong len, slong n, nmod_t mod);

FLINT_DLL void nmod_poly_power_sums(nmod_poly_t res, const nmod_poly_t poly, slong n);

FLINT_DLL void _nmod_poly_power_sums_to_poly_naive(ulong_ptr res, ulong_srcptr poly, slong len, nmod_t mod);

FLINT_DLL void nmod_poly_power_sums_to_poly_naive(nmod_poly_t res, const nmod_poly_t Q);

FLINT_DLL void _nmod_poly_power_sums_to_poly_schoenhage(ulong_ptr res, ulong_srcptr poly, slong len, nmod_t mod);

FLINT_DLL void nmod_poly_power_sums_to_poly_schoenhage(nmod_poly_t res, const nmod_poly_t Q);

FLINT_DLL void _nmod_poly_power_sums_to_poly(ulong_ptr res, ulong_srcptr poly, slong len, nmod_t mod);

FLINT_DLL void nmod_poly_power_sums_to_poly(nmod_poly_t res, const nmod_poly_t Q);

/* Transcendental functions **************************************************/

FLINT_DLL void _nmod_poly_atan_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_atan_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_tan_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_tan_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_asin_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_asin_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_sin_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_sin_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_cos_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_cos_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_asinh_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_asinh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_atanh_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_atanh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_sinh_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_sinh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_cosh_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_cosh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_tanh_series(ulong_ptr g, ulong_srcptr h, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_tanh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

FLINT_DLL void _nmod_poly_log_series_monomial_ui(ulong_ptr res, ulong coeff,
                ulong power, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_log_series_monomial_ui(nmod_poly_t res, ulong coeff,
                ulong power, slong n);

FLINT_DLL void _nmod_poly_log_series(ulong_ptr res, ulong_srcptr f, slong flen, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_log_series(nmod_poly_t res, const nmod_poly_t f, slong n);

FLINT_DLL void _nmod_poly_exp_series_monomial_ui(ulong_ptr res, ulong coeff,
                ulong power, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_exp_series_monomial_ui(nmod_poly_t res, ulong coeff,
                ulong power, slong n);

FLINT_DLL void _nmod_poly_exp_series_basecase(ulong_ptr f, ulong_srcptr h, slong hlen, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_exp_series_basecase(nmod_poly_t f, const nmod_poly_t h, slong n);
FLINT_DLL void  _nmod_poly_exp_expinv_series(ulong_ptr f, ulong_ptr g, ulong_srcptr h, slong hlen, slong n, nmod_t mod);
FLINT_DLL void _nmod_poly_exp_series(ulong_ptr f, ulong_srcptr h, slong hlen, slong n, nmod_t mod);
FLINT_DLL void nmod_poly_exp_series(nmod_poly_t f, const nmod_poly_t h, slong n);

/* Products  *****************************************************************/

FLINT_DLL void nmod_poly_product_roots_nmod_vec(nmod_poly_t poly, ulong_srcptr xs, slong n);

FLINT_DLL void _nmod_poly_product_roots_nmod_vec(ulong_ptr poly,
    ulong_srcptr xs, slong n, nmod_t mod);


FLINT_DLL void _nmod_poly_split_rabin(nmod_poly_t a, nmod_poly_t b,
                           const nmod_poly_t f, nmod_poly_t t, nmod_poly_t t2,
                                                       flint_rand_t randstate);

FLINT_DLL int nmod_poly_find_distinct_nonzero_roots(ulong_ptr roots,
                                                          const nmod_poly_t P);


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
    slong temp2loc; /* index of another tempory used in run */
    int good;   /* the moduli are good for CRT, essentially relatively prime */
} nmod_poly_multi_crt_struct;

typedef nmod_poly_multi_crt_struct nmod_poly_multi_crt_t[1];

FLINT_DLL void nmod_poly_multi_crt_init(nmod_poly_multi_crt_t CRT);

FLINT_DLL int nmod_poly_multi_crt_precompute(nmod_poly_multi_crt_t CRT,
                                   const nmod_poly_struct * moduli, slong len);

FLINT_DLL int nmod_poly_multi_crt_precompute_p(nmod_poly_multi_crt_t CRT,
                           const nmod_poly_struct * const * moduli, slong len);

FLINT_DLL void nmod_poly_multi_crt_precomp(nmod_poly_t output,
             const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * values);

FLINT_DLL void nmod_poly_multi_crt_precomp_p(nmod_poly_t output,
     const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * const * values);

FLINT_DLL int nmod_poly_multi_crt(nmod_poly_t output,
  const nmod_poly_struct * moduli, const nmod_poly_struct * values, slong len);

FLINT_DLL void nmod_poly_multi_crt_clear(nmod_poly_multi_crt_t CRT);

NMOD_POLY_INLINE
slong _nmod_poly_multi_crt_local_size(const nmod_poly_multi_crt_t CRT)
{
    return CRT->localsize;
}

FLINT_DLL void _nmod_poly_multi_crt_run(nmod_poly_struct * outputs,
             const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * inputs);

FLINT_DLL void _nmod_poly_multi_crt_run_p(nmod_poly_struct * outputs,
     const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * const * inputs);

/* Inflation and deflation ***************************************************/

FLINT_DLL ulong nmod_poly_deflation(const nmod_poly_t input);

FLINT_DLL void nmod_poly_deflate(nmod_poly_t result, const nmod_poly_t input,
    ulong deflation);

FLINT_DLL void nmod_poly_inflate(nmod_poly_t result, const nmod_poly_t input,
    ulong inflation);

/* Characteristic polynomial and minimal polynomial */

FLINT_DLL void _nmod_mat_charpoly_berkowitz(ulong_ptr p, const nmod_mat_t M, nmod_t mod);
FLINT_DLL void nmod_mat_charpoly_berkowitz(nmod_poly_t p, const nmod_mat_t M);
FLINT_DLL void nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t M);
FLINT_DLL void nmod_mat_charpoly(nmod_poly_t p, const nmod_mat_t M);

FLINT_DLL void nmod_mat_minpoly_with_gens(nmod_poly_t p, const nmod_mat_t X, ulong_ptr P);

FLINT_DLL void nmod_mat_minpoly(nmod_poly_t p, const nmod_mat_t M);

/* Berlekamp-Massey Algorithm - see nmod_poly/berlekamp_massey.c for more info ************/
typedef struct {
    slong npoints;
    nmod_poly_t R0, R1;
    nmod_poly_t V0, V1;
    nmod_poly_t qt, rt;
    nmod_poly_t points;
} nmod_berlekamp_massey_struct;
typedef nmod_berlekamp_massey_struct nmod_berlekamp_massey_t[1];

FLINT_DLL void nmod_berlekamp_massey_init(
                    nmod_berlekamp_massey_t B,
                    ulong p);

FLINT_DLL void nmod_berlekamp_massey_start_over(
                    nmod_berlekamp_massey_t B);

FLINT_DLL void nmod_berlekamp_massey_clear(
                    nmod_berlekamp_massey_t B);

FLINT_DLL void nmod_berlekamp_massey_set_prime(
                    nmod_berlekamp_massey_t B,
                    ulong p);

FLINT_DLL void nmod_berlekamp_massey_print(
                    const nmod_berlekamp_massey_t B);

FLINT_DLL void nmod_berlekamp_massey_add_points(
                    nmod_berlekamp_massey_t B,
                    ulong_srcptr a,
                    slong count);

FLINT_DLL void nmod_berlekamp_massey_add_zeros(
                    nmod_berlekamp_massey_t B,
                    slong count);

FLINT_DLL void nmod_berlekamp_massey_add_point(
                    nmod_berlekamp_massey_t B,
                    ulong a);

FLINT_DLL int nmod_berlekamp_massey_reduce(
                    nmod_berlekamp_massey_t B);

NMOD_POLY_INLINE ulong_srcptr nmod_berlekamp_massey_points(
                    const nmod_berlekamp_massey_t B)
{
    return B->points->coeffs;
}

NMOD_POLY_INLINE slong nmod_berlekamp_massey_point_count(
                    const nmod_berlekamp_massey_t B)
{
    return B->points->length;
}

NMOD_POLY_INLINE const nmod_poly_struct * nmod_berlekamp_massey_V_poly(
                    const nmod_berlekamp_massey_t B)
{
    return B->V1;
}

NMOD_POLY_INLINE const nmod_poly_struct * nmod_berlekamp_massey_R_poly(
                    const nmod_berlekamp_massey_t B)
{
    return B->R1;
}

#ifdef __cplusplus
}
#endif

#endif
