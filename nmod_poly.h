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

    Copyright (C) 2007, David Howden
    Copyright (C) 2010, 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Ashish Kedia

******************************************************************************/

#ifndef NMOD_POLY_H
#define NMOD_POLY_H

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "fmpz.h"

#ifdef __cplusplus
    extern "C" {
#endif

#define NMOD_DIVREM_DIVCONQUER_CUTOFF  300
#define NMOD_DIV_DIVCONQUER_CUTOFF     300 /* Must be <= NMOD_DIVREM_DIVCONQUER_CUTOFF */

#define NMOD_POLY_HGCD_CUTOFF  100      /* HGCD: Basecase -> Recursion      */
#define NMOD_POLY_GCD_CUTOFF  340       /* GCD:  Euclidean -> HGCD          */
#define NMOD_POLY_SMALL_GCD_CUTOFF 200  /* GCD (small n): Euclidean -> HGCD */

static __inline__
slong NMOD_DIVREM_BC_ITCH(slong lenA, slong lenB, nmod_t mod)
{
    const mp_bitcnt_t bits = 
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);
    
    if (bits <= FLINT_BITS)
        return lenA;
    else if (bits <= 2 * FLINT_BITS)
        return 2*(lenA + lenB - 1);
    else
        return 3*(lenA + lenB - 1);
}

static __inline__
slong NMOD_DIV_BC_ITCH(slong lenA, slong lenB, nmod_t mod)
{
    const mp_bitcnt_t bits = 
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);
    
    if (bits <= FLINT_BITS)
        return lenA - lenB + 1;
    else if (bits <= 2 * FLINT_BITS)
        return 2*lenA;
    else
        return 3*lenA;
}

static __inline__
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
    mp_ptr coeffs;
    slong alloc;
    slong length;
    nmod_t mod;
} nmod_poly_struct;

typedef nmod_poly_struct nmod_poly_t[1];

/* zn_poly helper functions  ************************************************

Copyright (C) 2007, 2008 David Harvey

*/

static __inline__
int signed_mpn_sub_n(mp_ptr res, mp_srcptr op1, mp_srcptr op2, slong n)
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

/* Memory management  ********************************************************/

void nmod_poly_init(nmod_poly_t poly, mp_limb_t n);

void nmod_poly_init_preinv(nmod_poly_t poly, mp_limb_t n, mp_limb_t ninv);

void nmod_poly_init2(nmod_poly_t poly, mp_limb_t n, slong alloc);

void nmod_poly_init2_preinv(nmod_poly_t poly, 
                                      mp_limb_t n, mp_limb_t ninv, slong alloc);

void nmod_poly_realloc(nmod_poly_t poly, slong alloc);

void nmod_poly_clear(nmod_poly_t poly);

void nmod_poly_fit_length(nmod_poly_t poly, slong alloc);

static __inline__
void _nmod_poly_set_length(nmod_poly_t poly, slong len)
{
    poly->length = len;
}

static __inline__
void _nmod_poly_normalise(nmod_poly_t poly)
{
    while (poly->length && (poly->coeffs[poly->length - 1] == WORD(0)))
        poly->length--;
}

/* Polynomial parameters  ****************************************************/

static __inline__
slong nmod_poly_length(const nmod_poly_t poly)
{
    return poly->length;
}

static __inline__
slong nmod_poly_degree(const nmod_poly_t poly)
{
    return poly->length - 1;
}

static __inline__
mp_limb_t nmod_poly_modulus(const nmod_poly_t poly)
{
    return poly->mod.n;
}

static __inline__
mp_bitcnt_t nmod_poly_max_bits(const nmod_poly_t poly)
{
    return _nmod_vec_max_bits(poly->coeffs, poly->length);
}

static __inline__
mp_ptr nmod_poly_lead(const nmod_poly_t poly)
{
    if (poly->length)
        return poly->coeffs + (poly->length - 1);
    else
        return NULL;
}

/* Assignment and basic manipulation  ****************************************/

static __inline__
void nmod_poly_set(nmod_poly_t a, const nmod_poly_t b)
{
    if (a != b)
    {
        nmod_poly_fit_length(a, b->length);
        flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
        a->length = b->length;
    }
}

static __inline__
void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)
{
    slong t;
    mp_ptr tp;

    t = poly1->alloc;
    poly1->alloc = poly2->alloc;
    poly2->alloc = t;

    t = poly1->length;
    poly1->length = poly2->length;
    poly2->length = t;

    tp = poly1->coeffs;
    poly1->coeffs = poly2->coeffs;
    poly2->coeffs = tp;
}

static __inline__
void nmod_poly_zero(nmod_poly_t res)
{
    res->length = 0;
}

static __inline__
void nmod_poly_one(nmod_poly_t res)
{
    nmod_poly_fit_length(res, 1);
    res->length = 1;
    res->coeffs[0] = 1;
}

static __inline__
void nmod_poly_truncate(nmod_poly_t poly, slong len)
{
    if (poly->length > len)
    {
        poly->length = len;
        _nmod_poly_normalise(poly);
    }
}

void _nmod_poly_reverse(mp_ptr output, mp_srcptr input, slong len, slong m);

void nmod_poly_reverse(nmod_poly_t output, const nmod_poly_t input, slong m);

/* Comparison  ***************************************************************/

static __inline__
int nmod_poly_equal(const nmod_poly_t a, const nmod_poly_t b)
{
    if (a->length != b->length)
        return 0;

    if (a != b)
        if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
            return 0;

   return 1;
}

static __inline__
int nmod_poly_is_zero(const nmod_poly_t poly)
{
    return (poly->length == 0);
}

static __inline__ int
nmod_poly_is_one(const nmod_poly_t poly)
{
    return (poly->length == 1) && (poly->coeffs[0] == 1);
}

/* Randomisation  ************************************************************/

void nmod_poly_randtest(nmod_poly_t poly, flint_rand_t state, slong len);

static __inline__ void
nmod_poly_randtest_not_zero(nmod_poly_t poly, flint_rand_t state, slong len)
{
    do {
        nmod_poly_randtest(poly, state, len);
    } while (nmod_poly_is_zero(poly));
}
      
void
nmod_poly_randtest_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);

void
nmod_poly_randtest_monic(nmod_poly_t poly, flint_rand_t state, slong len);

void
nmod_poly_randtest_monic_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);

void
nmod_poly_randtest_trinomial(nmod_poly_t poly, flint_rand_t state, slong len);

int
nmod_poly_randtest_trinomial_irreducible(nmod_poly_t poly, flint_rand_t state,
                                         slong len, slong max_attempts);

void
nmod_poly_randtest_pentomial(nmod_poly_t poly, flint_rand_t state, slong len);

int
nmod_poly_randtest_pentomial_irreducible(nmod_poly_t poly, flint_rand_t state,
                                         slong len, slong max_attempts);

void
nmod_poly_randtest_sparse_irreducible(nmod_poly_t poly, flint_rand_t state, slong len);

/* Getting and setting coefficients  *****************************************/

static __inline__
ulong nmod_poly_get_coeff_ui(const nmod_poly_t poly, slong j)
{
    return (j >= poly->length) ? 0 : poly->coeffs[j];
}

void nmod_poly_set_coeff_ui(nmod_poly_t poly, slong j, ulong c);

/* Input and output  *********************************************************/

char * nmod_poly_get_str(const nmod_poly_t poly);

int nmod_poly_set_str(nmod_poly_t poly, const char * s);

static __inline__
int nmod_poly_print(const nmod_poly_t a)
{
    size_t r;
    slong i;

    r = flint_printf("%wd %wu", a->length, a->mod.n);

    if (a->length == 0)
        return r;
    else
        if (r > 0)
            r = flint_printf(" ");

    for (i = 0; (r > 0) && (i < a->length); i++)
        r = flint_printf(" %wu", a->coeffs[i]);

    return (int) r;
}

int nmod_poly_fread(FILE * f, nmod_poly_t poly);

static __inline__
int nmod_poly_fprint(FILE * f, const nmod_poly_t poly)
{
    char *s;
    int r;

    s = nmod_poly_get_str(poly);
    r = fputs(s, f);
    flint_free(s);

    return (r < 0) ? r : 1;
}

static __inline__
int nmod_poly_read(nmod_poly_t poly)
{
    return nmod_poly_fread(stdin, poly);
}

/* Shifting  *****************************************************************/

void _nmod_poly_shift_left(mp_ptr res, mp_srcptr poly, slong len, slong k);

void nmod_poly_shift_left(nmod_poly_t res, const nmod_poly_t poly, slong k);

void _nmod_poly_shift_right(mp_ptr res, mp_srcptr poly, slong len, slong k);

void nmod_poly_shift_right(nmod_poly_t res, const nmod_poly_t poly, slong k);

/* Addition and subtraction  *************************************************/

void _nmod_poly_add(mp_ptr res, mp_srcptr poly1, slong len1, 
                                       mp_srcptr poly2, slong len2, nmod_t mod);

void nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1, 
                                                      const nmod_poly_t poly2);

void _nmod_poly_sub(mp_ptr res, mp_srcptr poly1, slong len1, 
                                       mp_srcptr poly2, slong len2, nmod_t mod);

void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, 
                                                      const nmod_poly_t poly2);

void nmod_poly_neg(nmod_poly_t res, const nmod_poly_t poly1);

/* Scalar multiplication and division  ***************************************/

void nmod_poly_scalar_mul_nmod(nmod_poly_t res, 
                                         const nmod_poly_t poly1, mp_limb_t c);

void _nmod_poly_make_monic(mp_ptr output, 
                                   mp_srcptr input, slong len, nmod_t mod);

void nmod_poly_make_monic(nmod_poly_t output, const nmod_poly_t input);

/* Bit packing and unpacking aand reduction  **********************************/

void _nmod_poly_KS2_pack1(mp_ptr res, mp_srcptr op, slong n, slong s,
                                                    ulong b, ulong k, slong r);

void _nmod_poly_KS2_pack(mp_ptr res, mp_srcptr op, slong n, slong s,
                                                    ulong b, ulong k, slong r);

void _nmod_poly_KS2_unpack1(mp_ptr res, mp_srcptr op, slong n, ulong b,
                                                                      ulong k);

void _nmod_poly_KS2_unpack2(mp_ptr res, mp_srcptr op, slong n, ulong b,
                                                                      ulong k);

void _nmod_poly_KS2_unpack3(mp_ptr res, mp_srcptr op, slong n, ulong b,
                                                                      ulong k);

void _nmod_poly_KS2_unpack(mp_ptr res, mp_srcptr op, slong n, ulong b,
                                                                      ulong k);

void _nmod_poly_KS2_reduce(mp_ptr res, slong s, mp_srcptr op, 
                                                 slong n, ulong w, nmod_t mod);

void _nmod_poly_KS2_recover_reduce1(mp_ptr res, slong s, mp_srcptr op1,
                                  mp_srcptr op2, slong n, ulong b, nmod_t mod);

void _nmod_poly_KS2_recover_reduce1(mp_ptr res, slong s, mp_srcptr op1,
                                  mp_srcptr op2, slong n, ulong b, nmod_t mod);

void _nmod_poly_KS2_recover_reduce2(mp_ptr res, slong s, mp_srcptr op1,
                                  mp_srcptr op2, slong n, ulong b, nmod_t mod);

void _nmod_poly_KS2_recover_reduce2b(mp_ptr res, slong s, mp_srcptr op1,
                                  mp_srcptr op2, slong n, ulong b, nmod_t mod);

void _nmod_poly_KS2_recover_reduce(mp_ptr res, slong s, mp_srcptr op1,
                                  mp_srcptr op2, slong n, ulong b, nmod_t mod);

void _nmod_poly_bit_pack(mp_ptr res, mp_srcptr poly, 
                                                  slong len, mp_bitcnt_t bits);

void _nmod_poly_bit_unpack(mp_ptr res, slong len, 
                                  mp_srcptr mpn, mp_bitcnt_t bits, nmod_t mod);

void nmod_poly_bit_pack(fmpz_t f, const nmod_poly_t poly,
                   mp_bitcnt_t bit_size);

void
nmod_poly_bit_unpack(nmod_poly_t poly, const fmpz_t f, mp_bitcnt_t bit_size);

void
_nmod_poly_KS2_pack(mp_ptr res, mp_srcptr op, slong n, slong s,
               ulong b, ulong k, slong r);

/* Multiplication  ***********************************************************/

void _nmod_poly_mul_classical(mp_ptr res, mp_srcptr poly1, slong len1, 
                                       mp_srcptr poly2, slong len2, nmod_t mod);

void nmod_poly_mul_classical(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow_classical(mp_ptr res, mp_srcptr poly1, slong len1, 
                           mp_srcptr poly2, slong len2, slong trunc, nmod_t mod);

void nmod_poly_mullow_classical(nmod_poly_t res, 
                  const nmod_poly_t poly1, const nmod_poly_t poly2, slong trunc);

void _nmod_poly_mulhigh_classical(mp_ptr res, mp_srcptr poly1, slong len1, 
                           mp_srcptr poly2, slong len2, slong start, nmod_t mod);

void nmod_poly_mulhigh_classical(nmod_poly_t res, 
                  const nmod_poly_t poly1, const nmod_poly_t poly2, slong start);

void _nmod_poly_mul_KS(mp_ptr out, mp_srcptr in1, slong len1, 
                        mp_srcptr in2, slong len2, mp_bitcnt_t bits, nmod_t mod);

void nmod_poly_mul_KS(nmod_poly_t res, 
             const nmod_poly_t poly1, const nmod_poly_t poly2, mp_bitcnt_t bits);

void _nmod_poly_mul_KS2(mp_ptr res, mp_srcptr op1, slong n1,
                                            mp_srcptr op2, slong n2, nmod_t mod);

void nmod_poly_mul_KS2(nmod_poly_t res,
                               const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mul_KS4(mp_ptr res, mp_srcptr op1, slong n1,
                                            mp_srcptr op2, slong n2, nmod_t mod);

void nmod_poly_mul_KS4(nmod_poly_t res,
                               const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow_KS(mp_ptr out, mp_srcptr in1, slong len1,
               mp_srcptr in2, slong len2, mp_bitcnt_t bits, slong n, nmod_t mod);

void nmod_poly_mullow_KS(nmod_poly_t res, const nmod_poly_t poly1, 
                             const nmod_poly_t poly2, mp_bitcnt_t bits, slong n);

void _nmod_poly_mul(mp_ptr res, mp_srcptr poly1, slong len1, 
                                       mp_srcptr poly2, slong len2, nmod_t mod);

void nmod_poly_mul(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow(mp_ptr res, mp_srcptr poly1, slong len1, 
                           mp_srcptr poly2, slong len2, slong trunc, nmod_t mod);

void nmod_poly_mullow(nmod_poly_t res, const nmod_poly_t poly1, 
                                          const nmod_poly_t poly2, slong trunc);

void _nmod_poly_mulhigh(mp_ptr res, mp_srcptr poly1, slong len1, 
                               mp_srcptr poly2, slong len2, slong n, nmod_t mod);

void nmod_poly_mulhigh(nmod_poly_t res, const nmod_poly_t poly1, 
                                              const nmod_poly_t poly2, slong n);

void _nmod_poly_mulmod(mp_ptr res, mp_srcptr poly1, slong len1, 
                             mp_srcptr poly2, slong len2, mp_srcptr f,
                            slong lenf, nmod_t mod);

void nmod_poly_mulmod(nmod_poly_t res,
    const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t f);

void _nmod_poly_mulmod_preinv(mp_ptr res, mp_srcptr poly1, slong len1, 
                              mp_srcptr poly2, slong len2, mp_srcptr f,
                              slong lenf, mp_srcptr finv, slong lenfinv, nmod_t mod);

void
nmod_poly_mulmod_preinv(nmod_poly_t res, const nmod_poly_t poly1,
                        const nmod_poly_t poly2, const nmod_poly_t f,
                        const nmod_poly_t finv);

int _nmod_poly_invmod(mp_limb_t *A, 
                      const mp_limb_t *B, slong lenB, 
                      const mp_limb_t *P, slong lenP, const nmod_t mod);

int nmod_poly_invmod(nmod_poly_t A, 
                     const nmod_poly_t B, const nmod_poly_t P);

/* Powering  *****************************************************************/

void _nmod_poly_pow_binexp(mp_ptr res, 
                              mp_srcptr poly, slong len, ulong e, nmod_t mod);

void nmod_poly_pow_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e);

void _nmod_poly_pow(mp_ptr res, mp_srcptr poly, slong len, ulong e, nmod_t mod);

void nmod_poly_pow(nmod_poly_t res, const nmod_poly_t poly, ulong e);

void _nmod_poly_pow_trunc_binexp(mp_ptr res, mp_srcptr poly, 
                                              ulong e, slong trunc, nmod_t mod);

void nmod_poly_pow_trunc_binexp(nmod_poly_t res, 
                                  const nmod_poly_t poly, ulong e, slong trunc);

void _nmod_poly_pow_trunc(mp_ptr res, mp_srcptr poly, 
                                              ulong e, slong trunc, nmod_t mod);

void nmod_poly_pow_trunc(nmod_poly_t res, 
                                  const nmod_poly_t poly, ulong e, slong trunc);

void
nmod_poly_powmod_ui_binexp(nmod_poly_t res, 
                           const nmod_poly_t poly, ulong e,
                           const nmod_poly_t f);

void
_nmod_poly_powmod_ui_binexp(mp_ptr res, mp_srcptr poly, 
                                ulong e, mp_srcptr f, slong lenf, nmod_t mod);

void
_nmod_poly_powmod_mpz_binexp(mp_ptr res, mp_srcptr poly, 
                                mpz_srcptr e, mp_srcptr f,
                                slong lenf, nmod_t mod);

void
nmod_poly_powmod_mpz_binexp(nmod_poly_t res, 
                           const nmod_poly_t poly, mpz_srcptr e,
                           const nmod_poly_t f);

void
_nmod_poly_powmod_ui_binexp_preinv (mp_ptr res, mp_srcptr poly,
                                    ulong e, mp_srcptr f, slong lenf,
                                    mp_srcptr finv, slong lenfinv, nmod_t mod);

void
nmod_poly_powmod_ui_binexp_preinv(nmod_poly_t res, 
                           const nmod_poly_t poly, ulong e,
                           const nmod_poly_t f, const nmod_poly_t finv);

void
_nmod_poly_powmod_x_ui_preinv (mp_ptr res, ulong e, mp_srcptr f, slong lenf,
                               mp_srcptr finv, slong lenfinv, nmod_t mod);

void
nmod_poly_powmod_x_ui_preinv(nmod_poly_t res, ulong e, const nmod_poly_t f,
                             const nmod_poly_t finv);

void
_nmod_poly_powmod_mpz_binexp_preinv (mp_ptr res, mp_srcptr poly,
                                    mpz_srcptr e, mp_srcptr f, slong lenf,
                                    mp_srcptr finv, slong lenfinv, nmod_t mod);

void
nmod_poly_powmod_mpz_binexp_preinv(nmod_poly_t res,
                           const nmod_poly_t poly, mpz_srcptr e,
                           const nmod_poly_t f, const nmod_poly_t finv);

/* Division  *****************************************************************/

void _nmod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, mp_ptr W,
                 mp_srcptr A, slong A_len, mp_srcptr B, slong B_len, nmod_t mod);

void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, 
                                     const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_divrem_divconquer_recursive(mp_ptr Q, mp_ptr BQ, mp_ptr W,  
                    mp_ptr V, mp_srcptr A, mp_srcptr B, slong lenB, nmod_t mod);

void _nmod_poly_divrem_divconquer(mp_ptr Q, mp_ptr R, 
                   mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_divrem_divconquer(nmod_poly_t Q, nmod_poly_t R,
                                     const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_divrem_q0(mp_ptr Q, mp_ptr R, 
                          mp_srcptr A, mp_srcptr B, slong lenA, nmod_t mod);

void _nmod_poly_divrem_q1(mp_ptr Q, mp_ptr R, 
                          mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                          nmod_t mod);

void _nmod_poly_divrem(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA, 
                                           mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R,
                                     const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_div_basecase(mp_ptr Q, mp_ptr W, mp_srcptr A, slong A_len, 
                                          mp_srcptr B, slong B_len, nmod_t mod);

void nmod_poly_div_basecase(nmod_poly_t Q, const nmod_poly_t A,
                                                          const nmod_poly_t B);

void _nmod_poly_div_divconquer_recursive(mp_ptr Q, mp_ptr W, mp_ptr V,
                              mp_srcptr A, mp_srcptr B, slong lenB, nmod_t mod);

void _nmod_poly_div_divconquer(mp_ptr Q, mp_srcptr A, slong lenA, 
                                           mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_div_divconquer(nmod_poly_t Q,
                                     const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_div(mp_ptr Q, mp_srcptr A, slong lenA, 
                                           mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_div(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_rem_basecase(mp_ptr R, mp_ptr W, mp_srcptr A, slong lenA, 
                                       mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_rem_basecase(nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_rem_q1(mp_ptr R, 
                       mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                       nmod_t mod);

void _nmod_poly_rem(mp_ptr R, mp_srcptr A, slong lenA, 
                              mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_rem(nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_inv_series_basecase(mp_ptr Qinv, 
                                              mp_srcptr Q, slong n, nmod_t mod);

void nmod_poly_inv_series_basecase(nmod_poly_t Qinv, 
                                                  const nmod_poly_t Q, slong n);

void _nmod_poly_inv_series_newton(mp_ptr Qinv, 
                                              mp_srcptr Q, slong n, nmod_t mod);

void nmod_poly_inv_series_newton(nmod_poly_t Qinv, 
                                                  const nmod_poly_t Q, slong n);

static __inline__
void _nmod_poly_inv_series(mp_ptr Qinv, mp_srcptr Q, slong n, nmod_t mod)
{
    _nmod_poly_inv_series_newton(Qinv, Q, n, mod);
}

static __inline__
void nmod_poly_inv_series(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)
{
    nmod_poly_inv_series_newton(Qinv, Q, n);
}

void _nmod_poly_div_series(mp_ptr Q, mp_srcptr A, mp_srcptr B, 
                                                          slong n, nmod_t mod);

void nmod_poly_div_series(nmod_poly_t Q, const nmod_poly_t A, 
                                                 const nmod_poly_t B, slong n);

void _nmod_poly_div_newton(mp_ptr Q, mp_srcptr A, slong Alen, 
                                          mp_srcptr B, slong Blen, nmod_t mod);

void nmod_poly_div_newton(nmod_poly_t Q, const nmod_poly_t A,
                                                         const nmod_poly_t B);

void _nmod_poly_divrem_newton(mp_ptr Q, mp_ptr R, 
                  mp_srcptr A, slong Alen, mp_srcptr B, slong Blen, nmod_t mod);

void nmod_poly_divrem_newton(nmod_poly_t Q, nmod_poly_t R, 
                                    const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_div_newton_n_preinv (mp_ptr Q, mp_srcptr A, slong lenA,
            mp_srcptr B, slong lenB, mp_srcptr Binv, slong lenBinv, nmod_t mod);

void nmod_poly_div_newton_n_preinv (nmod_poly_t Q, const nmod_poly_t A,
                                 const nmod_poly_t B, const nmod_poly_t Binv);

void _nmod_poly_divrem_newton_n_preinv (mp_ptr Q, mp_ptr R, mp_srcptr A,
 slong lenA, mp_srcptr B, slong lenB, mp_srcptr Binv, slong lenBinv, nmod_t mod);

void nmod_poly_divrem_newton_n_preinv(nmod_poly_t Q, nmod_poly_t R,
            const nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t Binv);

mp_limb_t
_nmod_poly_div_root(mp_ptr Q, mp_srcptr A, slong len, mp_limb_t c, nmod_t mod);

mp_limb_t
nmod_poly_div_root(nmod_poly_t Q, const nmod_poly_t A, mp_limb_t c);

/* Derivative  ***************************************************************/

void _nmod_poly_derivative(mp_ptr x_prime, mp_srcptr x, slong len, nmod_t mod);

void nmod_poly_derivative(nmod_poly_t x_prime, const nmod_poly_t x);

void _nmod_poly_integral(mp_ptr x_int, mp_srcptr x, slong len, nmod_t mod);

void nmod_poly_integral(nmod_poly_t x_int, const nmod_poly_t x);

/* Evaluation  ***************************************************************/
void
_nmod_poly_evaluate_fmpz(fmpz_t rop, const mp_srcptr poly, const slong len, const fmpz_t c);

void
nmod_poly_evaluate_fmpz(fmpz_t rop, const nmod_poly_t poly, const fmpz_t c);

mp_limb_t _nmod_poly_evaluate_nmod(mp_srcptr poly, 
                                           slong len, mp_limb_t c, nmod_t mod);

mp_limb_t nmod_poly_evaluate_nmod(const nmod_poly_t poly, mp_limb_t c);

void _nmod_poly_evaluate_nmod_vec(mp_ptr ys, mp_srcptr coeffs, slong len,
    mp_srcptr xs, slong n, nmod_t mod);

void nmod_poly_evaluate_nmod_vec(mp_ptr ys,
        const nmod_poly_t poly, mp_srcptr xs, slong n);

void nmod_poly_evaluate_nmod_vec(mp_ptr ys, const nmod_poly_t poly,
        mp_srcptr xs, slong n);

void _nmod_poly_evaluate_nmod_vec_iter(mp_ptr ys, mp_srcptr coeffs, slong len,
    mp_srcptr xs, slong n, nmod_t mod);

void nmod_poly_evaluate_nmod_vec_iter(mp_ptr ys,
    const nmod_poly_t poly, mp_srcptr xs, slong n);

void
_nmod_poly_evaluate_nmod_vec_fast_precomp(mp_ptr vs, mp_srcptr poly,
    slong plen, const mp_ptr * tree, slong len, nmod_t mod);

void _nmod_poly_evaluate_nmod_vec_fast(mp_ptr ys, mp_srcptr coeffs, slong len,
    mp_srcptr xs, slong n, nmod_t mod);

void nmod_poly_evaluate_nmod_vec_fast(mp_ptr ys,
        const nmod_poly_t poly, mp_srcptr xs, slong n);

void _nmod_poly_evaluate_mat(nmod_mat_t dest,
	const mp_srcptr poly, slong len, const nmod_mat_t c);

void nmod_poly_evaluate_mat(nmod_mat_t dest,
	const nmod_poly_t poly, const nmod_mat_t c);

/* Subproduct tree  **********************************************************/

mp_ptr * _nmod_poly_tree_alloc(slong len);

void _nmod_poly_tree_free(mp_ptr * tree, slong len);

void _nmod_poly_tree_build(mp_ptr * tree, mp_srcptr roots,
    slong len, nmod_t mod);

/* Interpolation  ************************************************************/

void _nmod_poly_interpolate_nmod_vec_newton(mp_ptr poly, mp_srcptr xs,
                        mp_srcptr ys, slong n, nmod_t mod);

void nmod_poly_interpolate_nmod_vec_newton(nmod_poly_t poly,
                        mp_srcptr xs, mp_srcptr ys, slong n);

void _nmod_poly_interpolate_nmod_vec_barycentric(mp_ptr poly, mp_srcptr xs,
                        mp_srcptr ys, slong n, nmod_t mod);

void nmod_poly_interpolate_nmod_vec_barycentric(nmod_poly_t poly,
                        mp_srcptr xs, mp_srcptr ys, slong n);

void _nmod_poly_interpolate_nmod_vec(mp_ptr poly, mp_srcptr xs,
                        mp_srcptr ys, slong n, nmod_t mod);

void nmod_poly_interpolate_nmod_vec(nmod_poly_t poly,
                        mp_srcptr xs, mp_srcptr ys, slong n);

void nmod_poly_interpolate_nmod_vec_fast(nmod_poly_t poly,
                                    mp_srcptr xs, mp_srcptr ys, slong n);

void _nmod_poly_interpolate_nmod_vec_fast(mp_ptr poly,
                            mp_srcptr xs, mp_srcptr ys, slong len, nmod_t mod);

void
_nmod_poly_interpolate_nmod_vec_fast_precomp(mp_ptr poly, mp_srcptr ys,
    const mp_ptr * tree, mp_srcptr weights, slong len, nmod_t mod);

void _nmod_poly_interpolation_weights(mp_ptr w, const mp_ptr * tree,
    slong len, nmod_t mod);

/* Composition  **************************************************************/

void _nmod_poly_compose_horner(mp_ptr res, mp_srcptr poly1, 
                            slong len1, mp_srcptr poly2, slong len2, nmod_t mod);

void nmod_poly_compose_horner(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_compose_divconquer(mp_ptr res, mp_srcptr poly1, slong len1, 
                                       mp_srcptr poly2, slong len2, nmod_t mod);
void nmod_poly_compose_divconquer(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_compose(mp_ptr res, mp_srcptr poly1, slong len1, 
                                       mp_srcptr poly2, slong len2, nmod_t mod);

void nmod_poly_compose(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

/* Taylor shift  *************************************************************/

void _nmod_poly_taylor_shift_horner(mp_ptr poly, mp_limb_t c,
    slong len, nmod_t mod);

void nmod_poly_taylor_shift_horner(nmod_poly_t g,
    const nmod_poly_t f, mp_limb_t c);

void _nmod_poly_taylor_shift_convolution(mp_ptr poly, mp_limb_t c,
    slong len, nmod_t mod);

void nmod_poly_taylor_shift_convolution(nmod_poly_t g,
    const nmod_poly_t f, mp_limb_t c);

void _nmod_poly_taylor_shift(mp_ptr poly, mp_limb_t c, slong len, nmod_t mod);

void nmod_poly_taylor_shift(nmod_poly_t g, const nmod_poly_t f, mp_limb_t c);

/* Modular composition  ******************************************************/

void
_nmod_poly_compose_mod_brent_kung(mp_ptr res, mp_srcptr f, slong lenf,
                            mp_srcptr g, mp_srcptr h, slong lenh, nmod_t mod);

void
nmod_poly_compose_mod_brent_kung(nmod_poly_t res, 
                    const nmod_poly_t f, const nmod_poly_t g,
                    const nmod_poly_t h);

void
_nmod_poly_reduce_matrix_mod_poly (nmod_mat_t A, const nmod_mat_t B,
                          const nmod_poly_t f);

void
_nmod_poly_precompute_matrix (nmod_mat_t A, mp_srcptr poly1, mp_srcptr poly2,
               slong len2, mp_srcptr poly2inv, slong len2inv, nmod_t mod);

void
nmod_poly_precompute_matrix (nmod_mat_t A, const nmod_poly_t poly1,
                          const nmod_poly_t poly2, const nmod_poly_t poly2inv);

void
_nmod_poly_compose_mod_brent_kung_precomp_preinv(mp_ptr res, mp_srcptr poly1,
                            slong len1, const nmod_mat_t A, mp_srcptr poly3,
                            slong len3, mp_srcptr poly3inv, slong len3inv,
                            nmod_t mod);

void
nmod_poly_compose_mod_brent_kung_precomp_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_mat_t A,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv);

void
_nmod_poly_compose_mod_brent_kung_preinv(mp_ptr res, mp_srcptr poly1, slong len1,
                            mp_srcptr poly2, mp_srcptr poly3, slong len3,
                            mp_srcptr poly3inv, slong len3inv, nmod_t mod);

void
nmod_poly_compose_mod_brent_kung_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv);

void
_nmod_poly_compose_mod_horner(mp_ptr res,
    mp_srcptr f, slong lenf, mp_srcptr g, mp_srcptr h, slong lenh, nmod_t mod);

void
nmod_poly_compose_mod_horner(nmod_poly_t res, 
                    const nmod_poly_t f, const nmod_poly_t g,
                    const nmod_poly_t h);

void
_nmod_poly_compose_mod(mp_ptr res, mp_srcptr f, slong lenf, 
                            mp_srcptr g,
                            mp_srcptr h, slong lenh, nmod_t mod);

void
nmod_poly_compose_mod(nmod_poly_t res, 
                    const nmod_poly_t f, const nmod_poly_t g,
                    const nmod_poly_t h);

/* Power series composition and reversion ************************************/

void
_nmod_poly_compose_series_horner(mp_ptr res, mp_srcptr poly1, slong len1, 
                            mp_srcptr poly2, slong len2, slong n, nmod_t mod);
void
nmod_poly_compose_series_horner(nmod_poly_t res, 
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

void
_nmod_poly_compose_series_brent_kung(mp_ptr res, mp_srcptr poly1, slong len1, 
                            mp_srcptr poly2, slong len2, slong n, nmod_t mod);
void
nmod_poly_compose_series_brent_kung(nmod_poly_t res, 
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

void
_nmod_poly_compose_series(mp_ptr res, mp_srcptr poly1, slong len1, 
                            mp_srcptr poly2, slong len2, slong n, nmod_t mod);
void
nmod_poly_compose_series(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n);

void
_nmod_poly_revert_series_lagrange(mp_ptr Qinv, mp_srcptr Q, slong n, nmod_t mod);

void
nmod_poly_revert_series_lagrange(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n);

void
_nmod_poly_revert_series_lagrange_fast(mp_ptr Qinv, mp_srcptr Q,
    slong n, nmod_t mod);

void
nmod_poly_revert_series_lagrange_fast(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n);

void
_nmod_poly_revert_series_newton(mp_ptr Qinv, mp_srcptr Q, slong n, nmod_t mod);

void
nmod_poly_revert_series_newton(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n);

void
_nmod_poly_revert_series(mp_ptr Qinv, mp_srcptr Q, slong n, nmod_t mod);

void
nmod_poly_revert_series(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n);

void
_nmod_poly_compose_series_divconquer(mp_ptr res, mp_srcptr poly1, slong len1, 
                                                 mp_srcptr poly2, slong len2, 
                                                 slong N, nmod_t mod);

void 
nmod_poly_compose_series_divconquer(nmod_poly_t res, 
    const nmod_poly_t poly1, const nmod_poly_t poly2, slong N);

/* Greatest common divisor  **************************************************/

slong _nmod_poly_gcd_euclidean(mp_ptr G, 
                   mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_gcd_euclidean(nmod_poly_t G, 
                                     const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_hgcd(mp_ptr *M, slong *lenM, 
                     mp_ptr A, slong *lenA, mp_ptr B, slong *lenB, 
                     mp_srcptr a, slong lena, mp_srcptr b, slong lenb, 
                     nmod_t mod);

slong _nmod_poly_gcd_hgcd(mp_ptr G, mp_srcptr A, slong lenA, 
                                   mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_gcd_hgcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_gcd(mp_ptr G, mp_srcptr A, slong lenA, 
                              mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_gcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_xgcd_euclidean(mp_ptr res, mp_ptr s, mp_ptr t, 
           mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, nmod_t mod);

void nmod_poly_xgcd_euclidean(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                                     const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_xgcd_hgcd(mp_ptr G, mp_ptr S, mp_ptr T, 
                          mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, 
                          nmod_t mod);

void nmod_poly_xgcd_hgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                         const nmod_poly_t A, const nmod_poly_t B);

slong _nmod_poly_xgcd(mp_ptr G, mp_ptr S, mp_ptr T, 
            mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, nmod_t mod);

void nmod_poly_xgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                                   const nmod_poly_t A, const nmod_poly_t B);

mp_limb_t 
_nmod_poly_resultant_euclidean(mp_srcptr poly1, slong len1, 
                               mp_srcptr poly2, slong len2, nmod_t mod);

mp_limb_t 
nmod_poly_resultant_euclidean(const nmod_poly_t f, const nmod_poly_t g);

static __inline__ mp_limb_t 
_nmod_poly_resultant(mp_srcptr poly1, slong len1, 
                     mp_srcptr poly2, slong len2, nmod_t mod)
{
    return _nmod_poly_resultant_euclidean(poly1, len1, poly2, len2, mod);
}

static __inline__ mp_limb_t 
nmod_poly_resultant(const nmod_poly_t f, const nmod_poly_t g)
{
    return nmod_poly_resultant_euclidean(f, g);
}

slong _nmod_poly_gcdinv(mp_limb_t *G, mp_limb_t *S, 
                        const mp_limb_t *A, slong lenA,
                        const mp_limb_t *B, slong lenB, 
                        const nmod_t mod);

void nmod_poly_gcdinv(nmod_poly_t G, nmod_poly_t S, 
                      const nmod_poly_t A, const nmod_poly_t B);

/* Square roots **************************************************************/

void _nmod_poly_invsqrt_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);

void nmod_poly_invsqrt_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_sqrt_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);

void nmod_poly_sqrt_series(nmod_poly_t g, const nmod_poly_t h, slong n);

int _nmod_poly_sqrt(mp_ptr s, mp_srcptr p, slong len, nmod_t mod);

int nmod_poly_sqrt(nmod_poly_t b, const nmod_poly_t a);

/* Transcendental functions **************************************************/

void _nmod_poly_atan_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_atan_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_tan_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_tan_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_asin_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_asin_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_sin_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_sin_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_cos_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_cos_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_asinh_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_asinh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_atanh_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_atanh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_sinh_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_sinh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_cosh_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_cosh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_tanh_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_tanh_series(nmod_poly_t g, const nmod_poly_t h, slong n);

void _nmod_poly_log_series_monomial_ui(mp_ptr res, mp_limb_t coeff,
                ulong power, slong n, nmod_t mod);
void nmod_poly_log_series_monomial_ui(nmod_poly_t res, mp_limb_t coeff,
                ulong power, slong n);

void _nmod_poly_log_series(mp_ptr res, mp_srcptr f, slong n, nmod_t mod);
void nmod_poly_log_series(nmod_poly_t res, const nmod_poly_t f, slong n);

void _nmod_poly_exp_series_monomial_ui(mp_ptr res, mp_limb_t coeff,
                ulong power, slong n, nmod_t mod);
void nmod_poly_exp_series_monomial_ui(nmod_poly_t res, mp_limb_t coeff,
                ulong power, slong n);


void
_nmod_poly_exp_series_basecase(mp_ptr f, mp_srcptr h,
                                    slong hlen, slong n, nmod_t mod);
void nmod_poly_exp_series_basecase(nmod_poly_t f, const nmod_poly_t h, slong n);

void  _nmod_poly_exp_expinv_series(mp_ptr f, mp_ptr g, mp_srcptr h, slong n, nmod_t mod);

void _nmod_poly_exp_series(mp_ptr f, mp_srcptr h, slong n, nmod_t mod);
void nmod_poly_exp_series(nmod_poly_t f, const nmod_poly_t h, slong n);

/* Products  *****************************************************************/

void
nmod_poly_product_roots_nmod_vec(nmod_poly_t poly, mp_srcptr xs, slong n);

void
_nmod_poly_product_roots_nmod_vec(mp_ptr poly,
    mp_srcptr xs, slong n, nmod_t mod);

/* Inflation and deflation ***************************************************/

ulong nmod_poly_deflation(const nmod_poly_t input);

void nmod_poly_deflate(nmod_poly_t result, const nmod_poly_t input,
    ulong deflation);

void nmod_poly_inflate(nmod_poly_t result, const nmod_poly_t input,
    ulong inflation);

#ifdef __cplusplus
    }
#endif

#include "nmod_poly_factor.h"

#endif
