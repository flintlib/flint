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

    Copyright (C) 2010 William Hart

******************************************************************************/

#ifndef NMOD_POLY_H
#define NMOD_POLY_H

#include <stdio.h>
#include <mpir.h>
#include "nmod_vec.h"
#include "ulong_extras.h"

typedef struct
{
    mp_ptr coeffs;
    long alloc;
    long length;
    nmod_t mod;
} nmod_poly_struct;

typedef nmod_poly_struct nmod_poly_t[1];

void nmod_poly_init(nmod_poly_t poly, mp_limb_t n);

void nmod_poly_init_preinv(nmod_poly_t poly, mp_limb_t n, mp_limb_t ninv);

void nmod_poly_init2(nmod_poly_t poly, mp_limb_t n, long alloc);

void nmod_poly_init2_preinv(nmod_poly_t poly, 
                                      mp_limb_t n, mp_limb_t ninv, long alloc);

void nmod_poly_realloc(nmod_poly_t poly, long alloc);

void nmod_poly_clear(nmod_poly_t poly);

void nmod_poly_fit_length(nmod_poly_t poly, long alloc);

static __inline__
mp_bitcnt_t nmod_poly_max_bits(nmod_poly_t poly)
{
    return _nmod_vec_max_bits(poly->coeffs, poly->length);
}

static __inline__
void _nmod_poly_normalise(nmod_poly_t poly)
{
    while (poly->length && (poly->coeffs[poly->length - 1] == 0L))
        poly->length--;
}

static __inline__
void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)
{
    long t;
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

void nmod_poly_randtest(nmod_poly_t poly, long len);

static __inline__
ulong nmod_poly_get_coeff_ui(nmod_poly_t poly, ulong j)
{
    return (j >= poly->length) ? 0 : poly->coeffs[j];
}

void nmod_poly_set_coeff_ui(nmod_poly_t poly, ulong j, ulong c);

static __inline__
void nmod_poly_set(nmod_poly_t a, nmod_poly_t b)
{
    if (a != b)
    {
        nmod_poly_fit_length(a, b->length);
        mpn_copyi(a->coeffs, b->coeffs, b->length);
        a->length = b->length;
   }
}

static __inline__
int nmod_poly_equal(nmod_poly_t a, nmod_poly_t b)
{
    if (a->length != b->length)
        return 0;

    if (a != b)
        if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
            return 0;

   return 1;
}

char * nmod_poly_to_string(nmod_poly_t poly);

int nmod_poly_from_string(char * s, nmod_poly_t poly);

static __inline__
void nmod_poly_print(nmod_poly_t a)
{
    long i;

    printf("%ld %lu", a->length, a->mod.n);

    if (a->length == 0)
        return;
    else
        printf(" ");

    for (i = 0; i < a->length; i++)
        printf(" %lu", a->coeffs[i]);
}

int nmod_poly_fread(FILE * f, nmod_poly_t poly);

static __inline__
void nmod_poly_fprint(FILE * f, nmod_poly_t poly)
{
   char * s = nmod_poly_to_string(poly);
   if (fputs(s, f) < 0) 
       printf("Error writing to file\n");
   free(s);
}

static __inline__
int nmod_poly_read(nmod_poly_t poly)
{
    return nmod_poly_fread(stdin, poly);
}

static __inline__
long nmod_poly_length(nmod_poly_t poly)
{
    return poly->length;
}

static __inline__
long nmod_poly_degree(nmod_poly_t poly)
{
    return poly->length - 1;
}

static __inline__
mp_limb_t nmod_poly_modulus(nmod_poly_t poly)
{
    return poly->mod.n;
}

static __inline__
int nmod_poly_is_zero(nmod_poly_t poly)
{
    return (poly->length == 0);
}

static __inline__
void nmod_poly_zero(nmod_poly_t res)
{
    res->length = 0;
}

static __inline__
void nmod_poly_truncate(nmod_poly_t poly, long len)
{
    if (poly->length > len)
    {
        poly->length = len;
        _nmod_poly_normalise(poly);
    }
}

void _nmod_poly_reverse(mp_ptr output, mp_srcptr input, long len, long m);

void nmod_poly_reverse(nmod_poly_t output, nmod_poly_t input, long m);

void nmod_poly_neg(nmod_poly_t res, const nmod_poly_t poly1);

void _nmod_poly_make_monic(mp_ptr output, 
                                   mp_srcptr input, long len, nmod_t mod);

void nmod_poly_make_monic(nmod_poly_t output, nmod_poly_t input);

void nmod_poly_shift_left(nmod_poly_t res, nmod_poly_t poly, long k);

void nmod_poly_shift_right(nmod_poly_t res, nmod_poly_t poly, long k);

void _nmod_poly_add(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1, 
                                                      const nmod_poly_t poly2);

void _nmod_poly_sub(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, 
                                                      const nmod_poly_t poly2);

void nmod_poly_scalar_mul(nmod_poly_t res, 
                                         const nmod_poly_t poly1, mp_limb_t c);

void _nmod_poly_mul_classical(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_mul_classical(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow_classical(mp_ptr res, mp_srcptr poly1, long len1, 
                           mp_srcptr poly2, long len2, long trunc, nmod_t mod);

void nmod_poly_mullow_classical(nmod_poly_t res, 
                 const nmod_poly_t poly1, const nmod_poly_t poly2, long trunc);

void _nmod_poly_mulhigh_classical(mp_ptr res, mp_srcptr poly1, long len1, 
                           mp_srcptr poly2, long len2, long start, nmod_t mod);

void nmod_poly_mulhigh_classical(nmod_poly_t res, 
                 const nmod_poly_t poly1, const nmod_poly_t poly2, long start);

void _nmod_poly_bit_pack(mp_ptr res, mp_srcptr poly, 
                                              long len, mp_bitcnt_t bits);

void _nmod_poly_bit_unpack(mp_ptr res, mp_srcptr mpn, 
                                       long len, mp_bitcnt_t bits, nmod_t mod);

void _nmod_poly_mul_KS(mp_ptr out, mp_srcptr in1, long len1, 
                       mp_srcptr in2, long len2, mp_bitcnt_t bits, nmod_t mod);

void nmod_poly_mul_KS(nmod_poly_t res, 
           const nmod_poly_t poly1, const nmod_poly_t poly2, mp_bitcnt_t bits);

void nmod_poly_mul(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2);

void nmod_poly_mullow_n(nmod_poly_t res, nmod_poly_t poly1, 
                                         nmod_poly_t poly2, long trunc);

void nmod_poly_mulhigh_n(nmod_poly_t res, nmod_poly_t poly1, 
                                          nmod_poly_t poly2, long n);

void _nmod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, 
                 mp_srcptr A, long A_len, mp_srcptr B, long B_len, nmod_t mod);

void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, 
                                                 nmod_poly_t A, nmod_poly_t B);

void _nmod_poly_derivative(mp_ptr x_prime, mp_srcptr x, long len, nmod_t mod);

void nmod_poly_derivative(nmod_poly_t x_prime, nmod_poly_t x);

mp_limb_t nmod_poly_evaluate(nmod_poly_t poly, mp_limb_t c);

void nmod_poly_compose(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2);

#endif






