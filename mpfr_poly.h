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

#ifndef MFPR_POLY_H
#define MPFR_POLY_H

#include <gmp.h>
#include <mpfr.h> 
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    __mpfr_struct * coeffs;
    long length;
    long alloc;
    mpfr_prec_t prec;
} mpfr_poly_struct;

/* fmpz_poly_t allows reference-like semantics for fmpz_poly_struct */
typedef mpfr_poly_struct mpfr_poly_t[1];

extern gmp_randstate_t mpfr_poly_randstate;

#define MUL_INPLACE_CUTOFF 1000

void mpfr_poly_init(mpfr_poly_t poly, mpfr_prec_t prec);

void mpfr_poly_init2(mpfr_poly_t poly, long alloc, mpfr_prec_t prec);

void mpfr_poly_realloc(mpfr_poly_t poly, long alloc);

void mpfr_poly_fit_length(mpfr_poly_t poly, long length);

void mpfr_poly_clear(mpfr_poly_t poly);

static __inline__
void _mpfr_poly_set_length(mpfr_poly_t poly, long length)
{
   poly->length = length;
}

static __inline__
void mpfr_poly_set_prec(mpfr_poly_t poly, mpfr_prec_t prec)
{
    long i;
    for (i = 0; i < poly->alloc; i++)
       mpfr_prec_round(poly->coeffs + i, prec, GMP_RNDN);
    poly->prec = prec;
}

void mpfr_poly_randinit(void);

void mpfr_poly_randclear(void);

void mpfr_poly_randtest(mpfr_poly_t poly, long length);

static __inline__
void mpfr_poly_swap(mpfr_poly_t poly1, mpfr_poly_t poly2)
{
    mpfr * tc;
    long t;
    mpfr_prec_t tp;

    tc = poly1->coeffs;
    poly1->coeffs = poly2->coeffs;
    poly2->coeffs = tc;
    
    t = poly1->length;
    poly1->length = poly2->length;
    poly2->length = t;
    
    t = poly1->alloc;
    poly1->alloc = poly2->alloc;
    poly2->alloc = t;
    
    tp = poly1->prec;
    poly1->prec = poly2->prec;
    poly2->prec = tp;
}

void _mpfr_poly_mul_classical(mpfr * res, mpfr * in1, long len1,
                              mpfr * in2, long len2, mpfr_prec_t prec);

void mpfr_poly_mul_classical(mpfr_poly_t res, mpfr_poly_t poly1, 
                                                    mpfr_poly_t poly2);

void _mpfr_poly_FHT(mpfr * coeffs, long n, mpfr_prec_t prec);

void _mpfr_poly_convolution_trans(mpfr * coeffs1, 
                             mpfr * coeffs2, long n, mpfr_prec_t prec);

void _mpfr_poly_revbin(mpfr * coeffs, long n);

void _mpfr_poly_scale(mpfr * coeffs, long n);

void _mpfr_poly_convolution_FHT(mpfr * coeffs1, 
                             mpfr * coeffs2, long n, mpfr_prec_t prec);

void mpfr_poly_mul_FHT(mpfr_poly_t res, mpfr_poly_t poly1, 
                                                    mpfr_poly_t poly2);

int _mpfr_poly_bound_newton(double * inter, double * slope, 
                              mpfr * poly, long len, mpfr_prec_t prec);

void mpfr_poly_mul(mpfr_poly_t res, mpfr_poly_t poly1, 
                                    mpfr_poly_t poly2, mpfr_prec_t fb);

#ifdef __cplusplus
}
#endif

#endif
