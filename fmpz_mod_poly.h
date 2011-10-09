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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#ifndef FMPZ_MOD_POLY_H
#define FMPZ_MOD_POLY_H

#undef ulong /* interferes with system includes */
#include <stdio.h>
#define ulong unsigned long

#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

/*  Type definitions *********************************************************/

typedef struct
{
    fmpz * coeffs;
    long alloc;
    long length;
    fmpz p;
} fmpz_mod_poly_struct;

typedef fmpz_mod_poly_struct fmpz_mod_poly_t[1];

/*  Initialisation and memory management *************************************/

void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_t p);

void fmpz_mod_poly_init2(fmpz_mod_poly_t poly, const fmpz_t p, long alloc);

void fmpz_mod_poly_clear(fmpz_mod_poly_t poly);

void fmpz_mod_poly_realloc(fmpz_mod_poly_t poly, long alloc);

void fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, long len);

/*  Normalisation and truncation *********************************************/

void _fmpz_mod_poly_normalise(fmpz_mod_poly_t poly);

static __inline__ 
void _fmpz_mod_poly_set_length(fmpz_mod_poly_t poly, long len)
{
    if (poly->length > len)
    {
        long i;

        for (i = len; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = len;
}

static __inline__ 
void fmpz_mod_poly_truncate(fmpz_mod_poly_t poly, long len)
{
    if (poly->length > len)
    {
        long i;

        for (i = len; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = len;
        _fmpz_mod_poly_normalise(poly);
    }  
}

/*  Randomisation ************************************************************/

void
fmpz_mod_poly_randtest(fmpz_mod_poly_t f, flint_rand_t state, long len);

void
fmpz_mod_poly_randtest_not_zero(fmpz_mod_poly_t f, 
                                flint_rand_t state, long len);

/*  Conversion ***************************************************************/

void fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f, const fmpz_poly_t g);

void fmpz_mod_poly_get_fmpz_poly(fmpz_poly_t f, const fmpz_mod_poly_t g);

/*  Reduction modulo p *******************************************************/

static __inline__ 
void _fmpz_mod_poly_reduce_coeffs(fmpz_mod_poly_t poly)
{
    long i;

    for (i = 0; i < poly->length; i++)
        fmpz_mod(poly->coeffs + i, poly->coeffs + i, &(poly->p));
    _fmpz_mod_poly_normalise(poly);
}

static __inline__ 
void _fmpz_poly_reduce_coeffs(fmpz_poly_t poly, const fmpz_t p)
{
    long i;

    for (i = 0; i < poly->length; i++)
        fmpz_mod(poly->coeffs + i, poly->coeffs + i, p);
}

/*  Attributes ***************************************************************/

static __inline__ 
long fmpz_mod_poly_degree(const fmpz_mod_poly_t poly)
{
    return poly->length - 1;
}

static __inline__ 
long fmpz_mod_poly_length(const fmpz_mod_poly_t poly)
{
    return poly->length;
}

static __inline__
fmpz * fmpz_mod_poly_lead(const fmpz_mod_poly_t poly)
{
    if (poly->length)
        return poly->coeffs + (poly->length - 1);
    else
        return NULL;
}

/*  Assignment and swap ******************************************************/

void fmpz_mod_poly_set(fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

void fmpz_mod_poly_swap(fmpz_mod_poly_t poly1, fmpz_mod_poly_t poly2);

static __inline__ 
void fmpz_mod_poly_zero(fmpz_mod_poly_t poly)
{
   _fmpz_mod_poly_set_length(poly, 0);
}

void fmpz_mod_poly_zero_coeffs(fmpz_mod_poly_t poly, long i, long j);

/*  Comparison ***************************************************************/

static __inline__ 
int fmpz_mod_poly_equal(const fmpz_mod_poly_t poly1, 
                        const fmpz_mod_poly_t poly2)
{
    return fmpz_poly_equal((fmpz_poly_struct *) poly1, 
                           (fmpz_poly_struct *) poly2);
}

static __inline__ 
int fmpz_mod_poly_is_zero(const fmpz_mod_poly_t poly)
{
    return (poly->length == 0);
}

/*  Getting and setting coefficients *****************************************/

void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, long n, const fmpz_t x);

static __inline__ 
void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly, long n)
{
    if (n < poly->length)
        fmpz_set(x, poly->coeffs + n);
    else
        fmpz_zero(x);
}

/*  Shifting *****************************************************************/

void _fmpz_mod_poly_shift_left(fmpz * res, const fmpz * poly, long len, long n);

void fmpz_mod_poly_shift_left(fmpz_mod_poly_t f, 
                              const fmpz_mod_poly_t g, long n);

void _fmpz_mod_poly_shift_right(fmpz * res, const fmpz * poly, long len, long n);

void fmpz_mod_poly_shift_right(fmpz_mod_poly_t f, 
                               const fmpz_mod_poly_t g, long n);

/*  Addition and subtraction *************************************************/

void _fmpz_mod_poly_add(fmpz *res, const fmpz *poly1, long len1, 
                                   const fmpz *poly2, long len2, const fmpz_t p);

void fmpz_mod_poly_add(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

void _fmpz_mod_poly_sub(fmpz *res, const fmpz *poly1, long len1, 
                                   const fmpz *poly2, long len2, const fmpz_t p);

void fmpz_mod_poly_sub(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

void _fmpz_mod_poly_neg(fmpz *res, const fmpz *poly, long len, const fmpz_t p);

void fmpz_mod_poly_neg(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

/*  Scalar multiplication ****************************************************/

void _fmpz_mod_poly_scalar_mul_fmpz(fmpz *res, const fmpz *poly, long len, 
                                    const fmpz_t p, const fmpz_t x);

void fmpz_mod_poly_scalar_mul_fmpz(fmpz_mod_poly_t res, 
    const fmpz_mod_poly_t poly, const fmpz_t x);

/*  Multiplication ***********************************************************/

void _fmpz_mod_poly_mul(fmpz *res, const fmpz *poly1, long len1, 
                                   const fmpz *poly2, long len2, const fmpz_t p);

void fmpz_mod_poly_mul(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

void _fmpz_mod_poly_mullow(fmpz *res, const fmpz *poly1, long len1, 
                                      const fmpz *poly2, long len2, 
                                      const fmpz_t p, long n);

void fmpz_mod_poly_mullow(fmpz_mod_poly_t res, 
    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, long n);

/*  Division *****************************************************************/

void _fmpz_mod_poly_divrem_basecase(fmpz * Q, fmpz * R, 
    const fmpz * A, long lenA, const fmpz * B, long lenB, 
    const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_divrem_basecase(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

void _fmpz_mod_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, fmpz * W, 
    const fmpz * A, const fmpz * B, long lenB, 
    const fmpz_t invB, const fmpz_t p);

void _fmpz_mod_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
    const fmpz * A, long lenA, const fmpz * B, long lenB, 
    const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_divrem_divconquer(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                                     const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

static __inline__
void _fmpz_mod_poly_divrem(fmpz *Q, fmpz *R, 
                           const fmpz *A, long lenA, const fmpz *B, long lenB, 
                           const fmpz_t invB, const fmpz_t p)
{
    _fmpz_mod_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, invB, p);
}

static __inline__ 
void fmpz_mod_poly_divrem(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_divrem_divconquer(Q, R, A, B);
}

/*  Greatest common divisor **************************************************/

void fmpz_mod_poly_make_monic(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

long _fmpz_mod_poly_gcd_euclidean(fmpz *G, const fmpz *A, long lenA, 
                                           const fmpz *B, long lenB, 
                                           const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_gcd_euclidean(fmpz_mod_poly_t G, 
                                 const fmpz_mod_poly_t A,
                                 const fmpz_mod_poly_t B);

/*  Input and output *********************************************************/

static __inline__ 
int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly)
{
    return fmpz_poly_fprint(file, (fmpz_poly_struct *) poly);
}

static __inline__ 
int fmpz_mod_poly_fprint_pretty(FILE * file, 
                                const fmpz_mod_poly_t poly, const char * x)
{
    return fmpz_poly_fprint_pretty(file, (fmpz_poly_struct *) poly, x);
}

static __inline__
int fmpz_mod_poly_print(const fmpz_mod_poly_t poly)
{
    return fmpz_poly_print((fmpz_poly_struct *) poly);
}

static __inline__
int fmpz_mod_poly_print_pretty(const fmpz_mod_poly_t poly, const char * x)
{
    return fmpz_poly_print_pretty((fmpz_poly_struct *) poly, x);
}

#endif

