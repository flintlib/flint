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

    Copyright (C) 2012 Andres Goens

******************************************************************************/

#ifndef FQ_POLY_H
#define FQ_POLY_H

#undef ulong /* interferes with system includes */
#include <stdio.h>
#define ulong unsigned long
#include "flint.h"
#include "fq.h"

/*  Type definitions *********************************************************/

typedef struct
{
    fq_struct *coeffs;
    fq_ctx_t ctx;
    long alloc;
    long length;
} fq_poly_struct;

typedef fq_poly_struct fq_poly_t[1];


/*  Memory management ********************************************************/

void fq_poly_init(fq_poly_t poly, fq_ctx_t ctx);

void fq_poly_clear(fq_poly_t);

void fq_poly_init2(fq_poly_t poly, long alloc);

void fq_poly_realloc(fq_poly_t poly, long alloc);

void fq_poly_fit_length(fq_poly_t poly, long len);

void fq_poly_clear(fq_poly_t poly);

void _fq_poly_normalise(fq_poly_t poly);

void _fq_poly_set_length(fq_poly_t poly, long newlen)

/*  Polynomial parameters  ***************************************************/

static __inline__
long fq_poly_length(const fq_poly_t poly)
{
    return poly->length;
}

static __inline__
long fq_poly_degree(const fq_poly_t poly)
{
    return poly->length - 1;
}

/*  Assignment and basic manipulation  ***************************************/

void fq_poly_set(fq_poly_t poly1, const fq_poly_t poly2);

static __inline__
void fq_poly_zero(fq_poly_t poly)
{
   _fq_poly_set_length(poly, 0);
}

void fq_poly_one(fq_poly_t poly);

void fq_poly_zero_coeffs(fq_poly_t poly, long i, long j);

void fq_poly_swap(fq_poly_t poly1, fq_poly_t poly2);

void fq_poly_change_ctx(fq_poly_t poly, const fq_ctx_t ctx);

/*  Randomisation  ***********************************************************/

void fq_poly_randtest(fq_poly_t f, const fq_ctx_t ctx,
                      flint_rand_t state, long len);
void fq_poly_randtest_not_zero(fq_poly_t f, const fq_ctx_t ctx,
                               flint_rand_t state, long len);


/*  Getting and setting coefficients  ****************************************/

void fq_poly_get_coeff(fq_t x, const fq_poly_t poly, long n);

void fq_poly_set_coeff(fq_poly_t poly, long n, fq_t x);

/*  Comparison  **************************************************************/

int fq_poly_equal(const fq_poly_t poly1, const fq_poly_t poly2);

static __inline__
int fq_poly_is_zero(fq_poly_t poly)
{
    return ((poly)->length == 0);
}

static __inline__
int fq_poly_is_one(const fq_poly_t op)
{
    return (op->length) == 1 && (fq_is_one((op->coeffs)));
}

static __inline__
int fq_poly_is_unit(const fq_poly_t op)
{
    return (op->length == 1) && (!(fq_is_zero(op->coeffs)));
}

static __inline__
int fq_poly_equal_fq(const fq_poly_t poly, const fq_t c)
{
	return ((poly->length == 0) && fq_is_zero(c)) ||
        ((poly->length == 1) && fq_equal(poly->coeffs, c));
}

/*  Addition and subtraction  ************************************************/


void fq_poly_add(fq_poly_t res, const fq_poly_t poly1, 
                                                   const fq_poly_t poly2);


void fq_poly_sub(fq_poly_t res, const fq_poly_t poly1, 
                                                   const fq_poly_t poly2);

void fq_poly_neg(fq_poly_t res, const fq_poly_t poly);


#endif

/*  Scalar multiplication and division  **************************************/


void fq_poly_scalar_mul_fq(fq_poly_t poly1, 
                               const fq_poly_t poly2, const fq_t x);

void fq_poly_scalar_addmul_fq(fq_poly_t poly1, 
                                   const fq_poly_t poly2, const fq_t x);

void fq_poly_scalar_submul_fq(fq_poly_t poly1, 
                                   const fq_poly_t poly2, const fq_t x);

void fq_poly_scalar_divexact_fq(fq_poly_t poly1, 
                                    const fq_poly_t poly2, const fq_t x);
/*  Multiplication  **********************************************************/

void fq_poly_mul_classical(fq_poly_t res, 
                          const fq_poly_t poly1, const fq_poly_t poly2);

/* Squaring ******************************************************************/

void fq_poly_sqr(fq_poly_t rop, const fq_poly_t op);

/*  Powering  ****************************************************************/


void fq_poly_pow(fq_poly_t res, const fq_poly_t poly, ulong e);

/*  Shifting  ****************************************************************/

void fq_poly_shift_left(fq_poly_t res, const fq_poly_t poly, long n);

void fq_poly_shift_right(fq_poly_t res, const fq_poly_t poly, long n);

/*  Norms  *******************************************************************/

long fq_poly_hamming_wt(fq_poly_t poly);

/*  Greatest common divisor  *************************************************/


void fq_poly_gcd(fq_poly_t res, const fq_poly_t poly1, 
                                                    const fq_poly_t poly2);

/*  Euclidean division  ******************************************************/

void fq_poly_divrem_basecase(fq_poly_t Q, fq_poly_t R, 
                                   const fq_poly_t A, const fq_poly_t B);

/*  Divisibility testing  ***************************************************/

int fq_poly_divides(fq_poly_t q, const fq_poly_t a, const fq_poly_t b);

/*  Derivative  **************************************************************/

void fq_poly_derivative(fq_poly_t res, const fq_poly_t poly);

/*  Evaluation  **************************************************************/

void fq_poly_evaluate_fq(fq_t res, const fq_poly_t f, const fq_t a);

/*  Composition  *************************************************************/

void fq_poly_compose(fq_poly_t res, const fq_poly_t poly1, 
                                                      const fq_poly_t poly2);
/*  Input and output  ********************************************************/

int fq_poly_fprint(FILE * file, const fq_poly_t poly);

int fq_poly_fprint_pretty(FILE * file, 
                                       const fq_poly_t poly, const char * x);

static __inline__
int fq_poly_print(const fq_poly_t poly)
{
    return fq_poly_fprint(stdout, poly);
}

static __inline__
int fq_poly_print_pretty(const fq_poly_t poly, const char * x)
{
    return fq_poly_fprint_pretty(stdout, poly, x);
}
