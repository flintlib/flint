/*=============================================================================

    fmpz_sparse.h: Sparse univariate Laurent polynomials with 
    fmpz coefficients and fmpz exponents.

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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#ifndef FMPZ_SPARSE_H
#define FMPZ_SPARSE_H

#ifdef FMPZ_SPARSE_INLINES_C
#define FMPZ_SPARSE_INLINE FLINT_DLL
#else
#define FMPZ_SPARSE_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#include <string.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_vec.h"
#include "nmod_poly.h"
#include "fmpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FMPZ_SPARSE_HEAP_XOVER (32)

/*  Type definitions *********************************************************/

typedef struct
{
    fmpz *coeffs;  /* Nonzero coefficients, in same order as the exponents. */
    fmpz *expons;  /* Term exponents, in reverse sorted order. */
    slong length;  /* Number of nonzero terms. */
    slong alloc;   /* Size that has been allocated for coeffs and expons. */
} fmpz_sparse_struct;

typedef fmpz_sparse_struct fmpz_sparse_t[1];

typedef struct
{
    fmpz q;
    fmpz * sample_points;
    fmpz * evaluations;
    slong length;
} fmpz_sparse_bp_interp_struct;

typedef fmpz_sparse_bp_interp_struct fmpz_sparse_bp_interp_t[1];

typedef struct
{
    mp_ptr shifts;
    nmod_t * cmods;
    nmod_t * emods;
    nmod_poly_struct * evaluations;
} fmpz_sparse_sp_interp_struct;

typedef fmpz_sparse_sp_interp_struct fmpz_sparse_sp_interp_t[1];

/*  Memory management ********************************************************/

FLINT_DLL void fmpz_sparse_init(fmpz_sparse_t poly);

/*  BEGINNING OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
FLINT_DLL void fmpz_sparse_init2(fmpz_sparse_t poly, slong alloc);

FLINT_DLL void fmpz_sparse_realloc(fmpz_sparse_t poly, slong alloc);
/*  ENDING OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

FLINT_DLL void fmpz_sparse_clear(fmpz_sparse_t poly);

FLINT_DLL void _fmpz_sparse_normalise(fmpz_sparse_t poly);

FLINT_DLL void _fmpz_sparse_reserve(fmpz_sparse_t poly, slong terms);

/*  Polynomial parameters  ***************************************************/

FMPZ_SPARSE_INLINE 
slong fmpz_sparse_terms(const fmpz_sparse_t poly)
{
    return poly->length;
}

FMPZ_SPARSE_INLINE 
void fmpz_sparse_degree(fmpz_t res, const fmpz_sparse_t poly)
{
    if (poly->length > 0) fmpz_set(res, poly->expons + 0);
    else fmpz_set_si(res, -1);
}

FMPZ_SPARSE_INLINE 
slong fmpz_sparse_degree_si(const fmpz_sparse_t poly)
{
    if (poly->length > 0) return fmpz_get_si(poly->expons + 0);
    else return -1;
}

FMPZ_SPARSE_INLINE 
void fmpz_sparse_lowdeg(fmpz_t res, const fmpz_sparse_t poly)
{
    if (poly->length > 0) fmpz_set(res, poly->expons + (poly->length-1));
    else fmpz_set_si(res, 1);
}

FMPZ_SPARSE_INLINE 
slong fmpz_sparse_lowdeg_si(const fmpz_sparse_t poly)
{
    if (poly->length > 0) return fmpz_get_si(poly->expons + (poly->length-1));
    else return 1;
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_is_poly(const fmpz_sparse_t poly)
{
    return poly->length == 0 || fmpz_sgn(poly->expons + (poly->length-1)) >= 0;
}

/*  Assignment and basic manipulation  ***************************************/

FLINT_DLL void fmpz_sparse_zero(fmpz_sparse_t poly);

FMPZ_SPARSE_INLINE
void fmpz_sparse_one(fmpz_sparse_t poly)
{
    fmpz_sparse_zero(poly);
    FLINT_ASSERT(poly->alloc >= 1);
    fmpz_init_set_si(poly->coeffs + 0, 1);
    fmpz_init_set_si(poly->expons + 0, 0);
    poly->length = 1;
}

FLINT_DLL void fmpz_sparse_set(fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2);

FMPZ_SPARSE_INLINE
void fmpz_sparse_set_fmpz_fmpz(fmpz_sparse_t poly, 
    const fmpz_t coeff, const fmpz_t expon)
{
    fmpz_sparse_zero(poly);
    if (!fmpz_is_zero(coeff)) {
        FLINT_ASSERT(poly->alloc >= 1);
        fmpz_init_set(poly->coeffs + 0, coeff);
        fmpz_init_set(poly->expons + 0, expon);
        poly->length = 1;
    }
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_set_fmpz_si(fmpz_sparse_t poly, 
    const fmpz_t coeff, slong expon)
{
    fmpz_sparse_zero(poly);
    if (!fmpz_is_zero(coeff)) {
        FLINT_ASSERT(poly->alloc >= 1);
        fmpz_init_set(poly->coeffs + 0, coeff);
        fmpz_init_set_si(poly->expons + 0, expon);
        poly->length = 1;
    }
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_set_si_fmpz(fmpz_sparse_t poly, 
    slong coeff, const fmpz_t expon)
{
    fmpz_sparse_zero(poly);
    if (coeff) {
        FLINT_ASSERT(poly->alloc >= 1);
        fmpz_init_set_si(poly->coeffs + 0, coeff);
        fmpz_init_set(poly->expons + 0, expon);
        poly->length = 1;
    }
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_set_si_si(fmpz_sparse_t poly, 
    slong coeff, slong expon)
{
    fmpz_sparse_zero(poly);
    if (coeff) {
        FLINT_ASSERT(poly->alloc >= 1);
        fmpz_init_set_si(poly->coeffs + 0, coeff);
        fmpz_init_set_si(poly->expons + 0, expon);
        poly->length = 1;
    }
}

/* FIXME */
FLINT_DLL void fmpz_sparse_set_fmpz_poly(fmpz_sparse_t poly1, 
    const fmpz_poly_t poly2);

FLINT_DLL void fmpz_sparse_get_fmpz_poly(fmpz_poly_t out, 
    const fmpz_sparse_t in);

FLINT_DLL int fmpz_sparse_set_str(fmpz_sparse_t poly, const char * str);

FLINT_DLL char * fmpz_sparse_get_str(const fmpz_sparse_t poly);

/* FIXME */
FLINT_DLL char * fmpz_sparse_get_str_pretty(const fmpz_sparse_t poly, const char * x);

FMPZ_SPARSE_INLINE
void fmpz_sparse_swap(fmpz_sparse_t poly1, fmpz_sparse_t poly2)
{
    FLINT_GENERIC_SWAP(fmpz *, poly1->coeffs, poly2->coeffs);
    FLINT_GENERIC_SWAP(fmpz *, poly1->expons, poly2->expons);
    FLINT_GENERIC_SWAP(slong, poly1->length, poly2->length);
    FLINT_GENERIC_SWAP(slong, poly1->alloc, poly2->alloc);
}

/* FIXME */
FLINT_DLL void fmpz_sparse_reverse(fmpz_sparse_t res, 
    const fmpz_sparse_t poly, slong n);

/* FIXME */
FLINT_DLL void fmpz_sparse_truncate(fmpz_sparse_t poly, const fmpz_t deg);

/* FIXME */
FLINT_DLL void fmpz_sparse_set_trunc(fmpz_sparse_t res, 
    const fmpz_sparse_t poly, const fmpz_t deg);

/* FIXME */
FLINT_DLL void fmpz_sparse_set_trunc_fmpz_poly(fmpz_poly_t res, 
    const fmpz_sparse_t poly, slong deg);

/*  Randomisation  ***********************************************************/

FLINT_DLL void fmpz_sparse_randtest(fmpz_sparse_t res, flint_rand_t state, 
     slong terms, const fmpz_t degree, mp_bitcnt_t bits);

/* FIXME */
FLINT_DLL void fmpz_sparse_randtest_unsigned(fmpz_sparse_t res, 
    flint_rand_t state, slong terms, const fmpz_t degree, mp_bitcnt_t bits);

/* FIXME */
FLINT_DLL void fmpz_sparse_randtest_not_zero(fmpz_sparse_t res, 
    flint_rand_t state, slong terms, const fmpz_t degree, mp_bitcnt_t bits);

/*  Getting and setting coefficients  ****************************************/

/* FIXME */
FLINT_DLL slong fmpz_sparse_get_coeff_si_si(const fmpz_sparse_t poly, 
    slong e);

/* FIXME */
FLINT_DLL slong fmpz_sparse_get_coeff_si_fmpz(const fmpz_sparse_t poly, 
    const fmpz_t e);

/* FIXME */
FLINT_DLL void fmpz_sparse_get_coeff_fmpz_si(fmpz_t res,
    const fmpz_sparse_t poly, slong e);

FLINT_DLL void fmpz_sparse_get_coeff(fmpz_t res,
    const fmpz_sparse_t poly, const fmpz_t e);

/* FIXME */
FLINT_DLL void fmpz_sparse_set_coeff_si_si(fmpz_sparse_t poly, 
    slong c, slong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_set_coeff_si_fmpz(fmpz_sparse_t poly, 
    slong c, const fmpz_t e);

/* FIXME */
FLINT_DLL void fmpz_sparse_set_coeff_fmpz_si(fmpz_sparse_t poly, 
    const fmpz_t c, slong e);

FLINT_DLL void fmpz_sparse_set_coeff(fmpz_sparse_t poly, 
    const fmpz_t c, const fmpz_t e);

/* FIXME */
FLINT_DLL fmpz* fmpz_sparse_get_coeff_ptr(fmpz_sparse_t poly, const fmpz_t e);

/* FIXME */
FLINT_DLL fmpz* fmpz_sparse_get_coeff_ptr_si(fmpz_sparse_t poly, slong e);

FMPZ_SPARSE_INLINE 
void fmpz_sparse_get_term(fmpz_t coeff, fmpz_t expon, 
    const fmpz_sparse_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    fmpz_set(coeff, poly->coeffs + i);
    fmpz_set(expon, poly->expons + i);
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_get_term_coeff(fmpz_t res, const fmpz_sparse_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    fmpz_set(res, poly->coeffs + i);
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_get_term_expon(fmpz_t res, const fmpz_sparse_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    fmpz_set(res, poly->expons + i);
}

FMPZ_SPARSE_INLINE
fmpz* fmpz_sparse_get_term_coeff_ptr(const fmpz_sparse_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    return poly->coeffs + i;
}

FMPZ_SPARSE_INLINE
fmpz* fmpz_sparse_get_term_expon_ptr(const fmpz_sparse_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    return poly->expons + i;
}

FMPZ_SPARSE_INLINE
slong fmpz_sparse_get_term_coeff_si(const fmpz_sparse_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    return fmpz_get_si(poly->coeffs + i);
}

FMPZ_SPARSE_INLINE
slong fmpz_sparse_get_term_expon_si(const fmpz_sparse_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    return fmpz_get_si(poly->expons + i);
}

FLINT_DLL slong _fmpz_sparse_index(const fmpz_sparse_t poly, const fmpz_t e);

/*  Comparison  **************************************************************/

FMPZ_SPARSE_INLINE 
int fmpz_sparse_equal(const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
    return
        poly1->length == poly2->length &&
        _fmpz_vec_equal(poly1->expons, poly2->expons, poly1->length) &&
        _fmpz_vec_equal(poly1->coeffs, poly2->coeffs, poly1->length);
}

/* FIXME */
FLINT_DLL int fmpz_sparse_equal_trunc(const fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2, const fmpz_t n);

FMPZ_SPARSE_INLINE 
int fmpz_sparse_is_zero(const fmpz_sparse_t poly)
{
    return poly->length == 0;
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_is_unit(const fmpz_sparse_t poly)
{
    return poly->length == 1 && 
        fmpz_is_zero(poly->expons + 0) && 
        fmpz_is_pm1(poly->coeffs + 0);
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_is_term(const fmpz_sparse_t poly, 
    const fmpz_t c, const fmpz_t e)
{
    return (poly->length == 0 && fmpz_is_zero(c)) ||
        (poly->length == 1 && 
         fmpz_equal(poly->coeffs+0, c) && 
         fmpz_equal(poly->expons+0, e)
        );
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_is_term_fmpz_si(const fmpz_sparse_t poly, 
    const fmpz_t c, slong e)
{
    return (poly->length == 0 && fmpz_is_zero(c)) ||
        (poly->length == 1 && 
         fmpz_equal(poly->coeffs+0, c) && 
         fmpz_equal_si(poly->expons+0, e)
        );
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_is_term_si_fmpz(const fmpz_sparse_t poly, 
    slong c, const fmpz_t e)
{
    return (poly->length == 0 && (c == 0)) ||
        (poly->length == 1 && 
         fmpz_equal_si(poly->coeffs+0, c) && 
         fmpz_equal(poly->expons+0, e)
        );
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_is_term_si_si(const fmpz_sparse_t poly, slong c, slong e)
{
    return (poly->length == 0 && (c == 0)) ||
        (poly->length == 1 && 
         fmpz_equal_si(poly->coeffs+0, c) && 
         fmpz_equal_si(poly->expons+0, e)
        );
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_equal_fmpz(const fmpz_sparse_t poly, const fmpz_t c)
{
    return fmpz_sparse_is_term_fmpz_si(poly, c, 0);
}

/* FIXME */
FLINT_DLL int fmpz_sparse_equal_fmpz_poly(const fmpz_sparse_t spoly, 
    const fmpz_poly_t dpoly);

/*  Addition and subtraction  ************************************************/

/*  BEGINNING OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
FLINT_DLL void _fmpz_sparse_new_add(fmpz * res_c, fmpz * res_e, slong * res_len, 
    const fmpz * poly1_c, const fmpz * poly1_e, slong len1, const fmpz * poly2_c, 
    const fmpz * poly2_e, slong len2);

FLINT_DLL void fmpz_sparse_new_add(fmpz_sparse_t res,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2);

FLINT_DLL void _fmpz_sparse_new_sub(fmpz * res_c, fmpz * res_e, slong * res_len,
        const fmpz * poly1_c, const fmpz * poly1_e, slong len1, const fmpz * poly2_c,
            const fmpz * poly2_e, slong len2);

FLINT_DLL void fmpz_sparse_new_sub(fmpz_sparse_t res,
        const fmpz_sparse_t poly1, const fmpz_sparse_t poly2);
/*  END OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

FLINT_DLL void fmpz_sparse_add(fmpz_sparse_t res,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2);

FLINT_DLL void fmpz_sparse_sub(fmpz_sparse_t res,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2);

FLINT_DLL void fmpz_sparse_neg(fmpz_sparse_t res, const fmpz_sparse_t poly);

/*  Scalar multiplication and division  **************************************/

/*  BEGINNING OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
FLINT_DLL void fmpz_sparse_scalar_mul_ui(fmpz_sparse_t res,
        const fmpz_sparse_t poly, ulong c);

FLINT_DLL void fmpz_sparse_scalar_mul_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong c);

FLINT_DLL void fmpz_sparse_scalar_mul(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c);
/*  END OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_addmul(fmpz_sparse_t poly1,
    const fmpz_sparse_t poly2, const fmpz_t c);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_submul(fmpz_sparse_t poly1,
    const fmpz_sparse_t poly2, const fmpz_t c);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_fdiv_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong c);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_fdiv(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_tdiv_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong c);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_tdiv(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_mul_2exp(fmpz_sparse_t res,
    const fmpz_sparse_t poly, ulong exp);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_fdiv_2exp(fmpz_sparse_t res,
    const fmpz_sparse_t poly, ulong exp);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_tdiv_2exp(fmpz_sparse_t res,
    const fmpz_sparse_t poly, ulong exp);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_mod(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c);

/* FIXME */
FLINT_DLL void fmpz_sparse_scalar_smod(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c);

/*  Bit packing  *************************************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_bit_pack(fmpz_t res,
    const fmpz_sparse_t poly, mp_bitcnt_t bit_size);

/* FIXME */
FLINT_DLL void fmpz_sparse_bit_unpack(fmpz_sparse_t res,
    const fmpz_t f, mp_bitcnt_t bit_size);

/*  Multiplication  **********************************************************/

/*  BEGINNING OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
FLINT_DLL void fmpz_sparse_new_mul_classical(fmpz_sparse_t res,
        const fmpz_sparse_t poly1, const fmpz_sparse_t poly2);
/*  END OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


FLINT_DLL void fmpz_sparse_mul_classical(fmpz_sparse_t res,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2);

FMPZ_SPARSE_INLINE
void fmpz_sparse_mul_heaps(fmpz_sparse_t res,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
    /* FIXME this is just a placeholder! */
    fmpz_sparse_mul_classical(res, poly1, poly2);
}

/* FIXME */
FLINT_DLL void fmpz_sparse_mul_interp(fmpz_sparse_t res,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2);

FMPZ_SPARSE_INLINE 
void fmpz_sparse_mul(fmpz_sparse_t res,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
    if (poly1->length + poly2->length < FMPZ_SPARSE_HEAP_XOVER)
    {
        fmpz_sparse_mul_classical(res, poly1, poly2);
    }
    else fmpz_sparse_mul_heaps(res, poly1, poly2);
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_sqr(fmpz_sparse_t res, const fmpz_sparse_t poly)
{
    fmpz_sparse_mul(res, poly, poly);
}

/*  Powering  ****************************************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_pow_recurrence(fmpz_sparse_t res,
    const fmpz_sparse_t poly, ulong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_pow_binomial(fmpz_sparse_t res,
    const fmpz_sparse_t poly, ulong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_pow_binexp(fmpz_sparse_t res,
    const fmpz_sparse_t poly, ulong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_pow(fmpz_sparse_t res,
    const fmpz_sparse_t poly, ulong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_pow_trunc(fmpz_sparse_t res,
    const fmpz_sparse_t poly, ulong e, const fmpz_t n);

/*  Shifting  ****************************************************************/

/*  BEGINNING OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
FMPZ_SPARSE_INLINE
void fmpz_sparse_shift_left(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t n)
{
    int i;
    for (i=0; i<poly->length; ++i)
    {
        fmpz_add(res->expons+i, poly->expons+i, n);
    }
}
/*  END OF WHITMAN'S WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* FIXME */
FLINT_DLL void fmpz_sparse_shift_left_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong n);

/* FIXME */
FLINT_DLL void fmpz_sparse_shift_right(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t n, int trunc_poly);

/* FIXME */
FLINT_DLL void fmpz_sparse_shift_right_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong n, int trunc_poly);

FLINT_DLL void _fmpz_sparse_vec_shift(fmpz_sparse_t poly, 
    slong start, slong end, slong dist);

FLINT_DLL void _fmpz_sparse_append_si(fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2, slong c, ulong e);

FLINT_DLL void _fmpz_sparse_append(fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2, const fmpz_t c, const fmpz_t e);

/*  Monomial multiplication and division *************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_mon_mul_si_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong c, slong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_mon_mul_si_fmpz(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong c, const fmpz_t e);

/* FIXME */
FLINT_DLL void fmpz_sparse_mon_mul_fmpz_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c, slong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_mon_mul_fmpz_fmpz(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c, const fmpz_t e);

/* FIXME */
FLINT_DLL void fmpz_sparse_mon_fdiv_si_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong c, slong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_mon_fdiv_si_fmpz(fmpz_sparse_t res,
    const fmpz_sparse_t poly, slong c, const fmpz_t e);

/* FIXME */
FLINT_DLL void fmpz_sparse_mon_fdiv_fmpz_si(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c, slong e);

/* FIXME */
FLINT_DLL void fmpz_sparse_mon_fdiv_fmpz_fmpz(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t c, const fmpz_t e);

/*  Norms  *******************************************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_2norm(fmpz_t res, const fmpz_sparse_t poly);

FMPZ_SPARSE_INLINE
ulong fmpz_sparse_max_limbs(const fmpz_sparse_t poly)
{
    return _fmpz_vec_max_limbs(poly->coeffs, poly->length);
}

FMPZ_SPARSE_INLINE
slong fmpz_sparse_max_bits(const fmpz_sparse_t poly)
{
    return _fmpz_vec_max_bits(poly->coeffs, poly->length);
}

FMPZ_SPARSE_INLINE
slong fmpz_sparse_max_ebits(const fmpz_sparse_t poly)
{
    fmpz * lead = poly->expons + 0;
    fmpz * trail = poly->expons + (poly->length-1);
    if (poly->length == 0) return 0;
    else if (poly->length == 1 || fmpz_cmpabs(lead, trail) >= 0)
    {
        return fmpz_bits(lead);
    }
    else return fmpz_bits(trail);
}

FMPZ_SPARSE_INLINE
ulong fmpz_sparse_max_elimbs(const fmpz_sparse_t poly)
{
    slong bc = fmpz_sparse_max_ebits(poly);
    if (bc >= 0) return (bc + FLINT_BITS - 1) / FLINT_BITS;
    else return (bc - FLINT_BITS + 1) / FLINT_BITS;
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_height(fmpz_t res, const fmpz_sparse_t poly)
{
    _fmpz_vec_height(res, poly->coeffs, poly->length);
}

/*  Euclidean division  ******************************************************/

FMPZ_SPARSE_INLINE
void fmpz_sparse_divrem(fmpz_sparse_t Q, fmpz_sparse_t R,
    const fmpz_sparse_t A, const fmpz_sparse_t B)
{
    FLINT_ASSERT(false);
    /* FIXME this is just a placeholder! */
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_div(fmpz_sparse_t Q, 
    const fmpz_sparse_t A, const fmpz_sparse_t B)
{
    /* TODO better efficiency */
    fmpz_sparse_t temp;
    fmpz_sparse_init(temp);
    fmpz_sparse_divrem(Q, temp, A, B);
    fmpz_sparse_clear(temp);
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_rem(fmpz_sparse_t R, fmpz_sparse_t A, fmpz_sparse_t B)
{
    /* TODO better efficiency */
    fmpz_sparse_t temp;
    fmpz_sparse_init(temp);
    fmpz_sparse_divrem(temp, R, A, B);
    fmpz_sparse_clear(temp);
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_divrem_dense(fmpz_sparse_t Q, fmpz_poly_t R,
    const fmpz_sparse_t A, const fmpz_poly_t B)
{
    /* FIXME this is just a placeholder! */
    FLINT_ASSERT(false);
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_div_dense(fmpz_sparse_t Q, fmpz_sparse_t A, fmpz_poly_t B)
{
    /* TODO better efficiency */
    fmpz_poly_t temp;
    fmpz_poly_init(temp);
    fmpz_sparse_divrem_dense(Q, temp, A, B);
    fmpz_poly_clear(temp);
}

/* FIXME */
FLINT_DLL void fmpz_sparse_rem_dense(fmpz_poly_t R, fmpz_sparse_t A, fmpz_poly_t B);

/* FIXME */
FLINT_DLL void fmpz_sparse_rem_cyc(fmpz_sparse_t res,
    const fmpz_sparse_t poly, const fmpz_t e);

FLINT_DLL void fmpz_sparse_rem_cyc_dense(fmpz_poly_t res,
    const fmpz_sparse_t poly, ulong e);

FLINT_DLL void fmpz_sparse_rem_cyc_nmod(nmod_poly_t res,
    const fmpz_sparse_t poly, ulong e, ulong q);

FLINT_DLL void fmpz_sparse_rem_cyc_mod_diverse(nmod_poly_t res,
    const fmpz_sparse_t poly, ulong a, ulong e, ulong q);

/*  Greatest common divisor  *************************************************/

FMPZ_SPARSE_INLINE
void fmpz_sparse_xgcd(fmpz_sparse_t r, 
    fmpz_sparse_t s, fmpz_sparse_t t,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
    /* FIXME this is just a placeholder! */
    FLINT_ASSERT(false);
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_gcd(fmpz_sparse_t res, 
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
    /* TODO better efficiency */
    fmpz_sparse_t s, t;
    fmpz_sparse_init(s);
    fmpz_sparse_init(t);
    fmpz_sparse_xgcd(res, s, t, poly1, poly2);
    fmpz_sparse_clear(s);
    fmpz_sparse_clear(t);
}

FMPZ_SPARSE_INLINE 
void fmpz_sparse_lcm(fmpz_sparse_t res, 
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
    /* TODO better efficiency */
    fmpz_sparse_gcd(res, poly1, poly2);
    fmpz_sparse_div(res, poly1, res);
    fmpz_sparse_mul(res, res, poly2);
}

/*  Gaussian content  ********************************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_content(fmpz_t res, const fmpz_sparse_t poly);

/* FIXME */
FLINT_DLL void fmpz_sparse_primitive_part(fmpz_sparse_t res, const fmpz_sparse_t poly);

/*  Sparse interpolation ****************************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_bp_interp_init(fmpz_sparse_bp_interp_t res,
    slong terms, const fmpz_t height, const fmpz_t degree);

/* FIXME */
FLINT_DLL void fmpz_sparse_bp_interp_clear(fmpz_sparse_bp_interp_t res);

/* FIXME */
FLINT_DLL void fmpz_sparse_bp_interp_eval(fmpz_sparse_bp_interp_t res,
    const fmpz_sparse_t poly);

/* FIXME */
FLINT_DLL void fmpz_sparse_bp_interp_mul(fmpz_sparse_bp_interp_t res,
    const fmpz_sparse_t poly);

/* FIXME */
FLINT_DLL void fmpz_sparse_bp_interp_add(fmpz_sparse_bp_interp_t res,
    const fmpz_t c, const fmpz_sparse_t poly);

/* FIXME */
FLINT_DLL void fmpz_sparse_bp_interp_pow(fmpz_sparse_bp_interp_t res, ulong pow);

/* FIXME */
FLINT_DLL void fmpz_sparse_bp_interp(fmpz_sparse_t res,
    const fmpz_sparse_bp_interp_t evals);

/* FIXME */
FLINT_DLL void fmpz_sparse_sp_interp_init(fmpz_sparse_sp_interp_t res,
    slong terms, const fmpz_t height, const fmpz_t degree);

/* FIXME */
FLINT_DLL void fmpz_sparse_sp_interp_clear(fmpz_sparse_sp_interp_t res);

/* FIXME */
FLINT_DLL void fmpz_sparse_sp_interp_eval(fmpz_sparse_sp_interp_t res,
    const fmpz_sparse_t poly);

/* FIXME */
FLINT_DLL void fmpz_sparse_sp_interp_mul(fmpz_sparse_sp_interp_t res,
    const fmpz_sparse_t poly);

/* FIXME */
FLINT_DLL void fmpz_sparse_sp_interp_add(fmpz_sparse_sp_interp_t res,
    const fmpz_t c, const fmpz_sparse_t poly);

/* FIXME */
FLINT_DLL void fmpz_sparse_sp_interp_pow(fmpz_sparse_sp_interp_t res, ulong pow);

/* FIXME */
FLINT_DLL void fmpz_sparse_sp_interp(fmpz_sparse_t res,
    const fmpz_sparse_sp_interp_t evals);

/*  Divisibility testing  ***************************************************/

FMPZ_SPARSE_INLINE
int fmpz_sparse_divides(fmpz_sparse_t q, 
    const fmpz_sparse_t a, const fmpz_sparse_t b)
{
    /* TODO better efficiency */
    int res;
    fmpz_sparse_t temp;
    fmpz_sparse_init(temp);
    fmpz_sparse_divrem(q, temp, a, b);
    res = fmpz_sparse_is_zero(temp);
    fmpz_sparse_clear(temp);
    return res;
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_divides_dense(fmpz_sparse_t q, 
    const fmpz_sparse_t a, const fmpz_poly_t b)
{
    /* TODO better efficiency */
    int res;
    fmpz_poly_t temp;
    fmpz_poly_init(temp);
    fmpz_sparse_divrem_dense(q, temp, a, b);
    res = fmpz_poly_is_zero(temp);
    fmpz_poly_clear(temp);
    return res;
}

/*  Derivative  **************************************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_derivative(fmpz_sparse_t res, const fmpz_sparse_t poly);

/*  Evaluation  **************************************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_evaluate(fmpz_t res,
    const fmpz_sparse_t f, const fmpz_t a);

/* FIXME */
FLINT_DLL ulong fmpz_sparse_evaluate_mod(const fmpz_sparse_t poly, 
    ulong a, ulong m);

/*  Composition  *************************************************************/

/* FIXME */
FLINT_DLL void fmpz_sparse_compose(fmpz_sparse_t res,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2);

/* FIXME */
FLINT_DLL void fmpz_sparse_compose_dense(fmpz_sparse_t res,
    const fmpz_poly_t poly1, const fmpz_sparse_t poly2);

/*  Input and output  ********************************************************/

FLINT_DLL int fmpz_sparse_fprint(FILE * file, const fmpz_sparse_t poly);

FLINT_DLL int fmpz_sparse_fprint_pretty(FILE * file, 
    const fmpz_sparse_t poly, const char *x);

FMPZ_SPARSE_INLINE
int fmpz_sparse_print(const fmpz_sparse_t poly)
{
    return fmpz_sparse_fprint(stdout, poly);
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_print_pretty(const fmpz_sparse_t poly, const char *x)
{
    return fmpz_sparse_fprint_pretty(stdout, poly, x);
}

FLINT_DLL int fmpz_sparse_fread(FILE * file, fmpz_sparse_t poly);

/* FIXME */
FMPZ_SPARSE_INLINE
int fmpz_sparse_fread_pretty(FILE * file,
    fmpz_sparse_t poly, char ** x)
{
    /* FIXME this is just a placeholder! */
    return -1;
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_read(fmpz_sparse_t poly)
{
    return fmpz_sparse_fread(stdin, poly);
}

FMPZ_SPARSE_INLINE
int fmpz_sparse_read_pretty(fmpz_sparse_t poly, char ** x)
{
    return fmpz_sparse_fread_pretty(stdin, poly, x);
}

FMPZ_SPARSE_INLINE
void fmpz_sparse_debug(const fmpz_sparse_t poly)
{
    flint_printf("(alloc = %wd, length = %wd,\n  coeffs = ", poly->alloc, poly->length);
    if (poly->coeffs) _fmpz_vec_print(poly->coeffs, poly->alloc);
    else flint_printf("NULL");
    flint_printf("\n  expons = ");
    if (poly->expons) _fmpz_vec_print(poly->expons, poly->alloc);
    else flint_printf("NULL");
    flint_printf("\n");
    fflush(stdout);
}

/*  OTHER  ******************************************************************/

/* TODO: CRT would be nice, but we would need a fmpz_sparse_nmod type first. */

/* TODO: Finding all integer roots of an fmpz_sparse polynomial */

/* TODO: Computing low-degree and cyclotomic factors of an fmpz_sparse polynomial */

#ifdef __cplusplus
}
#endif

#endif /* FMPZ_SPARSE_H */
