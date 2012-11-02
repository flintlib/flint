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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#ifndef GENERICS_H
#define GENERICS_H

#include <alloca.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_vec.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif

struct ring_struct;

typedef enum
{
    TYPE_FMPZ,
    TYPE_LIMB,
    TYPE_MOD,
    TYPE_FRAC,
    TYPE_POLY,
    TYPE_MAT,
    TYPE_COMPLEX
}
ring_type;

typedef void * elem_ptr;
typedef const void * elem_srcptr;

/* todo: make this a union */
typedef struct
{
    ring_type type;
    long size;
    void * parent;
    nmod_t nmod;
    void * modulus;
    void * numer;
    void * denom;
    long denom_offset;
}
ring_struct;

#define INDEX(_ptr, _i, _size) ((char *)(_ptr) + ((_i) * (_size)))
#define SRC_INDEX(_ptr, _i, _size) ((const char *)(_ptr) + ((_i) * (_size)))

typedef struct
{
    elem_ptr coeffs;
    long alloc;
    long length;
} elem_poly_struct;

typedef struct
{
    elem_ptr * rows;
    elem_ptr entries;
    long r;
    long c;
} elem_mat_struct;

typedef elem_poly_struct elem_poly_t[1];
typedef elem_mat_struct elem_mat_t[1];

typedef struct
{
    ring_struct * ring;
    elem_ptr elem;
} gen_struct;


typedef ring_struct ring_t[1];
typedef gen_struct gen_t[1];

#define RING_PARENT(ring) ((ring_struct *) ((ring)->parent))
#define RING_MODULUS(ring) ((elem_ptr) ((ring)->modulus))

#define RING_NUMER(ring) ((ring_struct *) ((ring)->numer))
#define RING_DENOM(ring) ((ring_struct *) ((ring)->denom))

#define NUMER(elem, ring) ((elem_ptr) (elem))
#define DENOM(elem, ring) ((elem_ptr) (((char *) (elem)) + (ring)->denom_offset))

#define REALPART(elem, ring) ((elem_ptr) (elem))
#define IMAGPART(elem, ring) ((elem_ptr) (((char *) (elem)) + (((ring)->size) >> 1)))

#define NOT_IMPLEMENTED(opname, ring) do { \
    printf("operation %s not implemented for ring ", opname); \
    ring_print(ring); printf("\n"); \
    abort(); \
  } while (0)

static __inline__ void
gen_swap(gen_t res, gen_t src)
{
    gen_struct t = *src;
    *src = *res;
    *res = t;
}

static __inline__ void
elem_poly_swap(elem_poly_struct * op1, elem_poly_struct * op2)
{
    elem_poly_struct t = *op1;
    *op1 = *op2;
    *op2 = t;
}

#define ELEM_TMP_MAX 64

#define ELEM_TMP_INIT(_x, _ring) \
    if ((_ring)->size > ELEM_TMP_MAX) \
        (_x) = flint_malloc((_ring)->size); \
    else \
        (_x) = alloca((_ring)->size); \
    elem_init((_x), (_ring));

#define ELEM_TMP_CLEAR(_x, _ring) \
    elem_clear((_x), (_ring)); \
    if ((_ring)->size > ELEM_TMP_MAX) \
        flint_free((_x));

void ring_init_fmpz(ring_t ring);
void ring_init_limb(ring_t ring);
void ring_init_mod(ring_t ring, const ring_t elem_ring, elem_srcptr modulus);
void ring_init_frac(ring_t ring, const ring_t numer_ring, const ring_t denom_ring);
void ring_init_poly(ring_t ring, const ring_t elem_ring);
void ring_init_mat(ring_t ring, const ring_t elem_ring);
void ring_init_complex(ring_t ring, const ring_t real_ring);
int ring_init_randtest(ring_t * R, flint_rand_t state, int maxdepth);
void ring_clear(ring_t ring);
void ring_print(const ring_t ring);

/* elem pointer versions */

void elem_init(elem_ptr elem, const ring_t ring);
void elem_clear(elem_ptr elem, const ring_t ring);
void elem_zero(elem_ptr x, const ring_t ring);
int elem_is_zero(elem_srcptr x, const ring_t ring);
void elem_one(elem_ptr x, const ring_t ring);
int elem_is_one(elem_srcptr x, const ring_t ring);
int elem_equal(elem_srcptr op1, elem_srcptr op2, const ring_t ring);
void elem_print(elem_srcptr elem, const ring_t ring);
void elem_swap(elem_ptr res, elem_ptr src, const ring_t ring);
void elem_set(elem_ptr res, elem_srcptr src, const ring_t ring);
void elem_set_si(elem_ptr elem, long v, const ring_t ring);
void elem_set_ui(elem_ptr elem, ulong v, const ring_t ring);
void elem_randtest(elem_ptr res, flint_rand_t state, const long * size, const ring_t ring);
void elem_randtest_not_zero(elem_ptr res, flint_rand_t state, const long * size, const ring_t ring);
void elem_neg(elem_ptr res, elem_srcptr src, const ring_t ring);
void elem_add(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring);
void elem_sub(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring);
void elem_mul(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring);
void elem_divexact(elem_ptr Q, elem_srcptr A, elem_srcptr B, const ring_t ring);
void elem_divrem(elem_ptr Q, elem_ptr R, elem_srcptr A, elem_srcptr B, const ring_t ring);
void elem_pow_ui(elem_ptr res, elem_srcptr op, ulong exp, const ring_t ring);

void elem_gcd(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring);

/* versions allowing the second operand to be in a parent ring
   (note: for gcd, the *output* is also in the parent ring) */

void elem_mul_parent(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring, const ring_t parent);
void elem_divexact_parent(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring, const ring_t parent);
void elem_gcd_parent(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring, const ring_t parent);

/* wrapped versions */

void gen_init(gen_t x, const ring_t ring);
void gen_clear(gen_t x);
void gen_zero(gen_t x);
int gen_is_zero(const gen_t x);
void gen_one(gen_t x);
int gen_is_one(const gen_t x);
int gen_equal(const gen_t op1, const gen_t op2);
void gen_print(gen_t x);
void gen_set(gen_t y, const gen_t x);
void gen_set_si(gen_t x, long v);
void gen_randtest(gen_t res, flint_rand_t state, const long * size);
void gen_randtest_not_zero(gen_t res, flint_rand_t state, const long * size);
void gen_neg(gen_t y, const gen_t x);
void gen_add(gen_t z, const gen_t x, const gen_t y);
void gen_sub(gen_t z, const gen_t x, const gen_t y);
void gen_mul(gen_t z, const gen_t x, const gen_t y);
void gen_divexact(gen_t Q, const gen_t A, const gen_t B);
void gen_divrem(gen_t Q, gen_t R, const gen_t A, const gen_t B);
void gen_pow_ui(gen_t z, const gen_t x, ulong exp);

void gen_pseudo_divrem(gen_t Q, gen_t R, ulong * d, const gen_t A, const gen_t B);


/* Fraction functions */

void elem_content_recursive(elem_ptr cont, elem_srcptr obj, const ring_t cont_ring, const ring_t obj_ring);
void elem_div_content_recursive(elem_ptr obj, elem_srcptr cont, const ring_t cont_ring, const ring_t obj_ring);
void _elem_frac_canonicalise(elem_ptr num, elem_ptr den, const ring_t num_ring, const ring_t den_ring);
void elem_frac_canonicalise(elem_srcptr x, const ring_t ring);

void elem_frac_add(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring);
void elem_frac_sub(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring);
void elem_frac_mul(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring);
void elem_frac_div(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring);

/* Vector functions */

elem_ptr _elem_vec_init(long len, const ring_t ring);
void _elem_vec_clear(elem_ptr vec, long len, const ring_t ring);
elem_ptr _elem_vec_realloc(elem_ptr vec, long old_len, long new_len, const ring_t ring);

void _elem_vec_zero(elem_ptr res, long len, const ring_t ring);
void _elem_vec_set(elem_ptr res, elem_srcptr src, long len, const ring_t ring);
void _elem_vec_neg(elem_ptr res, elem_srcptr src, long len, const ring_t ring);
void _elem_vec_scalar_mul(elem_ptr res, elem_srcptr vec, long len, elem_srcptr c, const ring_t ring);
void _elem_vec_scalar_addmul(elem_ptr res, elem_srcptr vec, long len, elem_srcptr c, const ring_t ring);
void _elem_vec_scalar_submul(elem_ptr res, elem_srcptr vec, long len, elem_srcptr c, const ring_t ring);
void _elem_vec_scalar_divexact(elem_ptr res, elem_srcptr vec, long len, elem_srcptr c, const ring_t ring);

/* Polynomial functions */

static __inline__ void
elem_poly_fit_length(elem_poly_t poly, long len, const ring_t poly_ring)
{
    long alloc = poly->alloc;

    if (len > alloc)
    {
        if (len < 2 * alloc)
            len = 2 * alloc;

        poly->coeffs = _elem_vec_realloc(poly->coeffs, alloc, len, RING_PARENT(poly_ring));
        poly->alloc = len;
    }
}

static __inline__ void
elem_poly_normalise(elem_poly_t poly, const ring_t poly_ring)
{
    long i, size = RING_PARENT(poly_ring)->size;
    elem_ptr ptr = poly->coeffs;

    i = poly->length - 1;

    for (i = poly->length - 1;
        (i >= 0) && elem_is_zero(INDEX(ptr, i, size), poly_ring->parent); i--);

    poly->length = i + 1;

/*
    if (poly->length < poly->alloc)
    {
        poly->coeffs = _elem_vec_realloc(poly->coeffs, poly->alloc, poly->length, poly_ring->parent);
        poly->alloc = poly->length;
    }
*/
}

static __inline__ void
elem_poly_set_length(elem_poly_t poly, long len, const ring_t poly_ring)
{
    poly->length = len;

/*
    if (poly->length < poly->alloc)
    {
        poly->coeffs = _elem_vec_realloc(poly->coeffs, poly->alloc, poly->length, poly_ring->parent);
        poly->alloc = poly->length;
    }
*/
}

void _elem_poly_print(elem_srcptr poly, long len, const ring_t ring);
void elem_poly_set_coeff_si(elem_poly_t elem, long index, long value, const ring_t ring);

void _elem_poly_add(elem_ptr res, elem_srcptr poly1, long len1, elem_srcptr poly2, long len2, const ring_t ring);
void elem_poly_add(elem_poly_t res, const elem_poly_t op1, const elem_poly_t op2, const ring_t ring);

void _elem_poly_sub(elem_ptr res, elem_srcptr poly1, long len1, elem_srcptr poly2, long len2, const ring_t ring);
void elem_poly_sub(elem_poly_t res, const elem_poly_t op1, const elem_poly_t op2, const ring_t ring);

void _elem_poly_mul(elem_ptr res, elem_srcptr poly1, long len1, elem_srcptr poly2, long len2, const ring_t ring);
void elem_poly_mul(elem_poly_t res, const elem_poly_t op1, const elem_poly_t op2, const ring_t ring);

void _elem_poly_divrem(elem_ptr Q, elem_ptr R, elem_srcptr A, long lenA, elem_srcptr B, long lenB, const ring_t ring);
void elem_poly_divrem(elem_poly_t Q, elem_poly_t R, const elem_poly_t A, const elem_poly_t B, const ring_t ring);

void _elem_poly_pseudo_divrem(elem_ptr Q, elem_ptr R, ulong * d, elem_srcptr A, long lenA, elem_srcptr B, long lenB, const ring_t ring);
void elem_poly_pseudo_divrem(elem_poly_t Q, elem_poly_t R, ulong * d, const elem_poly_t A, const elem_poly_t B, const ring_t ring);

void _elem_poly_pow_ui(elem_ptr res, elem_srcptr poly, long len, ulong e, const ring_t ring);
void elem_poly_pow_ui(elem_poly_t res, const elem_poly_t poly, ulong exp, const ring_t ring);

void _elem_poly_gcd_subresultant(elem_ptr res, elem_srcptr poly1, long len1, elem_srcptr poly2, long len2, const ring_t ring);
void elem_poly_gcd_subresultant(elem_poly_t res, const elem_poly_t poly1, const elem_poly_t poly2, const ring_t ring);

static __inline__ void
elem_poly_set(elem_poly_t res, const elem_poly_t src, const ring_t ring)
{
    long len = src->length;
    elem_poly_fit_length(res, len, ring);
    _elem_vec_set(res->coeffs, src->coeffs, len, ring->parent);
    elem_poly_set_length(res, len, ring);
}

static __inline__ void
elem_poly_neg(elem_poly_t res, const elem_poly_t src, const ring_t ring)
{
    long len = src->length;
    elem_poly_fit_length(res, len, ring);
    _elem_vec_neg(res->coeffs, src->coeffs, len, ring->parent);
    elem_poly_set_length(res, len, ring);
}

/* deprecate? */
void gen_set_coeff_si(gen_t x, long index, long value);

/* Matrix functions */

#define MAT_INDEX(rows, i, j, ring) (INDEX((rows)[i], (j), (ring)->size))
#define MAT_SRCINDEX(rows, i, j, ring) (SRC_INDEX((rows)[i], (j), (ring)->size))

void elem_mat_init(elem_mat_t mat, long rows, long cols, const ring_t ring);
void elem_mat_clear(elem_mat_t mat, const ring_t ring);
void elem_mat_set(elem_mat_t B, const elem_mat_t A, const ring_t ring);
void elem_mat_swap(elem_mat_t mat1, elem_mat_t mat2, const ring_t ring);
void elem_mat_zero(elem_mat_t mat, const ring_t ring);
void elem_mat_one(elem_mat_t mat, const ring_t ring);
int elem_mat_is_zero(const elem_mat_t mat, const ring_t ring);
int elem_mat_equal(const elem_mat_t A, const elem_mat_t B, const ring_t ring);
void elem_mat_neg(elem_mat_t B, const elem_mat_t A, const ring_t ring);
void elem_mat_add(elem_mat_t C, const elem_mat_t A, const elem_mat_t B, const ring_t ring);
void elem_mat_sub(elem_mat_t C, const elem_mat_t A, const elem_mat_t B, const ring_t ring);
void elem_mat_mul(elem_mat_t C, const elem_mat_t A, const elem_mat_t B, const ring_t ring);
void elem_mat_scalar_mul(elem_mat_t res, const elem_mat_t mat, elem_srcptr c, const ring_t ring);
void elem_mat_randtest(elem_mat_t mat, flint_rand_t state, const long * size, const ring_t ring);
void elem_mat_transpose(elem_mat_t B, const elem_mat_t A, const ring_t ring);
void elem_mat_print(const elem_mat_t mat, const ring_t ring);
long elem_mat_fflu(elem_mat_t B, elem_ptr den, long * perm, const elem_mat_t A, int rank_check, const ring_t ring);
int elem_mat_solve(elem_mat_t X, elem_ptr den, const elem_mat_t A, const elem_mat_t B, const ring_t ring);
void elem_mat_solve_fflu_precomp(elem_mat_t X, const long * perm, const elem_mat_t FFLU, const elem_mat_t B, const ring_t ring);
long elem_mat_rank(const elem_mat_t mat, const ring_t ring);
void elem_mat_det(elem_ptr det, const elem_mat_t mat, const ring_t ring);
long elem_mat_rref(elem_mat_t R, elem_ptr den, const elem_mat_t A, const ring_t ring);
long elem_mat_nullspace(elem_mat_t res, const elem_mat_t mat, const ring_t ring);
int elem_mat_inv(elem_mat_t X, elem_ptr den, const elem_mat_t A, const ring_t ring);

static __inline__ int elem_mat_is_empty(const elem_mat_t mat, const ring_t ring)
{
    return mat->r == 0 || mat->c == 0;
}

static __inline__ long elem_mat_nrows(const elem_mat_t mat, const ring_t ring)
{
    return mat->r;
}

static __inline__ long elem_mat_ncols(const elem_mat_t mat, const ring_t ring)
{
    return mat->c;
}

static __inline__ elem_ptr elem_mat_entry(const elem_mat_t mat, long i, long j, const ring_t ring)
{
    return MAT_INDEX(mat->rows, i, j, RING_PARENT(ring));
}

#ifdef __cplusplus
}
#endif

#endif

