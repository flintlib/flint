/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_ORE_POLY_H
#define GR_ORE_POLY_H

#ifdef GR_ORE_POLY_INLINES_C
#define GR_ORE_POLY_INLINE
#else
#define GR_ORE_POLY_INLINE static inline
#endif

#include "gr.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    ORE_ALGEBRA_DERIVATIVE,
    ORE_ALGEBRA_EULER_DERIVATIVE,
    ORE_ALGEBRA_FORWARD_SHIFT,
    ORE_ALGEBRA_FORWARD_DIFFERENCE,
}
ore_algebra_t;

#define ORE_POLY_NUM_ALGEBRAS 2

GR_ORE_POLY_INLINE
ore_algebra_t ore_algebra_randtest(flint_rand_t state)
{
   return (ore_algebra_t) n_randint(state, ORE_POLY_NUM_ALGEBRAS);
}

/* Compatible with gr_poly_struct */
typedef struct
{
    gr_ptr coeffs;
    slong alloc;
    slong length;
}
gr_ore_poly_struct;

typedef gr_ore_poly_struct gr_ore_poly_t[1];

/* Compatible with polynomial_ctx_t */
typedef struct
{
    gr_ctx_struct * base_ring;
    slong degree_limit;
    char * var;
    ore_algebra_t which_algebra;
    slong base_var;
    /* TODO: Function pointers to \sigma, \delta. */
}
_gr_ore_poly_ctx_struct;

typedef gr_ctx_struct gr_ore_poly_ctx_struct;

typedef gr_ore_poly_ctx_struct gr_ore_poly_ctx_t[1];

#define GR_ORE_POLY_CTX(ring_ctx) ((_gr_ore_poly_ctx_struct *)((ring_ctx)))
#define GR_ORE_POLY_ELEM_CTX(ring_ctx) (GR_ORE_POLY_CTX(ring_ctx)->base_ring)

/* Context object */

void gr_ctx_init_gr_ore_poly(gr_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_algebra_t which_algebra);

void gr_ore_poly_ctx_init(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_algebra_t which_algebra);
void gr_ore_poly_ctx_clear(gr_ore_poly_ctx_t ctx);

void gr_ore_poly_ctx_init_rand(gr_ore_poly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring);

WARN_UNUSED_RESULT int _gr_ore_poly_ctx_set_gen_name(gr_ctx_t ctx, const char * s);
WARN_UNUSED_RESULT int _gr_ore_poly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s);
WARN_UNUSED_RESULT int gr_ore_poly_gens_recursive(gr_vec_t vec, gr_ore_poly_ctx_t ctx);

int gr_ore_poly_ctx_write(gr_stream_t out, gr_ore_poly_ctx_t ctx);

truth_t gr_ore_poly_ctx_is_ring(gr_ore_poly_ctx_t ctx);
truth_t gr_ore_poly_ctx_is_zero_ring(gr_ore_poly_ctx_t ctx);
truth_t gr_ore_poly_ctx_is_commutative_ring(gr_ore_poly_ctx_t ctx);
truth_t gr_ore_poly_ctx_is_integral_domain(gr_ore_poly_ctx_t ctx);
truth_t gr_ore_poly_ctx_is_threadsafe(gr_ore_poly_ctx_t ctx);

/* Memory management */

void gr_ore_poly_init(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);
void gr_ore_poly_init2(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx);
void gr_ore_poly_clear(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);

GR_ORE_POLY_INLINE gr_ptr
gr_ore_poly_coeff_ptr(gr_ore_poly_t poly, slong i, gr_ore_poly_ctx_t ctx)
{
    return GR_ENTRY(poly->coeffs, i, GR_ORE_POLY_ELEM_CTX(ctx)->sizeof_elem);
}

GR_ORE_POLY_INLINE gr_srcptr
gr_ore_poly_coeff_srcptr(const gr_ore_poly_t poly, slong i, gr_ore_poly_ctx_t ctx)
{
    return GR_ENTRY(poly->coeffs, i, GR_ORE_POLY_ELEM_CTX(ctx)->sizeof_elem);
}

GR_ORE_POLY_INLINE slong gr_ore_poly_length(const gr_ore_poly_t poly, gr_ore_poly_ctx_t FLINT_UNUSED(ctx))
{
    return poly->length;
}

GR_ORE_POLY_INLINE void
gr_ore_poly_swap(gr_ore_poly_t poly1, gr_ore_poly_t poly2, gr_ore_poly_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(gr_ore_poly_struct, *poly1, *poly2);
}

GR_ORE_POLY_INLINE void
gr_ore_poly_set_shallow(gr_ore_poly_t res, const gr_ore_poly_t x, const gr_ore_poly_ctx_t ctx)
{
    *res = *x;
}

void gr_ore_poly_fit_length(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx);
void _gr_ore_poly_set_length(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx);
void _gr_ore_poly_normalise(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);

/* Basic manipulation */

WARN_UNUSED_RESULT int gr_ore_poly_set(gr_ore_poly_t res, const gr_ore_poly_t src, gr_ore_poly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ore_poly_truncate(gr_ore_poly_t poly, const gr_ore_poly_t src, slong newlen, gr_ore_poly_ctx_t ctx);

GR_ORE_POLY_INLINE WARN_UNUSED_RESULT int
gr_ore_poly_zero(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    _gr_ore_poly_set_length(poly, 0, ctx);
    return GR_SUCCESS;
}

WARN_UNUSED_RESULT int gr_ore_poly_one(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_neg_one(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_gen(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);

truth_t _gr_ore_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx);
truth_t gr_ore_poly_equal(const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx);

truth_t gr_ore_poly_is_zero(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);
truth_t gr_ore_poly_is_one(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);
truth_t gr_ore_poly_is_gen(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);

/* Input and output */

int _gr_ore_poly_write(gr_stream_t out, gr_srcptr poly, slong len, gr_ore_poly_ctx_t ctx);
int gr_ore_poly_write(gr_stream_t out, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);
int _gr_ore_poly_get_str(char ** res, gr_srcptr f, slong len, gr_ore_poly_ctx_t ctx);
int gr_ore_poly_get_str(char ** res, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);
int _gr_ore_poly_set_str(gr_ptr res, const char * s, slong len, gr_ore_poly_ctx_t ctx);
int gr_ore_poly_set_str(gr_ore_poly_t res, const char * s, gr_ore_poly_ctx_t ctx);
int gr_ore_poly_print(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);

/* Random generation */

int gr_ore_poly_randtest(gr_ore_poly_t poly, flint_rand_t state, slong len, gr_ore_poly_ctx_t ctx);

GR_ORE_POLY_INLINE WARN_UNUSED_RESULT int
_gr_ore_poly_randtest_default(gr_ore_poly_t res, flint_rand_t state, gr_ore_poly_ctx_t ctx)
{
    return gr_ore_poly_randtest(res, state, n_randint(state, 5), ctx);
}

/* Constants */

WARN_UNUSED_RESULT int gr_ore_poly_set_si(gr_ore_poly_t poly, slong x, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_set_ui(gr_ore_poly_t poly, ulong x, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_set_fmpz(gr_ore_poly_t poly, const fmpz_t x, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_set_fmpq(gr_ore_poly_t poly, const fmpq_t x, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_set_other(gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx);

/* Arithmetic */

WARN_UNUSED_RESULT int gr_ore_poly_neg(gr_ore_poly_t res, const gr_ore_poly_t src, gr_ore_poly_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ore_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_add(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ore_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_sub(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ore_poly_add_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_add_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_add_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_add_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_add_other(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ore_poly_sub_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_sub_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_sub_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_sub_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_sub_other(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ore_poly_mul_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_mul_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_mul_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_mul_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_other_mul(gr_ore_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
