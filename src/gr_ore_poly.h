/*
    Copyright (C) 2025 Ricardo Buring
    Copyright (C) 2025 Marc Mezzarobba

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

/* Ore polynomial rings */

typedef enum {
    ORE_ALGEBRA_CUSTOM,
    ORE_ALGEBRA_COMMUTATIVE,
    ORE_ALGEBRA_DERIVATIVE,
    ORE_ALGEBRA_EULER_DERIVATIVE,
    ORE_ALGEBRA_FORWARD_SHIFT,
    ORE_ALGEBRA_FORWARD_DIFFERENCE,
    ORE_ALGEBRA_BACKWARD_SHIFT,
    ORE_ALGEBRA_BACKWARD_DIFFERENCE,
    ORE_ALGEBRA_Q_SHIFT,
    ORE_ALGEBRA_MAHLER,
    /* todo: q-derivative? general substitution?... */
    ORE_ALGEBRA_FROBENIUS,

    ORE_POLY_NUM_ALGEBRAS
}
ore_algebra_t;

ore_algebra_t ore_algebra_randtest(flint_rand_t state);

/* Compatible with gr_poly_struct */
typedef struct
{
    gr_ptr coeffs;
    slong alloc;
    slong length;
}
gr_ore_poly_struct;

typedef gr_ore_poly_struct gr_ore_poly_t[1];

typedef gr_ctx_struct gr_ore_poly_ctx_struct;

typedef gr_ore_poly_ctx_struct gr_ore_poly_ctx_t[1];

typedef int (* gr_ore_poly_sigma_delta_t) (gr_ptr, gr_ptr, gr_srcptr, gr_ore_poly_ctx_struct *);

typedef struct
{
    slong base_var;
    gr_ptr sigma_x;
    union {
        gr_ptr q;
        slong mahler_base;
    };
}
gr_ore_poly_ore_data_t;

/* Compatible with polynomial_ctx_t */
typedef struct
{
    gr_ctx_struct * base_ring;
    slong degree_limit;
    char * var;
    ore_algebra_t which_algebra;
    gr_ore_poly_sigma_delta_t sigma_delta;
    void * ore_data;
}
_gr_ore_poly_ctx_struct;

#define GR_ORE_POLY_CTX(ring_ctx) ((_gr_ore_poly_ctx_struct *)((ring_ctx)))
#define GR_ORE_POLY_ELEM_CTX(ring_ctx) (GR_ORE_POLY_CTX(ring_ctx)->base_ring)
#define GR_ORE_POLY_ORE_DATA(ring_ctx) ((gr_ore_poly_ore_data_t *)(GR_ORE_POLY_CTX(ring_ctx)->ore_data))

GR_ORE_POLY_INLINE void *
gr_ore_poly_ctx_data_ptr(gr_ore_poly_ctx_t ctx)
{
    return GR_ORE_POLY_CTX(ctx)->ore_data;
}

/* Context object */

void gr_ctx_init_gr_ore_poly(gr_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_algebra_t which_algebra);

void gr_ore_poly_ctx_init(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_algebra_t which_algebra);
void gr_ore_poly_ctx_clear(gr_ore_poly_ctx_t ctx);

void gr_ore_poly_ctx_init_randtest(gr_ore_poly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring);
void gr_ore_poly_ctx_init_randtest2(gr_ctx_t base_ring, gr_ore_poly_ctx_t ctx, flint_rand_t state);

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

/* Action on the base ring */

WARN_UNUSED_RESULT int sigma_delta_unable(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_commutative(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_derivative(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_euler_derivative(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_forward_shift(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_backward_shift(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_forward_difference(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_backward_difference(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_q_shift(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_mahler(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
WARN_UNUSED_RESULT int sigma_delta_frobenius(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);

GR_ORE_POLY_INLINE WARN_UNUSED_RESULT int
gr_ore_poly_sigma_delta(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_t ctx)
{
    return GR_ORE_POLY_CTX(ctx)->sigma_delta(sigma, delta, a, ctx);
}

GR_ORE_POLY_INLINE WARN_UNUSED_RESULT int
gr_ore_poly_sigma(gr_ptr res, gr_srcptr a, gr_ore_poly_ctx_t ctx)
{
    return GR_ORE_POLY_CTX(ctx)->sigma_delta(res, NULL, a, ctx);
}

GR_ORE_POLY_INLINE WARN_UNUSED_RESULT int
gr_ore_poly_delta(gr_ptr res, gr_srcptr a, gr_ore_poly_ctx_t ctx)
{
    return GR_ORE_POLY_CTX(ctx)->sigma_delta(NULL, res, a, ctx);
}

extern const gr_ore_poly_sigma_delta_t _gr_ore_poly_default_sigma_delta[];

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
