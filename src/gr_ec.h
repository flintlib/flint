/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_EC_H
#define GR_EC_H

#ifdef GR_EC_INLINES_C
#define GR_EC_INLINE
#else
#define GR_EC_INLINE static inline
#endif

#include "gr.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Curve models */

typedef enum
{
    GR_EC_LONG_WEIERSTRASS,
    GR_EC_SHORT_WEIERSTRASS,

    GR_EC_NUM_MODELS
}
gr_ec_model_t;

/* Context object (the curve) */

typedef struct
{
    gr_ctx_struct * base_ring;
    /* a1, a2, a3, a4, a6, b2, b4, b6, b8, disc */
    gr_ptr coeffs;
    gr_ec_model_t model;
}
gr_ec_ctx_struct;

typedef gr_ec_ctx_struct gr_ec_ctx_t[1];

#define GR_EC_CTX_NUM_COEFFS 10

#define GR_EC_ELEM_CTX(ctx) ((ctx)->base_ring)
#define GR_EC_SIZEOF_ELEM(ctx) (GR_EC_ELEM_CTX(ctx)->sizeof_elem)

#define GR_EC_COEFF(ctx, i) GR_ENTRY((ctx)->coeffs, i, GR_EC_SIZEOF_ELEM(ctx))

#define GR_EC_A1(ctx) GR_EC_COEFF(ctx, 0)
#define GR_EC_A2(ctx) GR_EC_COEFF(ctx, 1)
#define GR_EC_A3(ctx) GR_EC_COEFF(ctx, 2)
#define GR_EC_A4(ctx) GR_EC_COEFF(ctx, 3)
#define GR_EC_A6(ctx) GR_EC_COEFF(ctx, 4)
#define GR_EC_B2(ctx) GR_EC_COEFF(ctx, 5)
#define GR_EC_B4(ctx) GR_EC_COEFF(ctx, 6)
#define GR_EC_B6(ctx) GR_EC_COEFF(ctx, 7)
#define GR_EC_B8(ctx) GR_EC_COEFF(ctx, 8)
#define GR_EC_DISC(ctx) GR_EC_COEFF(ctx, 9)

/* Points, homogeneous projective coordinates (X : Y : Z) */

typedef struct
{
    gr_ptr coords;      /* X, Y, Z */
}
gr_ec_point_struct;

typedef gr_ec_point_struct gr_ec_point_t[1];

/* Points, affine coordinates (x, y) */

typedef struct
{
    gr_ptr coords;      /* x, y */
    truth_t is_infinity;
}
gr_ec_aff_point_struct;

typedef gr_ec_aff_point_struct gr_ec_aff_point_t[1];

/* Points, Jacobian coordinates (X, Y, Z) with x = X/Z^2, y = Y/Z^3 */

typedef struct
{
    gr_ptr coords;      /* X, Y, Z */
    truth_t is_infinity;
}
gr_ec_jac_point_struct;

typedef gr_ec_jac_point_struct gr_ec_jac_point_t[1];

#define GR_EC_POINT_X(P, ctx) ((P)->coords)
#define GR_EC_POINT_Y(P, ctx) GR_ENTRY((P)->coords, 1, GR_EC_SIZEOF_ELEM(ctx))
#define GR_EC_POINT_Z(P, ctx) GR_ENTRY((P)->coords, 2, GR_EC_SIZEOF_ELEM(ctx))

#define GR_EC_AFF_POINT_X(P, ctx) ((P)->coords)
#define GR_EC_AFF_POINT_Y(P, ctx) GR_ENTRY((P)->coords, 1, GR_EC_SIZEOF_ELEM(ctx))

#define GR_EC_JAC_POINT_X(P, ctx) ((P)->coords)
#define GR_EC_JAC_POINT_Y(P, ctx) GR_ENTRY((P)->coords, 1, GR_EC_SIZEOF_ELEM(ctx))
#define GR_EC_JAC_POINT_Z(P, ctx) GR_ENTRY((P)->coords, 2, GR_EC_SIZEOF_ELEM(ctx))

/* Context object */

WARN_UNUSED_RESULT int gr_ec_ctx_init(gr_ec_ctx_t ctx, gr_ctx_t base_ring, gr_srcptr a1, gr_srcptr a2, gr_srcptr a3, gr_srcptr a4, gr_srcptr a6);
WARN_UNUSED_RESULT int gr_ec_ctx_init_si(gr_ec_ctx_t ctx, gr_ctx_t base_ring, slong a1, slong a2, slong a3, slong a4, slong a6);
WARN_UNUSED_RESULT int gr_ec_ctx_init_short_weierstrass(gr_ec_ctx_t ctx, gr_ctx_t base_ring, gr_srcptr a4, gr_srcptr a6);
WARN_UNUSED_RESULT int gr_ec_ctx_init_short_weierstrass_si(gr_ec_ctx_t ctx, gr_ctx_t base_ring, slong a4, slong a6);
WARN_UNUSED_RESULT int gr_ec_ctx_init_randtest(gr_ec_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring);

void gr_ec_ctx_clear(gr_ec_ctx_t ctx);

GR_EC_INLINE gr_ctx_struct * gr_ec_ctx_base_ring(gr_ec_ctx_t ctx)
{
    return GR_EC_ELEM_CTX(ctx);
}

GR_EC_INLINE gr_ec_model_t gr_ec_ctx_model(gr_ec_ctx_t ctx)
{
    return ctx->model;
}

GR_EC_INLINE truth_t gr_ec_ctx_is_over_field(gr_ec_ctx_t ctx)
{
    return gr_ctx_is_field(GR_EC_ELEM_CTX(ctx));
}

int gr_ec_ctx_write(gr_stream_t out, gr_ec_ctx_t ctx);
int gr_ec_ctx_get_str(char ** res, gr_ec_ctx_t ctx);
int gr_ec_ctx_print(gr_ec_ctx_t ctx);

/* Curve invariants */

GR_EC_INLINE gr_srcptr gr_ec_ctx_a_invariants_srcptr(gr_ec_ctx_t ctx)
{
    return GR_EC_A1(ctx);
}

GR_EC_INLINE gr_srcptr gr_ec_ctx_b_invariants_srcptr(gr_ec_ctx_t ctx)
{
    return GR_EC_B2(ctx);
}

WARN_UNUSED_RESULT int gr_ec_ctx_a_invariants(gr_ptr a1, gr_ptr a2, gr_ptr a3, gr_ptr a4, gr_ptr a6, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_ctx_b_invariants(gr_ptr b2, gr_ptr b4, gr_ptr b6, gr_ptr b8, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_ctx_c_invariants(gr_ptr c4, gr_ptr c6, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_ctx_discriminant(gr_ptr res, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_ctx_j_invariant(gr_ptr res, gr_ec_ctx_t ctx);

truth_t gr_ec_ctx_is_smooth(gr_ec_ctx_t ctx);

/*
    Projective points
*/

/* Memory management */

void gr_ec_point_init(gr_ec_point_t P, gr_ec_ctx_t ctx);
void gr_ec_point_clear(gr_ec_point_t P, gr_ec_ctx_t ctx);

GR_EC_INLINE void
gr_ec_point_swap(gr_ec_point_t P, gr_ec_point_t Q, gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(gr_ec_point_struct, *P, *Q);
}

GR_EC_INLINE gr_ptr gr_ec_point_x_ptr(gr_ec_point_t P, gr_ec_ctx_t FLINT_UNUSED(ctx)) { return (P)->coords; }
GR_EC_INLINE gr_ptr gr_ec_point_y_ptr(gr_ec_point_t P, gr_ec_ctx_t ctx) { return GR_EC_POINT_Y(P, ctx); }
GR_EC_INLINE gr_ptr gr_ec_point_z_ptr(gr_ec_point_t P, gr_ec_ctx_t ctx) { return GR_EC_POINT_Z(P, ctx); }

GR_EC_INLINE gr_srcptr gr_ec_point_x_srcptr(const gr_ec_point_t P, gr_ec_ctx_t FLINT_UNUSED(ctx)) { return (P)->coords; }
GR_EC_INLINE gr_srcptr gr_ec_point_y_srcptr(const gr_ec_point_t P, gr_ec_ctx_t ctx) { return GR_EC_POINT_Y(P, ctx); }
GR_EC_INLINE gr_srcptr gr_ec_point_z_srcptr(const gr_ec_point_t P, gr_ec_ctx_t ctx) { return GR_EC_POINT_Z(P, ctx); }

/* Basic manipulation */

WARN_UNUSED_RESULT int gr_ec_point_set(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_zero(gr_ec_point_t res, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ec_point_set_affine(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_set_affine(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_point_set_projective(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_set_projective(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_normalize(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_point_lift_x(gr_ec_point_t res, gr_srcptr x, gr_ec_ctx_t ctx);

/* Comparisons and properties */

truth_t gr_ec_point_is_inf(const gr_ec_point_t P, gr_ec_ctx_t ctx);
truth_t gr_ec_point_equal(const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx);
truth_t gr_ec_point_is_on_curve(const gr_ec_point_t P, gr_ec_ctx_t ctx);

/* Input and output */

int gr_ec_point_write(gr_stream_t out, const gr_ec_point_t P, gr_ec_ctx_t ctx);
int gr_ec_point_get_str(char ** res, const gr_ec_point_t P, gr_ec_ctx_t ctx);
int gr_ec_point_print(const gr_ec_point_t P, gr_ec_ctx_t ctx);

/* Random generation */

WARN_UNUSED_RESULT int gr_ec_point_randtest(gr_ec_point_t res, flint_rand_t state, gr_ec_ctx_t ctx);

/* Arithmetic */

WARN_UNUSED_RESULT int gr_ec_point_neg(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_add(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_sub(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_dbl(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_point_mul_ui(gr_ec_point_t res, const gr_ec_point_t P, ulong n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_mul_si(gr_ec_point_t res, const gr_ec_point_t P, slong n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_mul_fmpz(gr_ec_point_t res, const gr_ec_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);

/*
    Affine points

    All operations require the base ring to be a field.
*/

/* Memory management */

void gr_ec_aff_point_init(gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
void gr_ec_aff_point_clear(gr_ec_aff_point_t P, gr_ec_ctx_t ctx);

GR_EC_INLINE void
gr_ec_aff_point_swap(gr_ec_aff_point_t P, gr_ec_aff_point_t Q, gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(gr_ec_aff_point_struct, *P, *Q);
}

GR_EC_INLINE gr_ptr gr_ec_aff_point_x_ptr(gr_ec_aff_point_t P, gr_ec_ctx_t FLINT_UNUSED(ctx)) { return (P)->coords; }
GR_EC_INLINE gr_ptr gr_ec_aff_point_y_ptr(gr_ec_aff_point_t P, gr_ec_ctx_t ctx) { return GR_EC_AFF_POINT_Y(P, ctx); }

GR_EC_INLINE gr_srcptr gr_ec_aff_point_x_srcptr(const gr_ec_aff_point_t P, gr_ec_ctx_t FLINT_UNUSED(ctx)) { return (P)->coords; }
GR_EC_INLINE gr_srcptr gr_ec_aff_point_y_srcptr(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx) { return GR_EC_AFF_POINT_Y(P, ctx); }

/* Basic manipulation */

WARN_UNUSED_RESULT int gr_ec_aff_point_set(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_aff_point_zero(gr_ec_aff_point_t res, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ec_aff_point_set_affine(gr_ec_aff_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_aff_point_set_affine(gr_ec_aff_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_aff_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_aff_point_lift_x(gr_ec_aff_point_t res, gr_srcptr x, gr_ec_ctx_t ctx);

/* Comparisons and properties */

GR_EC_INLINE truth_t
gr_ec_aff_point_is_inf(const gr_ec_aff_point_t P, gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    return P->is_infinity;
}

truth_t gr_ec_aff_point_equal(const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
truth_t gr_ec_aff_point_is_on_curve(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);

/* Input and output */

int gr_ec_aff_point_write(gr_stream_t out, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
int gr_ec_aff_point_get_str(char ** res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
int gr_ec_aff_point_print(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);

/* Random generation */

WARN_UNUSED_RESULT int gr_ec_aff_point_randtest(gr_ec_aff_point_t res, flint_rand_t state, gr_ec_ctx_t ctx);

/* Arithmetic */

WARN_UNUSED_RESULT int gr_ec_aff_point_neg(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_aff_point_add(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_aff_point_sub(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_aff_point_dbl(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_aff_point_mul_ui(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, ulong n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_aff_point_mul_si(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, slong n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_aff_point_mul_fmpz(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);

/*
    Jacobian points
*/

/* Memory management */

void gr_ec_jac_point_init(gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
void gr_ec_jac_point_clear(gr_ec_jac_point_t P, gr_ec_ctx_t ctx);

GR_EC_INLINE void
gr_ec_jac_point_swap(gr_ec_jac_point_t P, gr_ec_jac_point_t Q, gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(gr_ec_jac_point_struct, *P, *Q);
}

GR_EC_INLINE gr_ptr gr_ec_jac_point_x_ptr(gr_ec_jac_point_t P, gr_ec_ctx_t FLINT_UNUSED(ctx)) { return (P)->coords; }
GR_EC_INLINE gr_ptr gr_ec_jac_point_y_ptr(gr_ec_jac_point_t P, gr_ec_ctx_t ctx) { return GR_EC_JAC_POINT_Y(P, ctx); }
GR_EC_INLINE gr_ptr gr_ec_jac_point_z_ptr(gr_ec_jac_point_t P, gr_ec_ctx_t ctx) { return GR_EC_JAC_POINT_Z(P, ctx); }

GR_EC_INLINE gr_srcptr gr_ec_jac_point_x_srcptr(const gr_ec_jac_point_t P, gr_ec_ctx_t FLINT_UNUSED(ctx)) { return (P)->coords; }
GR_EC_INLINE gr_srcptr gr_ec_jac_point_y_srcptr(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx) { return GR_EC_JAC_POINT_Y(P, ctx); }
GR_EC_INLINE gr_srcptr gr_ec_jac_point_z_srcptr(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx) { return GR_EC_JAC_POINT_Z(P, ctx); }

/* Basic manipulation */

WARN_UNUSED_RESULT int gr_ec_jac_point_set(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_zero(gr_ec_jac_point_t res, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ec_jac_point_set_affine(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_set_affine(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_set_jacobian(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_set_jacobian(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_jac_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_normalize(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);

/* Comparisons and properties */

GR_EC_INLINE truth_t
gr_ec_jac_point_is_inf(const gr_ec_jac_point_t P, gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    return P->is_infinity;
}

truth_t gr_ec_jac_point_equal(const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx);
truth_t gr_ec_jac_point_is_on_curve(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);

/* Input and output */

int gr_ec_jac_point_write(gr_stream_t out, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
int gr_ec_jac_point_get_str(char ** res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
int gr_ec_jac_point_print(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);

/* Random generation */

WARN_UNUSED_RESULT int gr_ec_jac_point_randtest(gr_ec_jac_point_t res, flint_rand_t state, gr_ec_ctx_t ctx);

/* Arithmetic */

WARN_UNUSED_RESULT int gr_ec_jac_point_neg(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_add(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_sub(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_dbl(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_jac_point_add_aff_point(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_sub_aff_point(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_jac_point_mul_ui(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, ulong n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_mul_si(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, slong n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_mul_fmpz(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);

/* Conversions between representations */

WARN_UNUSED_RESULT int gr_ec_point_set_aff_point(gr_ec_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_point_set_jac_point(gr_ec_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_aff_point_set_point(gr_ec_aff_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_aff_point_set_jac_point(gr_ec_aff_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int gr_ec_jac_point_set_point(gr_ec_jac_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ec_jac_point_set_aff_point(gr_ec_jac_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);

/* Model-specific implementations */

WARN_UNUSED_RESULT int _gr_ec_point_add_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_point_add_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_point_dbl_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_point_dbl_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_point_mul_fmpz_binary(gr_ec_point_t res, const gr_ec_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ec_aff_point_add_long_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_aff_point_add_short_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_aff_point_dbl_long_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_aff_point_dbl_short_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_aff_point_mul_fmpz_binary(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ec_jac_point_add_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_mul_fmpz_binary(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_mul_fmpz_naf(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);

/* Simultaneous normalization: one inversion for the whole vector */
WARN_UNUSED_RESULT int _gr_ec_jac_point_vec_get_aff_point_vec(gr_ec_aff_point_struct * res, const gr_ec_jac_point_struct * P, slong len, gr_ec_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* GR_EC_H */
