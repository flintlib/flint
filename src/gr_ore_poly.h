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
int gr_ore_poly_ctx_init_q_shift(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, gr_srcptr q);
int gr_ore_poly_ctx_init_mahler(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, slong mahler_base);
void gr_ore_poly_ctx_init_custom(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, const gr_ore_poly_sigma_delta_t sigma_delta, void * ore_data);
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
gr_ore_poly_set_shallow(gr_ore_poly_t res, const gr_ore_poly_t x, const gr_ore_poly_ctx_t FLINT_UNUSED(ctx))
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
WARN_UNUSED_RESULT int sigma_delta_compose(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx);
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

WARN_UNUSED_RESULT int gr_ore_poly_apply_custom(gr_ptr res, const gr_ore_poly_t P, gr_srcptr f, gr_srcptr d1, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_apply(gr_ptr res, const gr_ore_poly_t P, gr_srcptr f, gr_ore_poly_ctx_t ctx);

/* Conversions between operator types.

   These act on raw coefficient arrays and take the coefficient ring context
   (not the Ore polynomial ring context). The variable x is the generator of the
   coefficient ring of index var, i.e. entry var of gr_gens(ctx). Since x need
   not be invertible, ddx_to_euler multiplies the operator by x^(len-1);
   euler_to_ddx leaves it unchanged. The result res has the same length len as
   op. */
WARN_UNUSED_RESULT int _gr_ore_poly_ddx_to_euler(gr_ptr res, gr_srcptr op, slong len, slong var, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ore_poly_euler_to_ddx(gr_ptr res, gr_srcptr op, slong len, slong var, gr_ctx_t ctx);

/* Convert a length-len operator between the shift/difference algebra types
   src_alg and dst_alg. Each generator is a degree-one Laurent polynomial in the
   forward shift S, so on output *p and res satisfy
   S^p * res_in_dst = op_in_src (a power of S appears only when crossing between
   the forward side S, S-1 and the backward side S^-1, 1-S^-1; it then needs a
   GR_CTX_GR_POLY base ring). Returns GR_DOMAIN if either algebra is not a
   shift/difference type. The variable is the base ring generator of index var;
   only var == 0 over a univariate polynomial ring is implemented (else
   GR_UNABLE). */
WARN_UNUSED_RESULT int _gr_ore_poly_shift_convert(gr_ptr res, slong * p, gr_srcptr op, slong len, ore_algebra_t src_alg, ore_algebra_t dst_alg, slong var, gr_ctx_t ctx);

/* Exact bridge of the isomorphism theta <-> n, x <-> backward shift S^{-1}, on
   raw coefficient vectors over a GR_CTX_GR_POLY base ring (else GR_UNABLE). The
   caller allocates res to reslen = 1 + the maximum coefficient degree of op. */
WARN_UNUSED_RESULT int _gr_ore_poly_euler_to_backshift_univar(gr_ptr res, slong reslen, gr_srcptr op, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ore_poly_backshift_to_euler_univar(gr_ptr res, slong reslen, gr_srcptr op, slong len, gr_ctx_t ctx);

/* Convert between differential operators (op_ctx/res_ctx in {DERIVATIVE,
   EULER_DERIVATIVE}) and shift/difference operators (in the four shift types),
   via the generating-series isomorphism theta <-> n, x <-> backward shift.
   Both contexts must share a GR_CTX_GR_POLY base ring. On output *p is such that
   C^p * res corresponds to op under the isomorphism, where C is the companion
   generator on the res side (S for *_to_shift, x for *_to_differential).
   Returns GR_DOMAIN if a context algebra is outside its family. */
WARN_UNUSED_RESULT int gr_ore_poly_differential_to_shift(gr_ore_poly_t res, slong * p, const gr_ore_poly_t op, gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx);
WARN_UNUSED_RESULT int gr_ore_poly_shift_to_differential(gr_ore_poly_t res, slong * p, const gr_ore_poly_t op, gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx);

/* Convert an operator between two algebras of the same family (both differential
   or both shift/difference), dispatching to the conversions above. On output *p
   satisfies C^p * res = op, where the companion generator C depends on the pair
   (x for differential conversions, the forward shift S for shift ones).
   Returns GR_UNABLE for an unsupported pair (crossing the differential/shift
   boundary, or involving an algebra outside those two families). */
WARN_UNUSED_RESULT int gr_ore_poly_convert(gr_ore_poly_t res, slong * p, const gr_ore_poly_t op, gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx);

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

WARN_UNUSED_RESULT int _gr_ore_poly_lmul_gen(gr_ptr res, gr_srcptr poly, slong len, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_lmul_gen(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ore_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_mul(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ore_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr U, slong lenU, gr_srcptr V, slong lenV, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_divrem(gr_ore_poly_t Q, gr_ore_poly_t R, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_div(gr_ore_poly_t Q, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_ore_poly_rem(gr_ore_poly_t R, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
