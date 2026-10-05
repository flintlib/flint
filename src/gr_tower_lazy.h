/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_TOWER_LAZY_H
#define GR_TOWER_LAZY_H

/*
    Lazy exact fields over towers (gr_tower.h): the complex field of the
    numbers the towers can represent, its real and algebraic subfields,
    and their views. Elements carry a pointer to a shared tower which is
    created and extended on demand; the zero test is complete (under
    Schanuel's conjecture for termination).
*/

#include "fmpq_poly.h"
#include "gr_tower.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Elements of the lazy fields. The fields may be read (they are
    documented in gr_tower.rst) but not written: an element refers to
    a tower shared with other elements, and its data is converted and
    reduced lazily under the context's lock. The accessors below bring
    an element up to date before exposing its data. The layout may
    change between versions of FLINT.
*/
#define GR_TOWER_LAZY_REPR_FLAT 0       /* elem.flat: fmpz_mpoly_q in the polynomial context elem.flat.mctx */
#define GR_TOWER_LAZY_REPR_RATIONAL 1   /* elem.q: a rational number (element of the trivial tower) */
#define GR_TOWER_LAZY_REPR_DENSE 2      /* elem.dense: a polynomial in the first generator of the tower */

struct _gr_tower_dense_field_struct;

typedef struct
{
    gr_tower_flat_struct * F;       /* the tower (via F->T) and its flat machinery */
    slong level;                    /* prefix length, valid for structure_version of the tower */
    ulong structure_version;
    ulong reduced_version;          /* ideal version of F for which data is known to be reduced (0: unknown) */
    slong shallow;                  /* the data may be shared with another element (internal) */
    slong repr;                     /* GR_TOWER_LAZY_REPR_* */
#ifdef LAZY_UNION_DEBUG
    struct
#else
    union
#endif
    {
        fmpq q;                                    /* RATIONAL */
        struct
        {
            fmpz_mpoly_q_struct data;
            fmpz_mpoly_ctx_struct * mctx;          /* context in which data is stored */
        }
        flat;                                      /* FLAT */
        struct
        {
            fmpq_poly_struct poly;                 /* canonical, of length < the degree of the field */
            const struct _gr_tower_dense_field_struct * nf;   /* the field: generator and modulus (internal) */
        }
        dense;                                     /* DENSE */
    }
    elem;
}
gr_tower_lazy_elem_struct;

typedef gr_tower_lazy_elem_struct gr_tower_lazy_elem_t[1];

/* element printing (the GR_TOWER_PRINT_* flags of gr_tower.h) */
void gr_tower_lazy_ctx_set_print(gr_ctx_t ctx, int flags, slong digits);

/* choice of generators in lazy fields (gr_tower_lazy_ctx_set_gen_flags):
   by default, the square root of a negative rational number -A B^2 (A
   squarefree) is B sqrt(-A) with sqrt(-A) a single generator (a root of
   x^2 + A; sqrt(-1) = i and sqrt(-3) = 2 zeta_3 + 1 are roots of unity);
   with GR_TOWER_GENS_SPLIT_IMAGINARY it is i B sqrt(A), with sqrt(A) the
   product of the roots of the primes dividing A. Roots of unity of
   composite orders are products of roots of prime power orders, and
   square roots of positive rationals products of square roots of primes,
   unless GR_TOWER_GENS_COMPOSITE_ROOTS (zeta_N as one generator, a root
   of Phi_N, of the least common multiple of the orders requested, as in
   Calcium) or GR_TOWER_GENS_COMPOSITE_RADICALS (sqrt(A) with A the
   squarefree part as one generator) is set */
#define GR_TOWER_GENS_SPLIT_IMAGINARY 1
#define GR_TOWER_GENS_COMPOSITE_ROOTS 2
#define GR_TOWER_GENS_COMPOSITE_RADICALS 4
#define GR_TOWER_COMPOSITE_ROOTS_DEGREE 128   /* the cyclotomic degree limit set by GR_TOWER_GENS_COMPOSITE_ROOTS when none is set */
void gr_tower_lazy_ctx_set_gen_flags(gr_ctx_t ctx, int flags);
/* tuning options (GR_TOWER_OPT_*) of a lazy field, shared by its towers
   and read whenever they apply (a change takes effect at once, also for
   existing elements); GR_DOMAIN for a value out of range */
int gr_tower_lazy_ctx_set_option(gr_ctx_t ctx, slong option, slong value);
slong gr_tower_lazy_ctx_get_option(gr_ctx_t ctx, slong option);
int gr_tower_lazy_ctx_gen_flags(gr_ctx_t ctx);
int gr_tower_lazy_ctx_print_flags(gr_ctx_t ctx);
int gr_tower_lazy_ctx_field_flags(gr_ctx_t ctx);

/* Lazy field of algebraic and transcendental numbers over the base field
   (the rational numbers): elements
   carry a pointer to a shared tower, created and extended on demand. */
void gr_ctx_init_tower_lazy(gr_ctx_t ctx, gr_ctx_t base, int flags);
/* A view of the lazy context parent sharing its towers and elements,
   restricted to a subfield (GR_TOWER_LAZY_REAL, GR_TOWER_LAZY_ALGEBRAIC
   added to the restrictions of parent); parent may be cleared first. */
void gr_ctx_init_tower_lazy_view(gr_ctx_t ctx, gr_ctx_t parent, int field_flags);
int gr_tower_lazy_ctx_same_state(gr_ctx_t ctx1, gr_ctx_t ctx2);

/* The methods of the gr interface of the lazy fields, callable directly
   (the element functions of gr.h dispatch to them). Rational elements,
   and dense elements of the same field, take lock-free paths; other
   operations hold the context's lock. */
void gr_tower_lazy_init(gr_tower_lazy_elem_t x, gr_ctx_t ctx);
void gr_tower_lazy_clear(gr_tower_lazy_elem_t x, gr_ctx_t ctx);
void gr_tower_lazy_swap(gr_tower_lazy_elem_t x, gr_tower_lazy_elem_t y, gr_ctx_t ctx);
void gr_tower_lazy_set_shallow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_randtest(gr_tower_lazy_elem_t res, flint_rand_t state, gr_ctx_t ctx);
int gr_tower_lazy_write(gr_stream_t out, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_set_str(gr_tower_lazy_elem_t res, const char * s, gr_ctx_t ctx);
#ifdef FEXPR_H
int gr_tower_lazy_get_fexpr(fexpr_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
#endif
int gr_tower_lazy_zero(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int gr_tower_lazy_one(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int gr_tower_lazy_pi(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int gr_tower_lazy_i(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int gr_tower_lazy_gens(gr_vec_t vec, gr_ctx_t ctx);
int gr_tower_lazy_set(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_set_si(gr_tower_lazy_elem_t res, slong c, gr_ctx_t ctx);
int gr_tower_lazy_set_ui(gr_tower_lazy_elem_t res, ulong c, gr_ctx_t ctx);
int gr_tower_lazy_set_fmpz(gr_tower_lazy_elem_t res, const fmpz_t c, gr_ctx_t ctx);
int gr_tower_lazy_set_fmpq(gr_tower_lazy_elem_t res, const fmpq_t c, gr_ctx_t ctx);
int gr_tower_lazy_set_other(gr_tower_lazy_elem_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
int gr_tower_lazy_get_si(slong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_get_ui(ulong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_get_fmpz(fmpz_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_get_fmpq(fmpq_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
truth_t gr_tower_lazy_is_zero(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
truth_t gr_tower_lazy_is_one(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
truth_t gr_tower_lazy_equal(const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_neg(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_add(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_sub(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_mul(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_div(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_inv(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_add_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx);
int gr_tower_lazy_add_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx);
int gr_tower_lazy_add_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx);
int gr_tower_lazy_add_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx);
int gr_tower_lazy_sub_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx);
int gr_tower_lazy_sub_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx);
int gr_tower_lazy_sub_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx);
int gr_tower_lazy_sub_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx);
int gr_tower_lazy_mul_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx);
int gr_tower_lazy_mul_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx);
int gr_tower_lazy_mul_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx);
int gr_tower_lazy_mul_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx);
int gr_tower_lazy_div_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx);
int gr_tower_lazy_div_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx);
int gr_tower_lazy_div_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx);
int gr_tower_lazy_div_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx);
int gr_tower_lazy_pow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_sqrt(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_abs(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_conj(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_exp(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_log(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_poly_mullow(gr_tower_lazy_elem_struct * res, const gr_tower_lazy_elem_struct * p1, slong len1, const gr_tower_lazy_elem_struct * p2, slong len2, slong n, gr_ctx_t ctx);
int gr_tower_lazy_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx);
int gr_tower_lazy_poly_factor(gr_poly_t c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx);
int gr_tower_lazy_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int gr_tower_lazy_mat_det(gr_tower_lazy_elem_t res, const gr_mat_t A, gr_ctx_t ctx);
int gr_tower_lazy_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

int gr_tower_lazy_get_qqbar(qqbar_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);

/* Elementary functions of the lazy field, through exp, log and sqrt;
   inverse functions are checked against the principal branch
   numerically (GR_UNABLE on the branch cuts). Rounding and comparisons
   apply to real elements (GR_DOMAIN otherwise). */
int gr_tower_lazy_sin(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_cos(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_tan(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_sinh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_cosh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_tanh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_asin(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_acos(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_atan(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_asinh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_acosh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_atanh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_re(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_im(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_arg(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_sgn(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_csgn(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
truth_t gr_tower_lazy_is_real(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_cmp(int * res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_cmpabs(int * res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_floor(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_ceil(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_trunc(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_nint(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_get_d(double * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_root_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong n, gr_ctx_t ctx);
int gr_tower_lazy_get_acb(acb_t res, const gr_tower_lazy_elem_t x, slong prec, gr_ctx_t ctx);

/* Special functions of the lazy field (see lazy_special.c): values at
   canonical arguments become generators, the functional equations and
   special values being applied first. */
int gr_tower_lazy_gamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_rgamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_beta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int gr_tower_lazy_digamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_polygamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_erf(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_erfc(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_erfi(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_lambertw(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_lambertw_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t k, gr_ctx_t ctx);
int gr_tower_lazy_zeta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_hurwitz_zeta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t a, gr_ctx_t ctx);
int gr_tower_lazy_polylog(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_dilog(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_elliptic_k(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_elliptic_e(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_euler(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int gr_tower_lazy_catalan(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
void gr_tower_lazy_ctx_stats(gr_ctx_t ctx);
slong gr_tower_lazy_ctx_num_towers(gr_ctx_t ctx);
gr_tower_struct * gr_tower_lazy_get_tower(slong * level, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
const fmpz_mpoly_q_struct * gr_tower_lazy_get_data(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_repr(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int gr_tower_lazy_get_fmpq_poly(fmpq_poly_t res, fmpz_poly_t modulus, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
