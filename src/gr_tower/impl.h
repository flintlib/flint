/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_TOWER_IMPL_H
#define GR_TOWER_IMPL_H

#include "acb_poly.h"
#include "fmpq_mat.h"
#include "fmpq_poly.h"
#include "fmpz_mat.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

#define GR_TOWER_DEFAULT_PREC 64

void _gr_tower_fit_length(gr_tower_t T, slong len);
void _gr_tower_gen_set_status(gr_tower_gen_struct * g, int status);
int _gr_tower_adjoin_special_multi_flat_nocheck(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_struct * u, slong nargs, const fmpz_mpoly_ctx_t mctx, const char * name);

/* Arguments of generators (gen_args.c): arg (index 0, if present) and
   the additional arguments xargs (indices 1, 2, ...) */
slong _gr_tower_gen_num_args(const gr_tower_gen_struct * g);
gr_tower_flat_elem_struct * _gr_tower_gen_arg_ptr(const gr_tower_gen_struct * g, slong i);
void _gr_tower_gen_xargs_clear(gr_tower_gen_struct * g);
char * _gr_tower_gen_args_str(const gr_tower_gen_struct * g, gr_tower_t T);
int _gr_tower_gen_args_get_acb(acb_ptr res, const gr_tower_gen_struct * g, slong prec, gr_tower_flat_struct * F);
int _gr_tower_gen_copy_def_map(gr_tower_gen_struct * ng, const gr_tower_gen_struct * g, gr_tower_map_t map);

int _gr_tower_adjoin_special_flat_nocheck(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name);
gr_tower_gen_struct * _gr_tower_push_step(gr_tower_t T, const gr_poly_t m, const acb_t z, slong prec, int status, const char * name);

int _gr_tower_lazy_realify_locked(gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int _gr_tower_lazy_trig_locked(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, int kind, gr_ctx_t ctx);
int _gr_tower_lazy_pi_part_locked(fmpq_t r, gr_tower_lazy_elem_t y, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int _gr_tower_lazy_rational_repr_locked(fmpq_t c, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int _gr_tower_lazy_root_of_unity_angle_locked(fmpq_t r, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int _gr_tower_lazy_cyclotomic_eval(gr_tower_lazy_elem_t res, const fmpq_poly_t E, ulong m, gr_ctx_t ctx);
int _gr_tower_lazy_real_cyclotomic_eval(gr_tower_lazy_elem_t res, const fmpq_poly_t E, ulong n, gr_ctx_t ctx);
int _gr_tower_lazy_trig_pi_real_locked(gr_tower_lazy_elem_t res, const fmpq_t r, int which, gr_ctx_t ctx);
int _gr_tower_lazy_special_gen_locked(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, int kind, slong param, gr_ctx_t ctx);
int _gr_tower_lazy_special(gr_tower_lazy_elem_t res, int kind, slong param, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int _gr_tower_lazy_dirichlet_l_prim(gr_ptr res, gr_srcptr s, ulong q, ulong k, gr_ctx_t ctx);
int _gr_tower_lazy_hurwitz_rational(gr_ptr res, gr_srcptr s, slong p, slong q, gr_ctx_t ctx);
int _gr_tower_lazy_hurwitz_general(gr_ptr res, gr_srcptr s, gr_srcptr a, gr_ctx_t ctx);
/* functions of several arguments (lazy_multi.c, lazy_hypgeom.c): args[0] is
   the argument arg of the generator, as in gr_tower_adjoin_special_multi_flat */
int _gr_tower_lazy_special_multi(gr_tower_lazy_elem_t res, int kind, slong param, const gr_tower_lazy_elem_struct * args, slong nargs, gr_ctx_t ctx);
int _gr_tower_lazy_special_gen_multi_locked(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_struct * args, slong nargs, int kind, slong param, gr_ctx_t ctx);
int _gr_tower_lazy_special_alg_gen_multi_locked(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_struct * args, slong n, int kind, slong param, const gr_tower_lazy_elem_t W, ulong r, gr_ctx_t ctx);
int _gr_tower_lazy_elliptic_gen(gr_ptr res, gr_srcptr m, int kind, gr_ctx_t ctx);
/* the linking steps (Landen, quadratic transformations): nesting, and
   whether a generator K, E (or 2F1) is at the argument t numerically */
int _gr_tower_lazy_hyp_anchored(gr_ctx_t ctx, int delta);
int _gr_tower_lazy_hyp_anchor_present(const acb_t t, int elliptic_only, gr_ctx_t ctx);
int _gr_tower_lazy_elliptic_hypgeom(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t m, int kind, gr_ctx_t ctx);
int _gr_tower_lazy_modular_lambda(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx);
/* whether x is a rational number, exactly (a rational of small height
   close to x, verified): 1, 0 (not known to be), -1 (unknown) */
int _gr_tower_lazy_rational_recognize(fmpq_t c, gr_srcptr x, gr_ctx_t ctx);
int _gr_tower_lazy_elliptic_cm(gr_ptr res, int * done, gr_srcptr m, int kind, gr_ctx_t ctx);
int _gr_tower_lazy_jacobi_theta_j(gr_ptr res, slong j, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx);
int _gr_tower_lazy_hypgeom_args(gr_tower_lazy_elem_t res, slong param, const gr_tower_lazy_elem_struct * args, slong nargs, gr_ctx_t ctx);
truth_t _gr_tower_lazy_is_algebraic_repr_locked(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int _gr_tower_poly_norm(gr_poly_t res, const gr_poly_t M, slong j, gr_tower_t T);
int _gr_tower_base_poly_factor(gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t N, gr_tower_t T);
int _gr_tower_poly_factor_squarefree_trager(gr_vec_t fac, const gr_poly_t h, slong k, gr_tower_t T);
int _gr_tower_poly_get_acb_poly(acb_poly_t res, const gr_poly_t f, slong k, slong prec, gr_tower_t T);
int _gr_tower_newton_step(acb_t res, const acb_poly_t f, const acb_t z, slong wp);
int _gr_tower_certify_root(acb_t res, const gr_poly_t m, slong k, const acb_t z, slong prec, gr_tower_t T);
int _gr_tower_annihilating_poly(gr_poly_t res, gr_srcptr x, gr_tower_t T);
int _gr_tower_flat_grow(gr_tower_flat_t F);
void _gr_tower_move_to_front(gr_tower_t T, slong d);
void _gr_tower_move_gen(gr_tower_t T, slong d, slong p);
int _gr_tower_set_linear_gens(gr_tower_t T, slong n, const slong * d, const fmpz_mpoly_q_struct * v, const fmpz_mpoly_ctx_t mctx);
fmpz_mpoly_q_struct * _gr_tower_flat_stale_alloc(gr_tower_flat_t F);

/* Number fields of one generator in which elements of lazy fields have
   a dense form (lazy_dense.c): the generator (by gid) and its monic integral
   modulus at the time, immutable (dense elements refer to them without
   locking), kept in a list in F and freed with it. */
typedef struct _gr_tower_dense_field_struct
{
    struct _gr_tower_dense_field_struct * next;
    slong gid;
    slong degree;
    fmpz_poly_struct modulus;
    int ones;                     /* the modulus is 1 + x + ... + x^degree */
    const fmpz_poly_struct * src; /* the modulus of the flat tower it was last found for, */
    ulong version;                /* at this ideal version (a lookup without comparing moduli) */
}
_gr_tower_dense_field_struct;

void _gr_tower_dense_fields_clear(gr_tower_flat_t F);
/* (returns 0, leaving the tower unchanged, when the modulus cannot be used: its coefficients have denominators which cannot be rationalized, or the generator cannot be evaluated) */
int _gr_tower_make_algebraic(gr_tower_t T, slong d, const fmpz_mpoly_q_struct * m, slong len, const fmpz_mpoly_ctx_t mctx, int status);
int _gr_tower_flat_dense_applicable(gr_tower_flat_t F);
int _gr_tower_flat_mul_prefers_sparse(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_tower_flat_t F);
int _gr_tower_flat_mul_dense(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_tower_flat_t F);
int _gr_tower_flat_inv_dense(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t F);
void _gr_tower_flat_modp_clear(gr_tower_flat_t F);
int _gr_tower_flat_poly_mullow_dense(fmpz_mpoly_q_struct * const * res, fmpz_mpoly_q_struct * const * x, slong len1, fmpz_mpoly_q_struct * const * y, slong len2, slong n, gr_tower_flat_t F);
int _gr_tower_flat_mat_mul_dense(fmpz_mpoly_q_struct * const * C, fmpz_mpoly_q_struct * const * A, fmpz_mpoly_q_struct * const * B, slong r, slong s, slong c, gr_tower_flat_t F);
void _gr_tower_flat_transport(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t src, gr_tower_flat_t dst);
void _gr_tower_flat_transport_map(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t src, gr_tower_flat_t dst, const slong * order_map);
const fmpz_mpoly_struct * _gr_tower_flat_ideal_elem(gr_tower_flat_t F, slong k);
int _gr_tower_flat_find_layout(slong * cap, slong ** var_gid, const fmpz_mpoly_ctx_struct * mctx, gr_tower_flat_t F);

/* adds to the marked generators (by definition order) those their
   definitions involve, recursively */
void _gr_tower_involved_gens_closure(int * mark, gr_tower_t T);

/* marks the generators x involves, with those their definitions
   involve; whether these include a conjecturally transcendental one */
void _gr_tower_flat_involved(int * mark, const fmpz_mpoly_q_t x, gr_tower_flat_t F);
int _gr_tower_flat_involves_conjectural(const fmpz_mpoly_q_t x, gr_tower_flat_t F);

/* gr_tower_adjoin_root_ui with an enclosure zx of x (or NULL) */
int _gr_tower_adjoin_root_ui_enclosure(gr_tower_t T, gr_srcptr x, ulong n, const acb_t zx, int known, const char * name);
int _gr_tower_binomial_modular_evidence(gr_tower_t T, const fmpz_mpoly_q_t a, const fmpz_mpoly_ctx_t actx, ulong n, slong tries, ulong * power_prime);
slong _gr_tower_nested_num_coeffs(gr_srcptr x, gr_ctx_t ctx, slong limit);
char * _gr_tower_flat_get_str(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t mctx, gr_tower_t T);

/* Trial division of n by the first num_primes primes: 1 when complete;
   otherwise the cofactor is stored in cofactor (if not NULL) or appended
   to fac with exponent 1. */
int _gr_tower_fmpz_factor_trial(fmpz_factor_t fac, fmpz_t cofactor, const fmpz_t n, slong num_primes);

/* Whether the tower has generators whose transcendence is conjectural. */
int _gr_tower_has_conjectural(const gr_tower_t T);
int _gr_tower_has_conjectural_below(const gr_tower_t T, slong limit);

/* outcomes of the field zero test */
#define GR_TOWER_UNKNOWN 0
#define GR_TOWER_ZERO 1
#define GR_TOWER_NONZERO 2
#define GR_TOWER_FIELD_NONZERO 3

int _gr_tower_field_is_zero_at(gr_srcptr x, slong k, gr_tower_t T);
int _gr_tower_flat_nested_zero_code(const fmpz_mpoly_q_t x, gr_tower_flat_t F);

/* Certification that a nonzero element x of F_k (a nonzero rational
   function of the transcendental generators over the algebraic part)
   is a nonzero number, on a copy of the tower. */
truth_t _gr_tower_certify_nonzero_at(gr_srcptr x, slong k, gr_tower_t T);

/* Complete zero test of a flat element x (in the current context of the
   flat machinery F of a tower), searching for relations between the
   transcendental generators with definition order < limit and
   eliminating them in place (Richardson's algorithm). */
truth_t _gr_tower_decide_zero_flat(const fmpz_mpoly_q_t x, gr_tower_flat_t F, slong limit, slong depth);
truth_t _gr_tower_flat_num_is_zero_fixed(fmpz_mpoly_q_t x, gr_tower_flat_t F);

/* Searches for relations among all generators of F's tower at the given
   precision, eliminating those found; returns 1 if the tower changed. */
int _gr_tower_search_relations(gr_tower_flat_t F, slong prec);


/* Structured algebraic generators (roots of unity, roots of integers) */
int _gr_tower_gen_const_root(const gr_tower_gen_struct * g, fmpz_t p, ulong * n);
int _gr_tower_structural_proven(const gr_tower_gen_struct * g, const gr_tower_t U, slong limit);
int _gr_tower_unramified_radical(const gr_tower_gen_struct * g, const gr_tower_t U, slong limit, int require_proven);
void _gr_tower_express_downgraded(gr_tower_t T, slong from, slong to);
int _gr_tower_lazy_real_sign_locked(int * sign, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int _gr_tower_lazy_real_sign_fast(int * sgn, gr_srcptr x, gr_ctx_t ctx);
int _gr_tower_lazy_re_cmp(int * sgn, gr_srcptr x, const fmpq_t c, gr_ctx_t ctx);
int _gr_tower_lazy_im_sign(int * sgn, gr_srcptr x, gr_ctx_t ctx);
int _gr_tower_lazy_view_finish(int status, gr_ptr res, int real, int alg, gr_ctx_t ctx);
int _gr_tower_lazy_hurwitz_zeta_int(gr_ptr res, slong m, gr_srcptr a, gr_ctx_t ctx);
int _gr_tower_hypgeom_2f1_flags(const fmpz_mpoly_q_t a, const fmpz_mpoly_q_t b, const fmpz_mpoly_q_t c, const fmpz_mpoly_ctx_t mctx);
int _gr_tower_special_eval_multi_flags(acb_t res, int kind, slong param, acb_srcptr u, slong nargs, int flags, slong prec);
int _gr_tower_gen_hypgeom_flags(const gr_tower_gen_struct * g, gr_tower_flat_struct * F);
truth_t _gr_tower_lazy_is_real_exact(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
void _gr_tower_lazy_lock(gr_ctx_t ctx);
void _gr_tower_lazy_unlock(gr_ctx_t ctx);
int _gr_tower_lazy_outermost(gr_ctx_t ctx);
int _gr_tower_structured_image(fmpz_mpoly_q_t img, gr_tower_t U, const gr_tower_gen_struct * g, int gauss);
int _gr_tower_structured_present(gr_tower_t U, const gr_tower_gen_struct * g, int gauss);
int _gr_tower_structural_step_proven(gr_tower_t T, slong k);
int _gr_tower_gauss_sqrt_gen(gr_tower_t T, slong i);
int _gr_tower_cyclotomic_sqrt(fmpz_mpoly_q_t res, const fmpz_t c, gr_tower_flat_t F, slong limit);
int _gr_tower_all_proven(const gr_tower_t T);
void _gr_tower_gen_set_origin(gr_tower_gen_struct * g, const fmpz_poly_t p);
void _gr_tower_gen_rename(gr_tower_t T, gr_tower_gen_struct * g, const char * name);

int _gr_tower_get_acb_accurate(acb_t res, gr_srcptr x, slong prec, gr_tower_t T);
slong _gr_tower_select_fmpz_factor(const fmpz_poly_factor_t fac, int (*get_z)(acb_t, slong, void *), void * arg, gr_tower_t T);
slong _gr_tower_select_poly_factor(gr_vec_t fac, slong k, int (*get_z)(acb_t, slong, void *), void * arg, gr_tower_t T);
int _gr_tower_get_z_fixed(acb_t z, slong prec, void * arg);
int _gr_tower_select_candidate(int (*get)(acb_t, int, slong, void *), void * arg, gr_tower_t T);
int _gr_tower_principal_root_enclosure(acb_t res, gr_srcptr x, ulong n, const acb_t given, gr_tower_t T);

/* columns of the gamma relations at level q: x_k (as given by col_of_k),
   then log sin(pi k/q) (1 <= k < q/2), log n (2 <= n <= q; only primes
   are used), log pi */
#define GR_TOWER_GAMMA_CSIN(q) ((q) - 1)
#define GR_TOWER_GAMMA_CINT(q) ((q) - 1 + ((q) - 1) / 2)
#define GR_TOWER_GAMMA_CPI(q) (GR_TOWER_GAMMA_CINT(q) + (q) - 1)
#define GR_TOWER_GAMMA_NCOLS(q) (GR_TOWER_GAMMA_CPI(q) + 1)
slong _gr_tower_gamma_relations_rows(slong q);
slong _gr_tower_gamma_relations(fmpq_mat_t A, const slong * col_of_k, slong q);

/* columns of the Hurwitz zeta relations at level q (special.c): x_k (as
   given by col_of_k, 0 < k <= q), then E_k (1 <= k <= q/2), then Z */
#define GR_TOWER_HURWITZ_CE(q) (q)
#define GR_TOWER_HURWITZ_CZ(q) ((q) + (q) / 2)
#define GR_TOWER_HURWITZ_NCOLS(q) ((q) + (q) / 2 + 1)
slong _gr_tower_hurwitz_relations_rows(slong q);
slong _gr_tower_hurwitz_relations(fmpq_mat_t A, const slong * col_of_k, slong s, slong q);

/* the reduced row echelon form of the Hurwitz relations at weight s and
   level q (columns: the values zeta(s, k/q) in the order k_of_col, then
   E and Z as above), cached per thread */
typedef struct
{
    slong s, q, rank, ncols;
    slong * col_of_k;
    slong * k_of_col;
    slong * row_of_col;    /* the row with pivot in that column, or -1 */
    fmpq_mat_t B;
}
gr_tower_hurwitz_lattice_struct;

const gr_tower_hurwitz_lattice_struct * _gr_tower_hurwitz_lattice(slong s, slong q);
void _gr_tower_hurwitz_lattice_release(const gr_tower_hurwitz_lattice_struct * L);
slong _gr_tower_flat_max_dep(const fmpz_mpoly_q_t x, gr_tower_flat_t F);
int _gr_tower_flat_matches_gen(const fmpz_mpoly_q_t expr, slong d, gr_tower_flat_t F);
int _gr_tower_flat_matches_gen_bits(const fmpz_mpoly_q_t expr, slong d, slong bits, gr_tower_flat_t F);
int _gr_tower_flat_eliminate_gen(slong d, fmpz_mpoly_q_t expr, gr_tower_flat_t F);
void _gr_tower_cot_sum_cyclo(fmpq_poly_t E, const fmpq * ck, slong N, const fmpz_poly_t P, ulong m);
/* special.c: shared formulas of the special-function code */
void _gr_tower_cot_derivative_poly(fmpz_poly_t P, ulong j);
void _gr_tower_zeta_even_over_pi(fmpq_t res, ulong s);

/* certify.c helpers used by gamma_relations.c */
void _gr_tower_certify_mul_pow_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong e, gr_tower_flat_t F);
slong _gr_tower_certify_find_i(gr_tower_t T);
slong _gr_tower_certify_find_pi(gr_tower_t T, slong limit);
slong _gr_tower_certify_rational_radical_factor(fmpz_mpoly_q_struct * m, slong c, const fmpz_mpoly_q_t u, gr_tower_t T, slong j, gr_tower_flat_t F);
slong _gr_tower_certify_radical_factor(fmpz_mpoly_q_struct * m, slong c, const fmpz_mpoly_q_t u, gr_tower_t T, slong j, gr_tower_flat_t F);
slong _gr_tower_certify_lindep_all(fmpz_mat_t rel, acb_srcptr vec, slong n, slong prec);
int _gr_tower_gamma_round(gr_tower_flat_t F, slong limit, slong depth);
int _gr_tower_special_prepare_constants(gr_tower_t T, ulong m, int need_pi, slong before);
int _gr_tower_special_flat_cyclotomic_eval(fmpz_mpoly_q_t res, const fmpq_poly_t E, ulong m, gr_tower_flat_t F);
int _gr_tower_special_flat_zeta(fmpz_mpoly_q_t res, ulong m, gr_tower_flat_t F);
int _gr_tower_hurwitz_round(gr_tower_flat_t F, slong limit, slong depth);
int _gr_tower_hurwitz_line_round(gr_tower_flat_t F, slong limit, slong depth);
int _gr_tower_special_root_gen(slong * d, slong * mult, gr_tower_flat_t F, const fmpz_mpoly_q_t c, ulong delta, slong before);
int _gr_tower_elliptic_round(gr_tower_flat_t F, slong limit, slong depth);
int _gr_tower_special_trans_gen(slong * d, gr_tower_flat_t F, int kind, const fmpz_mpoly_q_t u, slong before);
int _gr_tower_dilog_round(gr_tower_flat_t F, slong limit, slong depth);
int _gr_tower_polylog_round(gr_tower_flat_t F, slong limit, slong depth);

#endif
