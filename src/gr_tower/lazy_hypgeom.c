/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Hypergeometric functions pFq(a; b; z) in the lazy fields.

    A value F(P; z), P = (a_1, ..., a_p; b_1, ..., b_q), is brought to a
    normal form in three layers:

    1. Degenerate parameters: cancellation of a_i = b_j, terminating
       series (a_i a nonpositive integer), poles (b_j a nonpositive
       integer), divergence (p > q + 1), z = 0, Gauss's sum at z = 1, and
       the reduction of the order when a_i - b_j is a nonnegative integer
       m: F = (theta + b_j)_m F_red / (b_j)_m, theta = z d/dz, where F_red
       omits a_i and b_j.

    2. Transformations of one term (Euler and Pfaff for 2F1, Kummer for
       1F1), which choose a canonical member of the orbit of (P, z).

    3. The contiguity module. With r = max(p, q + 1), the vector
       v(P) = (F, theta F, ..., theta^(r-1) F) satisfies
       v(P + e) = M v(P) for every unit shift e of a parameter, with M =
       I + C/a_i (raising a_i), I + C/(b_j - 1) (lowering b_j), where C is
       the companion matrix of theta modulo the hypergeometric operator
       theta prod (theta + b_j - 1) - z prod (theta + a_i). The shifts are
       invertible exactly when no a_i - b_j is an integer (det(I + C/a_i)
       is a multiple of prod_j (b_j - 1 - a_i)), so the values in the
       coset P + Z^(p+q) are rational combinations of the values of the r
       basis functions G_k = F(P0 + k (1, ..., 1)) at a canonical cell P0
       (parameters shifted to 0 < Re <= 1). When b_j - a_i is a positive
       integer (a reducible operator), the cell keeps b_j - a_i = 1 and
       the paths stay in the chamber b_j - a_i >= 1, where the shifts are
       invertible, and (theta + a_i) F = a_i F_red is a known linear
       relation between the basis values.

    Known values give more linear relations between the basis values: the
    identities of the table below, matched against every point of the
    coset (and of its images under the transformations of layer 2). The
    relations are reduced (preferring the basis values of the lowest
    index as the remaining generators); a basis value determined by them
    is replaced by its closed form, and the others become generators of
    kind GR_TOWER_HYPGEOM. Thus a single entry such as 0F1(; 1/2; z) =
    cosh(2 sqrt(z)) with 0F1(; 3/2; z) = sinh(2 sqrt(z)) / (2 sqrt(z))
    gives all of 0F1(; n + 1/2; z), and K, E give all of 2F1 with
    parameters in (1/2, 1/2, 1) + Z^3.

    The table is declarative: each entry is a pattern (affine parameters
    in symbols, a constant or free argument) and a value, written in a
    small expression language (see _hyp_eval_expr), with optional
    conditions of the form "re: expr" (Re(expr) > 0). The complete
    elliptic integrals K and E are entries of the table (and
    gr_tower_lazy_elliptic_k, gr_tower_lazy_elliptic_e are evaluated
    through it, so that their argument transformations and special values
    are those of 2F1: the imaginary-modulus transformation is Pfaff's,
    E(1) and K(1/2), E(1/2) are Gauss's sums).
*/

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "fmpq.h"
#include "fmpq_vec.h"
#include "arith.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

#define FLAGS(ctx) gr_tower_lazy_ctx_field_flags(ctx)
#define ALG(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_ALGEBRAIC) && _gr_tower_lazy_outermost(ctx))
#define REAL(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_REAL) && _gr_tower_lazy_outermost(ctx))

#define CHECK_PREC 64
/* nesting of the evaluation (order reduction, table entries in terms of
   other hypergeometric functions) */
#define HYP_DEPTH_LIMIT 8
/* terminating series up to this many terms */
#define HYP_TERM_LIMIT 2000
/* the reduction of the order for a_i - b_j up to this */
#define HYP_REDUCE_LIMIT 64
/* shifts of the parameters (total) along a path in the contiguity module */
#define HYP_SHIFT_LIMIT 200
/* the points of a coset matching a table entry: shifts of the pattern
   symbols in [-HYP_BOX, HYP_BOX] */
#define HYP_BOX 2
/* the largest order r handled by the contiguity module */
#define HYP_ORDER_LIMIT 6

#define ENTRY(v, i, sz) GR_ENTRY(v, i, sz)

/* -------------------------------------------------------------------- */
/* exact tests on lazy elements                                          */
/* -------------------------------------------------------------------- */

/* Whether x is a rational number (exactly), setting c. Rational
   representations first; otherwise a rational of small height close to
   x is verified exactly. Returns 1, 0 (not rational, or not known to be
   one), or -1 (unknown). */
int
_gr_tower_lazy_rational_recognize(fmpq_t c, gr_srcptr x, gr_ctx_t ctx)
{
    acb_t z;
    int res = 0;

    if (_gr_tower_lazy_rational_repr_locked(c, x, ctx))
        return 1;

    acb_init(z);
    if (gr_tower_lazy_get_acb(z, x, CHECK_PREC, ctx) == GR_SUCCESS &&
        arb_contains_zero(acb_imagref(z)) && arb_is_finite(acb_realref(z)) &&
        mag_cmp_2exp_si(arb_radref(acb_realref(z)), -20) < 0)
    {
        /* a candidate p/q with q <= 1000 close to x */
        arb_t t;
        fmpz_t n;
        slong q;
        arb_init(t);
        fmpz_init(n);
        for (q = 1; q <= 1000; q++)
        {
            arb_mul_ui(t, acb_realref(z), q, CHECK_PREC);
            arf_get_fmpz(n, arb_midref(t), ARF_RND_NEAR);
            if (arb_contains_fmpz(t, n))
            {
                gr_ptr u;
                truth_t eq;
                fmpz_set(fmpq_numref(c), n);
                fmpz_set_si(fmpq_denref(c), q);
                fmpq_canonicalise(c);
                GR_TMP_INIT(u, ctx);
                if (gr_set_fmpq(u, c, ctx) == GR_SUCCESS)
                {
                    eq = gr_equal(u, x, ctx);
                    if (eq == T_TRUE)
                        res = 1;
                    else if (eq == T_UNKNOWN)
                        res = -1;
                }
                GR_TMP_CLEAR(u, ctx);
                break;
            }
        }
        fmpz_clear(n);
        arb_clear(t);
    }
    acb_clear(z);
    return res;
}

#define _hyp_rational _gr_tower_lazy_rational_recognize

/* Whether x is an integer (then sets n, if it fits): 1, 0 (certainly not),
   -1 (unknown). */
static int
_hyp_integer(slong * n, gr_srcptr x, gr_ctx_t ctx)
{
    acb_t z;
    int res;

    acb_init(z);
    if (gr_tower_lazy_get_acb(z, x, CHECK_PREC, ctx) == GR_SUCCESS &&
        (!arb_contains_zero(acb_imagref(z)) || !arb_contains_int(acb_realref(z))))
    {
        acb_clear(z);
        return 0;
    }
    acb_clear(z);

    {
        fmpq_t c;
        fmpq_init(c);
        res = _hyp_rational(c, x, ctx);
        if (res == 1)
        {
            if (fmpz_is_one(fmpq_denref(c)) && fmpz_fits_si(fmpq_numref(c)))
                *n = fmpz_get_si(fmpq_numref(c));
            else if (fmpz_is_one(fmpq_denref(c)))
                res = -1;   /* (too large) */
            else
                res = 0;
        }
        else if (res == 0)
            res = -1;   /* (not known to be rational) */
        fmpq_clear(c);
    }

    return res;
}

/* n = the integer with 0 < Re(x) - n <= 1 */
static int
_hyp_cell_shift(slong * n, gr_srcptr x, gr_ctx_t ctx)
{
    acb_t z;
    arb_t f;
    fmpz_t a;
    int status;

    acb_init(z);
    arb_init(f);
    fmpz_init(a);

    status = gr_tower_lazy_get_acb(z, x, CHECK_PREC, ctx);
    if (status == GR_SUCCESS)
    {
        /* ceil(Re(x)) - 1 when Re(x) is not an integer */
        if (!arb_contains_int(acb_realref(z)) && arb_is_finite(acb_realref(z)))
        {
            arb_ceil(f, acb_realref(z), CHECK_PREC);
            if (arb_is_exact(f) && arf_is_int(arb_midref(f)))
            {
                arf_get_fmpz(a, arb_midref(f), ARF_RND_NEAR);
                if (fmpz_fits_si(a))
                {
                    *n = fmpz_get_si(a) - 1;
                    goto cleanup;
                }
            }
            status = GR_UNABLE;
        }
        else
        {
            /* Re(x) may be an integer: exactly */
            gr_ptr t;
            slong m;
            GR_TMP_INIT(t, ctx);
            status = gr_re(t, x, ctx);
            if (status == GR_SUCCESS)
            {
                int r = _hyp_integer(&m, t, ctx);
                if (r == 1)
                    *n = m - 1;
                else if (r == 0)
                {
                    /* Re(x) is not an integer although its enclosure
                       contains one: more precision */
                    acb_t w;
                    slong prec;
                    acb_init(w);
                    status = GR_UNABLE;
                    for (prec = 2 * CHECK_PREC; prec <= 4096; prec *= 2)
                    {
                        if (gr_tower_lazy_get_acb(w, t, prec, ctx) != GR_SUCCESS)
                            break;
                        if (!arb_contains_int(acb_realref(w)))
                        {
                            arb_ceil(f, acb_realref(w), prec);
                            if (arb_is_exact(f))
                            {
                                arf_get_fmpz(a, arb_midref(f), ARF_RND_NEAR);
                                if (fmpz_fits_si(a))
                                {
                                    *n = fmpz_get_si(a) - 1;
                                    status = GR_SUCCESS;
                                }
                            }
                            break;
                        }
                    }
                    acb_clear(w);
                }
                else
                    status = GR_UNABLE;
            }
            GR_TMP_CLEAR(t, ctx);
        }
    }

cleanup:
    acb_clear(z);
    arb_clear(f);
    fmpz_clear(a);
    return status;
}

/* res = x^e (principal), with e = e0 + n, 0 < Re(e0) <= 1, as x^e0 x^n
   for an irrational e, so that the powers at exponents differing by
   integers (contiguous parameters) share one generator */
static int
_hyp_pow(gr_ptr res, gr_srcptr x, gr_srcptr e, gr_ctx_t ctx)
{
    fmpq_t c;
    slong n;
    int status;

    fmpq_init(c);
    if (_gr_tower_lazy_rational_repr_locked(c, e, ctx))
    {
        fmpq_t xq;
        fmpq_init(xq);
        if (fmpq_sgn(c) < 0 && _gr_tower_lazy_rational_repr_locked(xq, x, ctx) && fmpq_sgn(xq) > 0)
        {
            /* (1/x)^(-c) for a positive rational x: the inverse of the
               base rather than of a power in a field of radicals (not on
               the negative axis, where the branches differ) */
            gr_ptr t;
            GR_TMP_INIT(t, ctx);
            fmpq_neg(c, c);
            fmpq_inv(xq, xq);
            status = gr_set_fmpq(t, xq, ctx);
            if (status == GR_SUCCESS)
                status = gr_pow_fmpq(res, t, c, ctx);
            GR_TMP_CLEAR(t, ctx);
        }
        else
            status = gr_pow_fmpq(res, x, c, ctx);
        fmpq_clear(xq);
    }
    else if (_hyp_cell_shift(&n, e, ctx) == GR_SUCCESS && n != 0)
    {
        gr_ptr t, u;
        GR_TMP_INIT2(t, u, ctx);
        status = gr_sub_si(t, e, n, ctx);
        status |= gr_pow(u, x, t, ctx);
        status |= gr_pow_si(t, x, n, ctx);
        status |= gr_mul(res, u, t, ctx);
        GR_TMP_CLEAR2(t, u, ctx);
    }
    else
        status = gr_pow(res, x, e, ctx);
    fmpq_clear(c);
    return status;
}

/* Numerical comparison of x and y (real parts, then imaginary parts) for
   a canonical order; 0 if they are equal (exactly) or the order cannot be
   decided (then the caller keeps its order). */
static int
_hyp_cmp(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    acb_t u, v;
    slong prec;
    int res = 0;

    if (gr_equal(x, y, ctx) == T_TRUE)
        return 0;

    acb_init(u);
    acb_init(v);
    for (prec = CHECK_PREC; prec <= 1024 && res == 0; prec *= 4)
    {
        if (gr_tower_lazy_get_acb(u, x, prec, ctx) != GR_SUCCESS ||
            gr_tower_lazy_get_acb(v, y, prec, ctx) != GR_SUCCESS)
            break;
        if (arb_lt(acb_realref(u), acb_realref(v)))
            res = -1;
        else if (arb_gt(acb_realref(u), acb_realref(v)))
            res = 1;
        else if (arb_overlaps(acb_realref(u), acb_realref(v)) &&
                 (arb_is_exact(acb_realref(u)) && arb_is_exact(acb_realref(v))))
        {
            if (arb_lt(acb_imagref(u), acb_imagref(v)))
                res = -1;
            else if (arb_gt(acb_imagref(u), acb_imagref(v)))
                res = 1;
        }
        else if (arb_overlaps(acb_realref(u), acb_realref(v)))
        {
            /* the real parts may be equal: compare them exactly */
            gr_ptr s, t;
            GR_TMP_INIT2(s, t, ctx);
            if (gr_re(s, x, ctx) == GR_SUCCESS && gr_re(t, y, ctx) == GR_SUCCESS &&
                gr_equal(s, t, ctx) == T_TRUE)
            {
                if (arb_lt(acb_imagref(u), acb_imagref(v)))
                    res = -1;
                else if (arb_gt(acb_imagref(u), acb_imagref(v)))
                    res = 1;
            }
            GR_TMP_CLEAR2(s, t, ctx);
        }
    }
    acb_clear(u);
    acb_clear(v);
    return res;
}

/* sorts v (n elements) by _hyp_cmp (insertion sort: n is small) */
static void
_hyp_sort(gr_ptr v, slong n, gr_ctx_t ctx)
{
    slong i, j, sz = ctx->sizeof_elem;
    for (i = 1; i < n; i++)
        for (j = i; j > 0 && _hyp_cmp(ENTRY(v, j - 1, sz), ENTRY(v, j, sz), ctx) > 0; j--)
            gr_swap(ENTRY(v, j - 1, sz), ENTRY(v, j, sz), ctx);
}

/* -------------------------------------------------------------------- */
/* parameter vectors                                                     */
/* -------------------------------------------------------------------- */

/* P = (a_1, ..., a_p; b_1, ..., b_q): stored as one vector of p + q
   elements, a first */
typedef struct
{
    slong p, q;
    gr_ptr v;
}
hyp_params_struct;

typedef hyp_params_struct hyp_params_t[1];

#define HA(P, i, sz) ENTRY((P)->v, (i), (sz))
#define HB(P, j, sz) ENTRY((P)->v, (P)->p + (j), (sz))

static void
_hyp_params_init(hyp_params_t P, slong p, slong q, gr_ctx_t ctx)
{
    P->p = p;
    P->q = q;
    P->v = gr_heap_init_vec(FLINT_MAX(p + q, 1), ctx);
}

static void
_hyp_params_clear(hyp_params_t P, gr_ctx_t ctx)
{
    gr_heap_clear_vec(P->v, FLINT_MAX(P->p + P->q, 1), ctx);
}

static int
_hyp_params_set(hyp_params_t R, const hyp_params_t P, gr_ctx_t ctx)
{
    if (R == P)
        return GR_SUCCESS;
    _hyp_params_clear(R, ctx);
    _hyp_params_init(R, P->p, P->q, ctx);
    return _gr_vec_set(R->v, P->v, P->p + P->q, ctx);
}

/* R = P without a_i and b_j (i or j = -1: none) */
static int
_hyp_params_drop(hyp_params_t R, const hyp_params_t P, slong i, slong j, gr_ctx_t ctx)
{
    slong k, n = 0, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    hyp_params_t S;

    _hyp_params_init(S, P->p - (i >= 0), P->q - (j >= 0), ctx);
    for (k = 0; k < P->p; k++)
        if (k != i)
            status |= gr_set(ENTRY(S->v, n++, sz), HA(P, k, sz), ctx);
    for (k = 0; k < P->q; k++)
        if (k != j)
            status |= gr_set(ENTRY(S->v, n++, sz), HB(P, k, sz), ctx);
    _hyp_params_clear(R, ctx);
    *R = *S;
    return status;
}

/* sorts the a and the b separately */
static void
_hyp_params_sort(hyp_params_t P, gr_ctx_t ctx)
{
    _hyp_sort(P->v, P->p, ctx);
    _hyp_sort(ENTRY(P->v, P->p, ctx->sizeof_elem), P->q, ctx);
}

/* -------------------------------------------------------------------- */
/* forward declarations                                                  */
/* -------------------------------------------------------------------- */

static int _hyp_eval(gr_ptr res, const hyp_params_t P, gr_srcptr z, int depth, int flags, gr_ctx_t ctx);

/* flags of _hyp_eval */
#define HYP_NO_ORBIT 1      /* the transformations of layer 2 were applied already */

/* -------------------------------------------------------------------- */
/* expressions of the table                                              */
/* -------------------------------------------------------------------- */

/*
    The values of the table entries are written in a small expression
    language: integers, the symbols of the pattern (single letters, z for
    the argument), + - * / ^, parentheses, pi, and the functions sqrt,
    exp, log, sin, cos, sinh, cosh, asin, atanh, gamma, rgamma, erf,
    elliptic_k_gen, elliptic_e_gen (the generators K(m), E(m) with no
    hypergeometric transformation: see _elliptic_gen in lazy_special.c),
    hyp0f1(b, z), hyp1f1(a, b, z), hyp2f1(a, b, c, z).
*/

/* the symbols a, ..., z */
#define HYP_NVARS_MAX 26

typedef struct
{
    const char * s;
    gr_srcptr * vars;     /* vars[c - 'a'] for the symbols, NULL if unbound */
    int depth;
    int numeric;          /* evaluated in a context of complex balls (arithmetic and sqrt only) */
    gr_ctx_struct * ctx;
}
hyp_expr_struct;

static int _hyp_ev_sum(gr_ptr res, hyp_expr_struct * E);

static void
_hyp_ev_space(hyp_expr_struct * E)
{
    while (*E->s == ' ')
        E->s++;
}

static int
_hyp_ev_args(gr_ptr args, slong n, hyp_expr_struct * E)
{
    slong i, sz = E->ctx->sizeof_elem;
    int status = GR_SUCCESS;

    _hyp_ev_space(E);
    if (*E->s != '(')
        return GR_UNABLE;
    E->s++;
    for (i = 0; i < n && status == GR_SUCCESS; i++)
    {
        status = _hyp_ev_sum(ENTRY(args, i, sz), E);
        _hyp_ev_space(E);
        if (status == GR_SUCCESS)
        {
            if (*E->s == ((i == n - 1) ? ')' : ','))
                E->s++;
            else
                status = GR_UNABLE;
        }
    }
    return status;
}

/* whether b is a nonpositive integer */
static int
_hyp_ev_lower_pole(gr_srcptr b, gr_ctx_t ctx)
{
    slong n;
    return gr_get_si(&n, b, ctx) == GR_SUCCESS && n <= 0;
}

static int
_hyp_ev_atom(gr_ptr res, hyp_expr_struct * E)
{
    gr_ctx_struct * ctx = E->ctx;
    int status = GR_SUCCESS;

    _hyp_ev_space(E);

    if (*E->s == '(')
    {
        E->s++;
        status = _hyp_ev_sum(res, E);
        _hyp_ev_space(E);
        if (*E->s != ')')
            return GR_UNABLE;
        E->s++;
        return status;
    }

    if (isdigit((unsigned char) *E->s))
    {
        slong n = 0;
        while (isdigit((unsigned char) *E->s))
            n = 10 * n + (*E->s++ - '0');
        return gr_set_si(res, n, ctx);
    }

    if (isalpha((unsigned char) *E->s))
    {
        char name[32];
        slong len = 0;

        while ((isalnum((unsigned char) *E->s) || *E->s == '_') && len < 31)
            name[len++] = *E->s++;
        name[len] = '\0';

        if (len == 1 && name[0] >= 'a' && name[0] <= 'z')
        {
            gr_srcptr v = E->vars[name[0] - 'a'];
            return (v == NULL) ? GR_UNABLE : gr_set(res, v, ctx);
        }

        if (strcmp(name, "pi") == 0)
            return gr_pi(res, ctx);

        {
            gr_ptr args;
            slong n = (strcmp(name, "hyp2f1") == 0) ? 4 : (strcmp(name, "hyp1f1") == 0) ? 3 :
                      (strcmp(name, "hyp0f1") == 0) ? 2 : 1;
            slong sz = ctx->sizeof_elem;

            args = gr_heap_init_vec(n, ctx);
            status = _hyp_ev_args(args, n, E);

            if (status == GR_SUCCESS)
            {
                if (strcmp(name, "sqrt") == 0) status = gr_sqrt(res, args, ctx);
                else if (strcmp(name, "exp") == 0) status = gr_exp(res, args, ctx);
                else if (strcmp(name, "log") == 0) status = gr_log(res, args, ctx);
                else if (strcmp(name, "sin") == 0) status = gr_sin(res, args, ctx);
                else if (strcmp(name, "cos") == 0) status = gr_cos(res, args, ctx);
                else if (strcmp(name, "sinh") == 0) status = gr_sinh(res, args, ctx);
                else if (strcmp(name, "cosh") == 0) status = gr_cosh(res, args, ctx);
                else if (strcmp(name, "asin") == 0) status = gr_asin(res, args, ctx);
                else if (strcmp(name, "atanh") == 0) status = gr_atanh(res, args, ctx);
                else if (strcmp(name, "gamma") == 0) status = gr_gamma(res, args, ctx);
                else if (strcmp(name, "rgamma") == 0) status = gr_rgamma(res, args, ctx);
                else if (strcmp(name, "erf") == 0) status = gr_erf(res, args, ctx);
                else if (strcmp(name, "elliptic_k_gen") == 0) status = _gr_tower_lazy_elliptic_gen(res, args, GR_TOWER_ELLIPTIC_K, ctx);
                else if (strcmp(name, "elliptic_e_gen") == 0) status = _gr_tower_lazy_elliptic_gen(res, args, GR_TOWER_ELLIPTIC_E, ctx);
                else if (n >= 2 && _hyp_ev_lower_pole(ENTRY(args, n - 2, sz), ctx))
                {
                    /* the identities hold by continuation in the
                       parameters: at a nonpositive integer lower
                       parameter the function on the right is the limit,
                       not the terminating series of the convention
                       (2F1(1/4, 3/4; -1/2; z) through 2F1(1/2, -1; -2; w),
                       with c = -1/2 in entry "a, a+1/2 | c") */
                    status = GR_UNABLE;
                }
                else if (n >= 2)
                {
                    /* a hypergeometric function of another family */
                    hyp_params_t P;
                    slong i;
                    _hyp_params_init(P, n - 2, 1, ctx);
                    for (i = 0; i < n - 1; i++)
                        status |= gr_set(ENTRY(P->v, i, sz), ENTRY(args, i, sz), ctx);
                    if (status == GR_SUCCESS)
                        status = _hyp_eval(res, P, ENTRY(args, n - 1, sz), E->depth + 1, 0, ctx);
                    _hyp_params_clear(P, ctx);
                }
                else
                    status = GR_UNABLE;
            }

            gr_heap_clear_vec(args, n, ctx);
            return status;
        }
    }

    return GR_UNABLE;
}

static int
_hyp_ev_power(gr_ptr res, hyp_expr_struct * E)
{
    int status = _hyp_ev_atom(res, E);
    _hyp_ev_space(E);
    if (status == GR_SUCCESS && *E->s == '^')
    {
        gr_ptr t;
        fmpq_t c;
        E->s++;
        GR_TMP_INIT(t, E->ctx);
        _hyp_ev_space(E);
        if (*E->s == '-')
        {
            E->s++;
            status = _hyp_ev_power(t, E);
            status |= gr_neg(t, t, E->ctx);
        }
        else
            status = _hyp_ev_power(t, E);
        fmpq_init(c);
        if (status == GR_SUCCESS && E->numeric)
        {
            slong n;
            if (gr_get_si(&n, t, E->ctx) == GR_SUCCESS)
                status = gr_pow_si(res, res, n, E->ctx);
            else
                status = gr_pow(res, res, t, E->ctx);
        }
        else if (status == GR_SUCCESS)
        {
            /* (rational exponents: roots, powers) */
            if (_gr_tower_lazy_rational_repr_locked(c, t, E->ctx))
                status = gr_pow_fmpq(res, res, c, E->ctx);
            else
                status = gr_pow(res, res, t, E->ctx);
        }
        fmpq_clear(c);
        GR_TMP_CLEAR(t, E->ctx);
    }
    return status;
}

static int
_hyp_ev_unary(gr_ptr res, hyp_expr_struct * E)
{
    _hyp_ev_space(E);
    if (*E->s == '-')
    {
        int status;
        E->s++;
        status = _hyp_ev_unary(res, E);
        return status | gr_neg(res, res, E->ctx);
    }
    return _hyp_ev_power(res, E);
}

static int
_hyp_ev_product(gr_ptr res, hyp_expr_struct * E)
{
    int status = _hyp_ev_unary(res, E);
    gr_ptr t;

    GR_TMP_INIT(t, E->ctx);
    for (;;)
    {
        char op;
        _hyp_ev_space(E);
        op = *E->s;
        if (status != GR_SUCCESS || (op != '*' && op != '/'))
            break;
        E->s++;
        status = _hyp_ev_unary(t, E);
        if (status == GR_SUCCESS)
            status = (op == '*') ? gr_mul(res, res, t, E->ctx) : gr_div(res, res, t, E->ctx);
    }
    GR_TMP_CLEAR(t, E->ctx);
    return status;
}

static int
_hyp_ev_sum(gr_ptr res, hyp_expr_struct * E)
{
    int status = _hyp_ev_product(res, E);
    gr_ptr t;

    GR_TMP_INIT(t, E->ctx);
    for (;;)
    {
        char op;
        _hyp_ev_space(E);
        op = *E->s;
        if (status != GR_SUCCESS || (op != '+' && op != '-'))
            break;
        E->s++;
        status = _hyp_ev_product(t, E);
        if (status == GR_SUCCESS)
            status = (op == '+') ? gr_add(res, res, t, E->ctx) : gr_sub(res, res, t, E->ctx);
    }
    GR_TMP_CLEAR(t, E->ctx);
    return status;
}

/* evaluates the expression s with the given symbols */
static int
_hyp_eval_expr(gr_ptr res, const char * s, gr_srcptr * vars, int depth, gr_ctx_t ctx)
{
    hyp_expr_struct E;
    int status;

    E.s = s;
    E.vars = vars;
    E.depth = depth;
    E.numeric = 0;
    E.ctx = ctx;

    status = _hyp_ev_sum(res, &E);
    _hyp_ev_space(&E);
    if (status == GR_SUCCESS && *E.s != '\0')
        status = GR_UNABLE;
    return status;
}

/* the expression s (of len characters; arithmetic and sqrt), numerically,
   with the symbols bound in vars */
static int
_hyp_eval_expr_acb(acb_t res, const char * s, slong len, gr_srcptr * vars, gr_ctx_t ctx)
{
    gr_ctx_t CC;
    hyp_expr_struct E;
    acb_ptr v;
    gr_srcptr nvars[HYP_NVARS_MAX];
    char * t;
    slong i;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, 128);
    v = _acb_vec_init(HYP_NVARS_MAX);
    for (i = 0; i < HYP_NVARS_MAX && status == GR_SUCCESS; i++)
    {
        nvars[i] = NULL;
        if (vars[i] != NULL)
        {
            status = gr_tower_lazy_get_acb(v + i, vars[i], 128, ctx);
            nvars[i] = v + i;
        }
    }

    t = flint_malloc(len + 1);
    memcpy(t, s, len);
    t[len] = '\0';
    if (status == GR_SUCCESS)
    {
        E.s = t;
        E.vars = nvars;
        E.depth = 0;
        E.numeric = 1;
        E.ctx = CC;
        status = _hyp_ev_sum(res, &E);
        _hyp_ev_space(&E);
        if (status == GR_SUCCESS && *E.s != '\0')
            status = GR_UNABLE;
    }
    flint_free(t);
    _acb_vec_clear(v, HYP_NVARS_MAX);
    gr_ctx_clear(CC);
    return status;
}

/* -------------------------------------------------------------------- */
/* affine patterns                                                       */
/* -------------------------------------------------------------------- */

/* an affine form c_0 + sum_v c_v v in the symbols a, ..., y (v = 1..25) */
#define HYP_NSYM 26
/* the symbols a, ..., z of the value expressions (z: the argument) */
#define HYP_NVARS 26

typedef struct
{
    fmpq c[HYP_NSYM];    /* c[0]: the constant; c[1 + k]: the symbol 'a' + k */
}
hyp_affine_struct;

static void
_hyp_affine_init(hyp_affine_struct * A)
{
    slong i;
    for (i = 0; i < HYP_NSYM; i++)
        fmpq_init(A->c + i);
}

static void
_hyp_affine_clear(hyp_affine_struct * A)
{
    slong i;
    for (i = 0; i < HYP_NSYM; i++)
        fmpq_clear(A->c + i);
}

static int _hyp_aff_sum(hyp_affine_struct * A, const char ** s);

static void
_hyp_aff_space(const char ** s)
{
    while (**s == ' ')
        (*s)++;
}

static int
_hyp_aff_is_const(const hyp_affine_struct * A)
{
    slong i;
    for (i = 1; i < HYP_NSYM; i++)
        if (!fmpq_is_zero(A->c + i))
            return 0;
    return 1;
}

static int
_hyp_aff_atom(hyp_affine_struct * A, const char ** s)
{
    slong i;

    _hyp_aff_space(s);
    for (i = 0; i < HYP_NSYM; i++)
        fmpq_zero(A->c + i);

    if (**s == '(')
    {
        (*s)++;
        if (!_hyp_aff_sum(A, s))
            return 0;
        _hyp_aff_space(s);
        if (**s != ')')
            return 0;
        (*s)++;
        return 1;
    }
    if (isdigit((unsigned char) **s))
    {
        slong n = 0;
        while (isdigit((unsigned char) **s))
            n = 10 * n + (*(*s)++ - '0');
        fmpq_set_si(A->c, n, 1);
        return 1;
    }
    if (islower((unsigned char) **s) && **s != 'z' && !isalnum((unsigned char) (*s)[1]))
    {
        fmpq_one(A->c + 1 + (**s - 'a'));
        (*s)++;
        return 1;
    }
    return 0;
}

static int
_hyp_aff_unary(hyp_affine_struct * A, const char ** s)
{
    _hyp_aff_space(s);
    if (**s == '-')
    {
        slong i;
        (*s)++;
        if (!_hyp_aff_unary(A, s))
            return 0;
        for (i = 0; i < HYP_NSYM; i++)
            fmpq_neg(A->c + i, A->c + i);
        return 1;
    }
    return _hyp_aff_atom(A, s);
}

static int
_hyp_aff_product(hyp_affine_struct * A, const char ** s)
{
    hyp_affine_struct B;
    int ok = _hyp_aff_unary(A, s);
    slong i;

    _hyp_affine_init(&B);
    for (;;)
    {
        char op;
        _hyp_aff_space(s);
        op = **s;
        if (!ok || (op != '*' && op != '/'))
            break;
        (*s)++;
        ok = _hyp_aff_unary(&B, s);
        if (!ok)
            break;
        if (op == '*')
        {
            /* one side must be constant */
            if (_hyp_aff_is_const(&B))
            {
                for (i = 0; i < HYP_NSYM; i++)
                    fmpq_mul(A->c + i, A->c + i, B.c);
            }
            else if (_hyp_aff_is_const(A))
            {
                fmpq_t k;
                fmpq_init(k);
                fmpq_set(k, A->c);
                for (i = 0; i < HYP_NSYM; i++)
                    fmpq_mul(A->c + i, B.c + i, k);
                fmpq_clear(k);
            }
            else
                ok = 0;
        }
        else
        {
            if (!_hyp_aff_is_const(&B) || fmpq_is_zero(B.c))
                ok = 0;
            else
                for (i = 0; i < HYP_NSYM; i++)
                    fmpq_div(A->c + i, A->c + i, B.c);
        }
    }
    _hyp_affine_clear(&B);
    return ok;
}

static int
_hyp_aff_sum(hyp_affine_struct * A, const char ** s)
{
    hyp_affine_struct B;
    int ok = _hyp_aff_product(A, s);
    slong i;

    _hyp_affine_init(&B);
    for (;;)
    {
        char op;
        _hyp_aff_space(s);
        op = **s;
        if (!ok || (op != '+' && op != '-'))
            break;
        (*s)++;
        ok = _hyp_aff_product(&B, s);
        if (ok)
            for (i = 0; i < HYP_NSYM; i++)
            {
                if (op == '+')
                    fmpq_add(A->c + i, A->c + i, B.c + i);
                else
                    fmpq_sub(A->c + i, A->c + i, B.c + i);
            }
    }
    _hyp_affine_clear(&B);
    return ok;
}

/* the value of A at the symbols vals (all those occurring bound) */
static int
_hyp_aff_eval(gr_ptr res, const hyp_affine_struct * A, gr_srcptr * vals, gr_ctx_t ctx)
{
    gr_ptr t;
    slong i;
    int status;

    GR_TMP_INIT(t, ctx);
    status = gr_set_fmpq(res, A->c, ctx);
    for (i = 1; i < HYP_NSYM && status == GR_SUCCESS; i++)
    {
        if (fmpq_is_zero(A->c + i))
            continue;
        if (vals[i - 1] == NULL)
            status = GR_UNABLE;
        else
        {
            status = gr_mul_fmpq(t, vals[i - 1], A->c + i, ctx);
            status |= gr_add(res, res, t, ctx);
        }
    }
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* the table                                                             */
/* -------------------------------------------------------------------- */

/*
    Each entry: "p q | a_1, ..., a_p | b_1, ..., b_q | z | value" with an
    optional " | conditions", comma-separated: "re: expr" (Re(expr) > 0),
    "ne: expr" (expr != 0), "disc: expr" (|expr| < 1), "cut: expr" (expr
    off the cut [1, inf)), and "anchor: expr", which makes the entry a
    link: it is used only when a generator 2F1, K or E of the context is
    at the argument expr (numerically, up to Pfaff's transformation), and
    not inside the evaluation of another link.
    The z field is z (free) or a rational constant. The parameters are
    affine in the symbols, each symbol occurring alone (with coefficient
    +-1) in some parameter. Every entry is checked numerically by the
    tests (t-hypgeom).
*/
static const char * const _hyp_table[] = {
    /* elementary 0F1 */
    "0 1 | | 1/2 | z | cosh(2*sqrt(z))",
    "0 1 | | 3/2 | z | sinh(2*sqrt(z))/(2*sqrt(z))",

    /* 1F1: the error function (b - a = 1: with the reduction
       (theta + 1/2) F = exp(z)/2, all of 1F1(n + 1/2; m + 3/2; z),
       m >= n); Kummer's second formula (Bessel functions as 0F1) */
    "1 1 | 1/2 | 3/2 | z | sqrt(pi)*erf(sqrt(-z))/(2*sqrt(-z))",
    "1 1 | a | 2*a | z | exp(z/2)*hyp0f1(a+1/2, z^2/16)",

    /* 2F1: complete elliptic integrals */
    "2 1 | 1/2, 1/2 | 1 | z | 2*elliptic_k_gen(z)/pi",
    "2 1 | -1/2, 1/2 | 1 | z | 2*elliptic_e_gen(z)/pi",

    /* 2F1: elementary (z off the cut [1, inf), principal branches) */
    "2 1 | 1, 1 | 2 | z | -log(1-z)/z",
    "2 1 | 1/2, 1 | 3/2 | z | atanh(sqrt(z))/sqrt(z)",
    "2 1 | 1/2, 1/2 | 3/2 | z | asin(sqrt(z))/sqrt(z)",
    "2 1 | a, a+1/2 | 1/2 | z | ((1+sqrt(z))^(-2*a)+(1-sqrt(z))^(-2*a))/2",
    "2 1 | a, a+1/2 | 3/2 | z | ((1+sqrt(z))^(1-2*a)-(1-sqrt(z))^(1-2*a))/(2*(1-2*a)*sqrt(z)) | ne: 1-2*a",

    /* 2F1: summation theorems */
    "2 1 | a, b | 1+a-b | -1 | gamma(1+a-b)*gamma(1+a/2)/(gamma(1+a)*gamma(1+a/2-b))",
    "2 1 | a, b | (a+b+1)/2 | 1/2 | sqrt(pi)*gamma((a+b+1)/2)/(gamma((a+1)/2)*gamma((b+1)/2))",
    "2 1 | a, 1-a | c | 1/2 | gamma(c/2)*gamma((1+c)/2)/(gamma((c+a)/2)*gamma((1+c-a)/2))",

    /* 3F2 at 1: Dixon, Watson, Whipple */
    "3 2 | a, b, c | 1+a-b, 1+a-c | 1 | gamma(1+a/2)*gamma(1+a-b)*gamma(1+a-c)*gamma(1+a/2-b-c)/(gamma(1+a)*gamma(1+a/2-b)*gamma(1+a/2-c)*gamma(1+a-b-c)) | re: 1+a/2-b-c",
    "3 2 | a, b, c | (a+b+1)/2, 2*c | 1 | sqrt(pi)*gamma(c+1/2)*gamma((a+b+1)/2)*gamma(c-(a+b-1)/2)/(gamma((a+1)/2)*gamma((b+1)/2)*gamma(c-(a-1)/2)*gamma(c-(b-1)/2)) | re: 1+2*c-a-b",
    "3 2 | a, 1-a, c | e, 1+2*c-e | 1 | pi*2^(1-2*c)*gamma(e)*gamma(1+2*c-e)/(gamma((a+e)/2)*gamma((a+1+2*c-e)/2)*gamma((1-a+e)/2)*gamma((2-a+2*c-e)/2)) | re: c",

    /* 2F1: quadratic transformations, as links to the generators of the
       context (used when a generator 2F1, K or E is at the argument of
       the anchor; one step, not nested): from the parameters (a, b; 2b),
       (a, b; 1+a-b), (a, b; (a+b+1)/2) down to (a/2, a/2+1/2; c),
       (a/2, b/2; c), and back up. The regions of validity: the
       transformations hold off the cut [1, inf), on the unit disk (the
       circle is mapped to the cut by 4z/(1+z)^2), and on Re(z) < 1/2 */
    "2 1 | a, b | 2*b | z | (1-z/2)^(-a)*hyp2f1(a/2, a/2+1/2, b+1/2, z^2/(2-z)^2) | cut: z, anchor: z^2/(2-z)^2",
    "2 1 | a, b | 1+a-b | z | (1+z)^(-a)*hyp2f1(a/2, a/2+1/2, 1+a-b, 4*z/(1+z)^2) | disc: z, anchor: 4*z/(1+z)^2",
    "2 1 | a, b | (a+b+1)/2 | z | hyp2f1(a/2, b/2, (a+b+1)/2, 4*z*(1-z)) | re: 1/2-z, anchor: 4*z*(1-z)",
    "2 1 | a, a+1/2 | c | z | (1+sqrt(z))^(-2*a)*hyp2f1(2*a, c-1/2, 2*c-1, 2*sqrt(z)/(1+sqrt(z))) | cut: z, anchor: 2*sqrt(z)/(1+sqrt(z))",
    "2 1 | a, a+1/2 | c | z | ((1+sqrt(1-z))/2)^(-2*a)*hyp2f1(2*a, 1+2*a-c, c, (1-sqrt(1-z))/(1+sqrt(1-z))) | cut: z, anchor: (1-sqrt(1-z))/(1+sqrt(1-z))",
    "2 1 | a, b | a+b+1/2 | z | hyp2f1(2*a, 2*b, a+b+1/2, (1-sqrt(1-z))/2) | cut: z, anchor: (1-sqrt(1-z))/2",

    NULL
};

slong
_gr_tower_hypgeom_table_length(void)
{
    slong n = 0;
    while (_hyp_table[n] != NULL)
        n++;
    return n;
}

const char *
_gr_tower_hypgeom_table_entry(slong i)
{
    return _hyp_table[i];
}

/* a parsed entry */
typedef struct
{
    slong p, q;
    hyp_affine_struct * par;    /* p + q affine forms */
    int z_free;
    fmpq_t z0;
    const char * value;
    const char * cond;          /* NULL if none */
    const char * anchor;        /* the expression of an "anchor:" condition, NULL if none */
    slong anchor_len;
    slong nsym;                 /* symbols occurring */
    int used[HYP_NSYM - 1];
}
hyp_entry_struct;

static void
_hyp_entry_clear(hyp_entry_struct * R)
{
    slong i;
    for (i = 0; i < R->p + R->q; i++)
        _hyp_affine_clear(R->par + i);
    flint_free(R->par);
    fmpq_clear(R->z0);
}

/* parses an entry (the value and conditions are kept as pointers into s) */
static int
_hyp_entry_parse(hyp_entry_struct * R, const char * s)
{
    slong i, k;
    const char * t = s;

    R->p = strtol(t, (char **) &t, 10);
    R->q = strtol(t, (char **) &t, 10);
    /* (the table is internal: a malformed entry is a bug) */
    if (R->p < 0 || R->q < 0 || R->p + R->q > 16)
        flint_throw(FLINT_ERROR, "(%s): malformed entry %s\n", __func__, s);
    R->par = flint_malloc(sizeof(hyp_affine_struct) * FLINT_MAX(R->p + R->q, 1));
    for (i = 0; i < R->p + R->q; i++)
        _hyp_affine_init(R->par + i);
    fmpq_init(R->z0);
    R->cond = NULL;

    for (k = 0; k < 2; k++)
    {
        slong n = (k == 0) ? R->p : R->q, off = (k == 0) ? 0 : R->p;
        _hyp_aff_space(&t);
        if (*t != '|')
            return 0;
        t++;
        for (i = 0; i < n; i++)
        {
            if (!_hyp_aff_sum(R->par + off + i, &t))
                return 0;
            _hyp_aff_space(&t);
            if (i < n - 1)
            {
                if (*t != ',')
                    return 0;
                t++;
            }
        }
    }

    _hyp_aff_space(&t);
    if (*t != '|')
        return 0;
    t++;
    _hyp_aff_space(&t);
    if (*t == 'z')
    {
        R->z_free = 1;
        t++;
    }
    else
    {
        hyp_affine_struct Z;
        _hyp_affine_init(&Z);
        R->z_free = 0;
        if (!_hyp_aff_sum(&Z, &t) || !_hyp_aff_is_const(&Z))
        {
            _hyp_affine_clear(&Z);
            return 0;
        }
        fmpq_set(R->z0, Z.c);
        _hyp_affine_clear(&Z);
    }
    _hyp_aff_space(&t);
    if (*t != '|')
        return 0;
    t++;
    R->value = t;

    {
        const char * c = strchr(t, '|');
        if (c != NULL)
            R->cond = c + 1;
    }

    R->anchor = NULL;
    R->anchor_len = 0;
    if (R->cond != NULL)
    {
        const char * c = strstr(R->cond, "anchor:");
        if (c != NULL)
        {
            const char * end;
            c += 7;
            end = strchr(c, ',');
            R->anchor = c;
            R->anchor_len = (end != NULL) ? (end - c) : (slong) strlen(c);
        }
    }

    for (k = 0; k < HYP_NSYM - 1; k++)
    {
        R->used[k] = 0;
        for (i = 0; i < R->p + R->q; i++)
            if (!fmpq_is_zero(R->par[i].c + 1 + k))
                R->used[k] = 1;
    }

    return 1;
}

/* the value of the entry (up to the separator of the conditions) */
static int
_hyp_entry_value(gr_ptr res, const hyp_entry_struct * R, gr_srcptr * vars, int depth, gr_ctx_t ctx)
{
    char * s;
    slong len;
    int status;

    len = (R->cond != NULL) ? (R->cond - 1 - R->value) : (slong) strlen(R->value);
    s = flint_malloc(len + 1);
    memcpy(s, R->value, len);
    s[len] = '\0';
    status = _hyp_eval_expr(res, s, vars, depth, ctx);
    flint_free(s);
    return status;
}

static int _hyp_off_cut(gr_srcptr z, gr_ctx_t ctx);

/* whether the conditions of the entry hold: 1, 0, or -1 (unknown) */
static int
_hyp_entry_conditions(const hyp_entry_struct * R, gr_srcptr * vars, int depth, gr_ctx_t ctx)
{
    const char * c = R->cond;
    int res = 1;

    while (c != NULL && res == 1)
    {
        const char * end = strchr(c, ',');
        char * s;
        slong len;
        int is_re, is_ne, is_disc, is_cut;
        gr_ptr t;

        while (*c == ' ')
            c++;
        is_re = (strncmp(c, "re:", 3) == 0);
        is_ne = (strncmp(c, "ne:", 3) == 0);
        is_disc = (strncmp(c, "disc:", 5) == 0);
        is_cut = (strncmp(c, "cut:", 4) == 0);
        if (strncmp(c, "anchor:", 7) == 0)
        {
            /* (tested when the entry is matched) */
            c = (end != NULL) ? end + 1 : NULL;
            continue;
        }
        if (!is_re && !is_ne && !is_disc && !is_cut)
            return -1;
        c += is_disc ? 5 : is_cut ? 4 : 3;
        len = (end != NULL) ? (end - c) : (slong) strlen(c);
        s = flint_malloc(len + 1);
        memcpy(s, c, len);
        s[len] = '\0';

        GR_TMP_INIT(t, ctx);
        if (_hyp_eval_expr(t, s, vars, depth, ctx) != GR_SUCCESS)
            res = -1;
        else if (is_re)
        {
            int sgn;
            if (_gr_tower_lazy_re_cmp(&sgn, t, NULL, ctx) != GR_SUCCESS)
                res = -1;
            else
                res = (sgn > 0);
        }
        else if (is_disc)
        {
            /* |t| < 1 (numerically: a point on the circle is not decided) */
            acb_t w;
            arb_t a;
            acb_init(w);
            arb_init(a);
            if (gr_tower_lazy_get_acb(w, t, CHECK_PREC, ctx) != GR_SUCCESS)
                res = -1;
            else
            {
                acb_abs(a, w, CHECK_PREC);
                arb_sub_ui(a, a, 1, CHECK_PREC);
                res = arb_is_negative(a) ? 1 : arb_is_nonnegative(a) ? 0 : -1;
            }
            acb_clear(w);
            arb_clear(a);
        }
        else if (is_cut)
            res = _hyp_off_cut(t, ctx);
        else
        {
            truth_t z = gr_is_zero(t, ctx);
            res = (z == T_FALSE) ? 1 : (z == T_TRUE) ? 0 : -1;
        }
        GR_TMP_CLEAR(t, ctx);
        flint_free(s);

        c = (end != NULL) ? end + 1 : NULL;
    }

    return res;
}

/* -------------------------------------------------------------------- */
/* the contiguity module                                                 */
/* -------------------------------------------------------------------- */

static slong
_hyp_order(slong p, slong q)
{
    return FLINT_MAX(p, q + 1);
}

/*
    The companion matrix C (r x r, row-major) of theta acting on
    (F, theta F, ..., theta^(r-1) F) for the parameters P at z, from the
    operator sum_k c_k theta^k = theta prod (theta + b_j - 1) - z prod
    (theta + a_i). GR_DOMAIN if the leading coefficient vanishes (z = 1
    when p = q + 1).
*/
static int
_hyp_companion(gr_ptr C, const hyp_params_t P, gr_srcptr z, gr_ctx_t ctx)
{
    slong p = P->p, q = P->q, r = _hyp_order(p, q), i, j, sz = ctx->sizeof_elem;
    gr_ptr A, B, c, t;
    int status = GR_SUCCESS;

    A = gr_heap_init_vec(r + 2, ctx);   /* x prod (x + b_j - 1): degree q + 1 */
    B = gr_heap_init_vec(r + 2, ctx);   /* prod (x + a_i): degree p */
    c = gr_heap_init_vec(r + 1, ctx);
    GR_TMP_INIT(t, ctx);

    /* A = x */
    status |= gr_one(ENTRY(A, 1, sz), ctx);
    for (j = 0; j < q; j++)
    {
        /* A *= (x + b_j - 1) */
        status |= gr_sub_ui(t, HB(P, j, sz), 1, ctx);
        for (i = j + 2; i >= 1; i--)
        {
            gr_ptr u;
            GR_TMP_INIT(u, ctx);
            status |= gr_mul(u, ENTRY(A, i, sz), t, ctx);
            status |= gr_add(ENTRY(A, i, sz), u, ENTRY(A, i - 1, sz), ctx);
            GR_TMP_CLEAR(u, ctx);
        }
        status |= gr_mul(ENTRY(A, 0, sz), ENTRY(A, 0, sz), t, ctx);
    }

    status |= gr_one(ENTRY(B, 0, sz), ctx);
    for (i = 0; i < p; i++)
    {
        for (j = i + 1; j >= 1; j--)
        {
            gr_ptr u;
            GR_TMP_INIT(u, ctx);
            status |= gr_mul(u, ENTRY(B, j, sz), HA(P, i, sz), ctx);
            status |= gr_add(ENTRY(B, j, sz), u, ENTRY(B, j - 1, sz), ctx);
            GR_TMP_CLEAR(u, ctx);
        }
        status |= gr_mul(ENTRY(B, 0, sz), ENTRY(B, 0, sz), HA(P, i, sz), ctx);
    }

    for (i = 0; i <= r; i++)
    {
        status |= gr_mul(t, ENTRY(B, i, sz), z, ctx);
        status |= gr_sub(ENTRY(c, i, sz), ENTRY(A, i, sz), t, ctx);
    }

    if (status == GR_SUCCESS)
    {
        truth_t zero = gr_is_zero(ENTRY(c, r, sz), ctx);
        if (zero == T_TRUE)
            status = GR_DOMAIN;
        else if (zero == T_UNKNOWN)
            status = GR_UNABLE;
    }

    if (status == GR_SUCCESS)
    {
        for (i = 0; i < r * r; i++)
            status |= gr_zero(ENTRY(C, i, sz), ctx);
        for (i = 0; i + 1 < r; i++)
            status |= gr_one(ENTRY(C, i * r + i + 1, sz), ctx);
        for (j = 0; j < r; j++)
        {
            status |= gr_div(t, ENTRY(c, j, sz), ENTRY(c, r, sz), ctx);
            status |= gr_neg(ENTRY(C, (r - 1) * r + j, sz), t, ctx);
        }
    }

    gr_heap_clear_vec(A, r + 2, ctx);
    gr_heap_clear_vec(B, r + 2, ctx);
    gr_heap_clear_vec(c, r + 1, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* M = I + C / d (r x r) */
static int
_hyp_step_matrix(gr_mat_t M, const hyp_params_t P, gr_srcptr z, gr_srcptr d, gr_ctx_t ctx)
{
    slong r = gr_mat_nrows(M, ctx), i, j, sz = ctx->sizeof_elem;
    gr_ptr C;
    int status;

    C = gr_heap_init_vec(r * r, ctx);
    status = _hyp_companion(C, P, z, ctx);
    for (i = 0; i < r && status == GR_SUCCESS; i++)
        for (j = 0; j < r && status == GR_SUCCESS; j++)
        {
            status = gr_div(GR_MAT_ENTRY(M, i, j, sz), ENTRY(C, i * r + j, sz), d, ctx);
            if (i == j)
                status |= gr_add_ui(GR_MAT_ENTRY(M, i, j, sz), GR_MAT_ENTRY(M, i, j, sz), 1, ctx);
        }
    gr_heap_clear_vec(C, r * r, ctx);
    return status;
}

/*
    T = the matrix with v(Q) = T v(P0), for Q = P0 + n (n the integer
    shifts, applied to the parameters of P0 in order). The path raises the
    b first, then moves the a, then lowers the b, which keeps the
    chambers b_j - a_i >= 1 of reducible pairs (see the header).
*/
static int
_hyp_path(gr_mat_t T, const hyp_params_t P0, const slong * n, gr_srcptr z, gr_ctx_t ctx)
{
    slong p = P0->p, q = P0->q, r = _hyp_order(p, q), sz = ctx->sizeof_elem, i, k, phase;
    hyp_params_t P;
    gr_mat_t M, U;
    gr_ptr d;
    int status = GR_SUCCESS;

    _hyp_params_init(P, p, q, ctx);
    status |= _gr_vec_set(P->v, P0->v, p + q, ctx);
    gr_mat_init(M, r, r, ctx);
    gr_mat_init(U, r, r, ctx);
    GR_TMP_INIT(d, ctx);
    status |= gr_mat_one(T, ctx);

    for (phase = 0; phase < 3 && status == GR_SUCCESS; phase++)
    {
        slong lo = (phase == 1) ? 0 : p, hi = (phase == 1) ? p : p + q;

        for (i = lo; i < hi && status == GR_SUCCESS; i++)
        {
            int is_a = (i < p);
            slong m = n[i];

            if (!is_a && ((phase == 0 && m <= 0) || (phase == 2 && m >= 0)))
                continue;

            for (k = 0; k < FLINT_ABS(m) && status == GR_SUCCESS; k++)
            {
                gr_ptr x = ENTRY(P->v, i, sz);

                if (is_a && m > 0)
                {
                    /* raise a_i: v' = (I + C/a_i) v */
                    status |= _hyp_step_matrix(M, P, z, x, ctx);
                    status |= gr_add_ui(x, x, 1, ctx);
                    if (status == GR_SUCCESS)
                    {
                        status = gr_mat_mul(U, M, T, ctx);
                        gr_mat_swap(U, T, ctx);
                    }
                }
                else if (is_a)
                {
                    /* lower a_i: v = (I + C(P')/(a_i - 1)) v' */
                    status |= gr_sub_ui(x, x, 1, ctx);
                    status |= _hyp_step_matrix(M, P, z, x, ctx);
                    if (status == GR_SUCCESS)
                    {
                        status = gr_mat_nonsingular_solve(U, M, T, ctx);
                        gr_mat_swap(U, T, ctx);
                    }
                }
                else if (m < 0)
                {
                    /* lower b_j: v' = (I + C/(b_j - 1)) v */
                    status |= gr_sub_ui(d, x, 1, ctx);
                    status |= _hyp_step_matrix(M, P, z, d, ctx);
                    status |= gr_set(x, d, ctx);
                    if (status == GR_SUCCESS)
                    {
                        status = gr_mat_mul(U, M, T, ctx);
                        gr_mat_swap(U, T, ctx);
                    }
                }
                else
                {
                    /* raise b_j: v = (I + C(P')/b_j) v' */
                    status |= gr_set(d, x, ctx);
                    status |= gr_add_ui(x, x, 1, ctx);
                    status |= _hyp_step_matrix(M, P, z, d, ctx);
                    if (status == GR_SUCCESS)
                    {
                        status = gr_mat_nonsingular_solve(U, M, T, ctx);
                        gr_mat_swap(U, T, ctx);
                    }
                }
            }
        }
    }

    /* (a singular step: not expected within a chamber) */
    if (status == GR_DOMAIN)
        status = GR_UNABLE;

    _hyp_params_clear(P, ctx);
    gr_mat_clear(M, ctx);
    gr_mat_clear(U, ctx);
    GR_TMP_CLEAR(d, ctx);
    return status;
}

/*
    S (r x r) with theta^k F(P0) = sum_j S[k][j] G_j, G_j = F(P0 + j e),
    e = (1, ..., 1): theta^k = sum_j S2(k, j) z^j D^j and D^j F =
    (a)_j / (b)_j F(P0 + j e).
*/
static int
_hyp_basis_matrix(gr_mat_t S, const hyp_params_t P0, gr_srcptr z, gr_ctx_t ctx)
{
    slong p = P0->p, q = P0->q, r = _hyp_order(p, q), sz = ctx->sizeof_elem, k, j, i;
    gr_ptr f, t;
    fmpz_t s2;
    int status = GR_SUCCESS;

    f = gr_heap_init_vec(r, ctx);    /* f[j] = z^j (a)_j / (b)_j */
    GR_TMP_INIT(t, ctx);
    fmpz_init(s2);

    status |= gr_one(f, ctx);
    for (j = 1; j < r; j++)
    {
        status |= gr_mul(ENTRY(f, j, sz), ENTRY(f, j - 1, sz), z, ctx);
        for (i = 0; i < p; i++)
        {
            status |= gr_add_si(t, HA(P0, i, sz), j - 1, ctx);
            status |= gr_mul(ENTRY(f, j, sz), ENTRY(f, j, sz), t, ctx);
        }
        for (i = 0; i < q; i++)
        {
            status |= gr_add_si(t, HB(P0, i, sz), j - 1, ctx);
            status |= gr_div(ENTRY(f, j, sz), ENTRY(f, j, sz), t, ctx);
        }
    }

    for (k = 0; k < r; k++)
        for (j = 0; j < r; j++)
        {
            if (j > k)
                status |= gr_zero(GR_MAT_ENTRY(S, k, j, sz), ctx);
            else
            {
                arith_stirling_number_2(s2, k, j);
                status |= gr_mul_fmpz(GR_MAT_ENTRY(S, k, j, sz), ENTRY(f, j, sz), s2, ctx);
            }
        }

    gr_heap_clear_vec(f, r, ctx);
    GR_TMP_CLEAR(t, ctx);
    fmpz_clear(s2);
    return status;
}

/* w (row vector of length r) := w * M */
static int
_hyp_row_mul(gr_ptr w, const gr_mat_t M, gr_ctx_t ctx)
{
    slong r = gr_mat_nrows(M, ctx), i, j, sz = ctx->sizeof_elem;
    gr_ptr u, t;
    int status = GR_SUCCESS;

    u = gr_heap_init_vec(r, ctx);
    GR_TMP_INIT(t, ctx);
    for (j = 0; j < r; j++)
        for (i = 0; i < r; i++)
        {
            status |= gr_mul(t, ENTRY(w, i, sz), GR_MAT_ENTRY(M, i, j, sz), ctx);
            status |= gr_add(ENTRY(u, j, sz), ENTRY(u, j, sz), t, ctx);
        }
    status |= _gr_vec_set(w, u, r, ctx);
    gr_heap_clear_vec(u, r, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* transformations of one term                                           */
/* -------------------------------------------------------------------- */

/*
    The transformations F(P; z) = f(P, z) F(g(P); g(z)), as involutions:
    1F1: Kummer (b - a; b; -z) with f = exp(z); 2F1: Euler
    (c - a, c - b; c; z) with f = (1 - z)^(c - a - b), Pfaff
    (a, c - b; c; z/(z - 1)) with f = (1 - z)^(-a) and the same with a, b
    exchanged (valid off the cut z >= 1).
*/
#define HYP_T_ID 0
#define HYP_T_KUMMER 1
#define HYP_T_EULER 2
#define HYP_T_PFAFF_A 3
#define HYP_T_PFAFF_B 4

/* applies the transformation: R, w = g(P), g(z), and f (unless f is NULL) */
static int
_hyp_transform(hyp_params_t R, gr_ptr w, gr_ptr f, int g, const hyp_params_t P, gr_srcptr z, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= _hyp_params_set(R, P, ctx);

    if (g == HYP_T_ID)
    {
        status |= gr_set(w, z, ctx);
        if (f != NULL)
            status |= gr_one(f, ctx);
    }
    else if (g == HYP_T_KUMMER)
    {
        status |= gr_sub(HA(R, 0, sz), HB(P, 0, sz), HA(P, 0, sz), ctx);
        status |= gr_neg(w, z, ctx);
        if (f != NULL)
            status |= gr_exp(f, z, ctx);
    }
    else
    {
        gr_srcptr a = HA(P, 0, sz), b = HA(P, 1, sz), c = HB(P, 0, sz);
        gr_ptr one_minus_z;
        GR_TMP_INIT(one_minus_z, ctx);
        status |= gr_sub_ui(one_minus_z, z, 1, ctx);
        status |= gr_neg(one_minus_z, one_minus_z, ctx);

        if (g == HYP_T_EULER)
        {
            status |= gr_sub(HA(R, 0, sz), c, a, ctx);
            status |= gr_sub(HA(R, 1, sz), c, b, ctx);
            status |= gr_set(w, z, ctx);
            status |= gr_sub(t, c, a, ctx);
            status |= gr_sub(t, t, b, ctx);
        }
        else
        {
            gr_srcptr keep = (g == HYP_T_PFAFF_A) ? a : b;
            gr_srcptr other = (g == HYP_T_PFAFF_A) ? b : a;
            status |= gr_set(HA(R, 0, sz), keep, ctx);
            status |= gr_sub(HA(R, 1, sz), c, other, ctx);
            status |= gr_sub_ui(t, z, 1, ctx);
            status |= gr_div(w, z, t, ctx);
            status |= gr_neg(t, keep, ctx);
        }

        if (status == GR_SUCCESS && f != NULL)
        {
            truth_t zero = gr_is_zero(t, ctx);
            if (zero == T_TRUE)
                status = gr_one(f, ctx);
            else
                status = _hyp_pow(f, one_minus_z, t, ctx);
        }

        GR_TMP_CLEAR(one_minus_z, ctx);
    }

    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* whether z is off the cut [1, inf) (1, 0, -1 unknown) */
static int
_hyp_off_cut(gr_srcptr z, gr_ctx_t ctx)
{
    acb_t w;
    int res;

    acb_init(w);
    if (gr_tower_lazy_get_acb(w, z, CHECK_PREC, ctx) != GR_SUCCESS)
        res = -1;
    else if (!arb_contains_zero(acb_imagref(w)))
        res = 1;
    else
    {
        arb_t t;
        arb_init(t);
        arb_sub_ui(t, acb_realref(w), 1, CHECK_PREC);
        if (arb_is_negative(t))
            res = 1;
        else
        {
            /* real exactly? (then on the cut if >= 1) */
            truth_t real = gr_tower_lazy_is_real(z, ctx);
            if (real == T_FALSE)
                res = 1;
            else if (real == T_TRUE)
            {
                int sgn;
                gr_ptr u;
                GR_TMP_INIT(u, ctx);
                if (gr_sub_ui(u, z, 1, ctx) == GR_SUCCESS && _gr_tower_lazy_re_cmp(&sgn, u, NULL, ctx) == GR_SUCCESS)
                    res = (sgn < 0);
                else
                    res = -1;
                GR_TMP_CLEAR(u, ctx);
            }
            else
                res = -1;
        }
        arb_clear(t);
    }
    acb_clear(w);
    return res;
}

/* the transformations applicable to (P, z): a list ending with -1 */
static void
_hyp_transforms(int * gs, const hyp_params_t P, gr_srcptr z, gr_ctx_t ctx)
{
    slong n = 0;
    gs[n++] = HYP_T_ID;
    if (P->p == 1 && P->q == 1)
        gs[n++] = HYP_T_KUMMER;
    else if (P->p == 2 && P->q == 1 && _hyp_off_cut(z, ctx) == 1)
    {
        gs[n++] = HYP_T_EULER;
        gs[n++] = HYP_T_PFAFF_A;
        gs[n++] = HYP_T_PFAFF_B;
    }
    gs[n] = -1;
}

/* -------------------------------------------------------------------- */
/* matching the table                                                    */
/* -------------------------------------------------------------------- */

/* next permutation of perm[0..n-1] (lexicographic); 0 when done */
static int
_hyp_next_perm(slong * perm, slong n)
{
    slong i, j, t;
    if (n <= 1)
        return 0;
    i = n - 2;
    while (i >= 0 && perm[i] >= perm[i + 1])
        i--;
    if (i < 0)
        return 0;
    j = n - 1;
    while (perm[j] <= perm[i])
        j--;
    t = perm[i]; perm[i] = perm[j]; perm[j] = t;
    for (i = i + 1, j = n - 1; i < j; i++, j--)
    {
        t = perm[i]; perm[i] = perm[j]; perm[j] = t;
    }
    return 1;
}

/*
    Matches the entry R against the coset of the parameters Pm at the
    argument w: for every assignment of the parameters to the slots of the
    pattern, the symbols are solved for (from slots where a symbol occurs
    alone with coefficient +-1); the other slots must then differ from
    the parameters by rationals which shifts of the symbols make integral.
    Every point of the coset obtained by shifting the symbols by at most
    HYP_BOX (0 with exact_only) is passed to found(), with its offsets from
    Pm (in the order of Pm) and the values of the symbols there (z = w).
*/
typedef int (* hyp_found_func)(const slong * off, gr_srcptr * vars, void * data);

static int
_hyp_match_entry(const hyp_entry_struct * R, const hyp_params_t Pm, gr_srcptr w,
    hyp_found_func found, void * data, int exact_only, gr_ctx_t ctx)
{
    slong p = Pm->p, q = Pm->q, sz = ctx->sizeof_elem, n = p + q, i, k;
    slong * pa, * pb, * slot_of;   /* slot_of[slot] = index of the parameter in Pm */
    gr_srcptr vars[HYP_NVARS];
    gr_ptr symvals, t;
    int status = GR_SUCCESS;

    if (R->p != p || R->q != q)
        return GR_SUCCESS;

    /* the argument */
    if (!R->z_free)
    {
        gr_ptr c;
        truth_t eq;
        GR_TMP_INIT(c, ctx);
        status = gr_set_fmpq(c, R->z0, ctx);
        eq = (status == GR_SUCCESS) ? gr_equal(c, w, ctx) : T_UNKNOWN;
        GR_TMP_CLEAR(c, ctx);
        if (eq != T_TRUE)
            return GR_SUCCESS;
    }

    pa = flint_malloc(sizeof(slong) * FLINT_MAX(p, 1));
    pb = flint_malloc(sizeof(slong) * FLINT_MAX(q, 1));
    slot_of = flint_malloc(sizeof(slong) * FLINT_MAX(n, 1));
    symvals = gr_heap_init_vec(HYP_NVARS, ctx);
    GR_TMP_INIT(t, ctx);

    for (i = 0; i < p; i++)
        pa[i] = i;

    do
    {
        for (i = 0; i < q; i++)
            pb[i] = i;
        do
        {
            int ok = 1;
            slong sym;
            int bound[HYP_NVARS];
            fmpq * dres;   /* the residues: value of the slot - parameter */

            for (i = 0; i < p; i++)
                slot_of[i] = pa[i];
            for (i = 0; i < q; i++)
                slot_of[p + i] = p + pb[i];

            for (sym = 0; sym < HYP_NVARS; sym++)
            {
                bound[sym] = 0;
                vars[sym] = NULL;
            }

            /* bind the symbols from slots where one occurs alone, with
               coefficient +-1: sym = +-(param - rest) */
            {
                int progress = 1;
                while (progress && ok)
                {
                    progress = 0;
                    for (k = 0; k < n && ok; k++)
                    {
                        const hyp_affine_struct * A = R->par + k;
                        slong unb = -1, cnt = 0;
                        hyp_affine_struct B;

                        for (sym = 0; sym < HYP_NSYM - 1; sym++)
                            if (!fmpq_is_zero(A->c + 1 + sym) && !bound[sym])
                            {
                                cnt++;
                                unb = sym;
                            }
                        if (cnt != 1)
                            continue;
                        if (!(fmpz_is_pm1(fmpq_numref(A->c + 1 + unb)) && fmpz_is_one(fmpq_denref(A->c + 1 + unb))))
                            continue;

                        _hyp_affine_init(&B);
                        for (i = 0; i < HYP_NSYM; i++)
                            fmpq_set(B.c + i, A->c + i);
                        fmpq_zero(B.c + 1 + unb);
                        status = _hyp_aff_eval(t, &B, vars, ctx);
                        _hyp_affine_clear(&B);

                        if (status == GR_SUCCESS)
                            status = gr_sub(ENTRY(symvals, unb, sz), ENTRY(Pm->v, slot_of[k], sz), t, ctx);
                        if (status == GR_SUCCESS && !fmpz_is_one(fmpq_numref(A->c + 1 + unb)))
                            status = gr_neg(ENTRY(symvals, unb, sz), ENTRY(symvals, unb, sz), ctx);
                        if (status != GR_SUCCESS)
                        {
                            ok = 0;
                            break;
                        }
                        vars[unb] = ENTRY(symvals, unb, sz);
                        bound[unb] = 1;
                        progress = 1;
                    }
                }
            }

            for (sym = 0; sym < HYP_NSYM - 1 && ok; sym++)
                if (R->used[sym] && !bound[sym])
                    ok = 0;

            dres = _fmpq_vec_init(n);
            for (k = 0; k < n && ok; k++)
            {
                status = _hyp_aff_eval(t, R->par + k, vars, ctx);
                status |= gr_sub(t, t, ENTRY(Pm->v, slot_of[k], sz), ctx);
                if (status != GR_SUCCESS || _hyp_rational(dres + k, t, ctx) != 1)
                    ok = 0;
            }
            status = GR_SUCCESS;

            if (ok)
            {
                /* the shifts m of the bound symbols in the box for which
                   all the residues become integers */
                slong nb = 0, symlist[HYP_NVARS], m[HYP_NVARS], total, idx, base;
                slong * off = flint_malloc(sizeof(slong) * n);
                gr_ptr shifted = gr_heap_init_vec(HYP_NVARS, ctx);
                gr_srcptr svars[HYP_NVARS];

                for (sym = 0; sym < HYP_NSYM - 1; sym++)
                    if (bound[sym])
                        symlist[nb++] = sym;

                base = exact_only ? 1 : 2 * HYP_BOX + 1;
                total = 1;
                for (i = 0; i < nb; i++)
                    total *= base;

                for (idx = 0; idx < total && status == GR_SUCCESS; idx++)
                {
                    slong rem = idx;
                    int good = 1;
                    fmpq_t s, u;

                    for (i = 0; i < nb; i++)
                    {
                        m[i] = exact_only ? 0 : (rem % base) - HYP_BOX;
                        rem /= base;
                    }

                    fmpq_init(s);
                    fmpq_init(u);
                    for (k = 0; k < n && good; k++)
                    {
                        fmpq_set(s, dres + k);
                        for (i = 0; i < nb; i++)
                        {
                            fmpq_mul_si(u, R->par[k].c + 1 + symlist[i], m[i]);
                            fmpq_add(s, s, u);
                        }
                        if (!fmpz_is_one(fmpq_denref(s)) || !fmpz_fits_si(fmpq_numref(s)) ||
                            FLINT_ABS(fmpz_get_si(fmpq_numref(s))) > HYP_SHIFT_LIMIT)
                            good = 0;
                        else
                            off[slot_of[k]] = fmpz_get_si(fmpq_numref(s));
                    }
                    fmpq_clear(s);
                    fmpq_clear(u);

                    if (!good)
                        continue;

                    for (sym = 0; sym < HYP_NVARS; sym++)
                        svars[sym] = NULL;
                    for (i = 0; i < nb; i++)
                    {
                        status = gr_add_si(ENTRY(shifted, symlist[i], sz), vars[symlist[i]], m[i], ctx);
                        svars[symlist[i]] = ENTRY(shifted, symlist[i], sz);
                    }
                    svars['z' - 'a'] = w;

                    if (status == GR_SUCCESS)
                        status = found(off, svars, data);
                }

                gr_heap_clear_vec(shifted, HYP_NVARS, ctx);
                flint_free(off);
            }

            _fmpq_vec_clear(dres, n);
        }
        while (status == GR_SUCCESS && _hyp_next_perm(pb, q));
    }
    while (status == GR_SUCCESS && _hyp_next_perm(pa, p));

    flint_free(pa);
    flint_free(pb);
    flint_free(slot_of);
    gr_heap_clear_vec(symvals, HYP_NVARS, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* layer 3: the coset                                                    */
/* -------------------------------------------------------------------- */

/*
    The cell of P: P0 = P - n with 0 < Re(x) <= 1 for every parameter,
    except that for a reducible pair (b_j - a_i a positive integer, with
    a_i in its cell) b_j is placed at a_i + 1. Sets red[j] = i for such a
    b_j (-1 otherwise).
*/
static int
_hyp_cell(hyp_params_t P0, slong * n, slong * red, const hyp_params_t P, gr_ctx_t ctx)
{
    slong p = P->p, q = P->q, sz = ctx->sizeof_elem, i, j;
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= _hyp_params_set(P0, P, ctx);

    for (i = 0; i < p && status == GR_SUCCESS; i++)
    {
        status = _hyp_cell_shift(n + i, HA(P, i, sz), ctx);
        if (status == GR_SUCCESS)
            status = gr_sub_si(HA(P0, i, sz), HA(P, i, sz), n[i], ctx);
    }

    for (j = 0; j < q && status == GR_SUCCESS; j++)
    {
        red[j] = -1;
        for (i = 0; i < p && status == GR_SUCCESS; i++)
        {
            slong d;
            int r;
            status = gr_sub(t, HB(P, j, sz), HA(P, i, sz), ctx);
            r = (status == GR_SUCCESS) ? _hyp_integer(&d, t, ctx) : 0;
            if (r == -1)
                status = GR_UNABLE;
            else if (r == 1 && d >= 1)
            {
                /* b_j = a_i + d: the cell has b_j = a_i0 + 1 */
                red[j] = i;
                n[p + j] = n[i] + d - 1;
                status = gr_add_ui(HB(P0, j, sz), HA(P0, i, sz), 1, ctx);
                break;
            }
        }
        if (status == GR_SUCCESS && red[j] < 0)
        {
            status = _hyp_cell_shift(n + p + j, HB(P, j, sz), ctx);
            if (status == GR_SUCCESS)
                status = gr_sub_si(HB(P0, j, sz), HB(P, j, sz), n[p + j], ctx);
        }
    }

    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* whether the point P0 + off of the coset is usable: no parameter a
   nonpositive integer, the chambers of the reducible pairs kept */
static int
_hyp_point_ok(const hyp_params_t P0, const slong * off, const slong * red, gr_ctx_t ctx)
{
    slong p = P0->p, q = P0->q, sz = ctx->sizeof_elem, i, j, tot = 0;
    gr_ptr t;
    int ok = 1;

    for (i = 0; i < p + q; i++)
        tot += FLINT_ABS(off[i]);
    if (tot > HYP_SHIFT_LIMIT)
        return 0;

    for (j = 0; j < q && ok; j++)
        if (red[j] >= 0 && off[p + j] - off[red[j]] < 0)
            ok = 0;

    GR_TMP_INIT(t, ctx);
    for (i = 0; i < p + q && ok; i++)
    {
        slong m;
        int r;
        if (gr_add_si(t, ENTRY(P0->v, i, sz), off[i], ctx) != GR_SUCCESS)
            ok = 0;
        else
        {
            r = _hyp_integer(&m, t, ctx);
            if (r == -1 || (r == 1 && m <= 0))
                ok = 0;
        }
    }
    GR_TMP_CLEAR(t, ctx);
    return ok;
}

/*
    A candidate relation for the basis values of a coset: a reducible pair
    (theta + a_i) F(P0) = a_i F_red, or a point P0 + off of the coset whose
    value is known from a table entry (matched through a transformation g:
    the entry gives F(Q'; w) at the point Q' = g(P0 + off), and F(P0 +
    off; z) = F(Q'; w) / f(Q', w)). The values are only computed for the
    relations which are used.
*/
typedef struct
{
    int type;            /* 0: reducible pair, 1: table point, 2: a positive integer numerator */
    slong i, j;          /* the pair */
    slong * off;
    slong dist;
    int zconst;          /* the entry has a constant argument */
    slong entry;
    gr_ptr vars;         /* HYP_NVARS values of the symbols */
    int bound[HYP_NVARS];
    int g;
    gr_ptr Qp;           /* the matched point (p + q values) */
    gr_ptr w;
    int anchored;        /* from an entry with an anchor (a linking step) */
}
hyp_cand_struct;

typedef struct
{
    hyp_cand_struct * c;
    slong len, alloc;
    slong p, q;
    const hyp_params_struct * P0;   /* the cell */
    const slong * red;
    int g;                          /* the transformation being matched */
    const hyp_params_struct * Pt;   /* g(P0) */
    gr_srcptr w;                    /* g(z) */
    slong entry;
    int zconst;
    int anchored;
    slong limit;
    gr_ctx_struct * ctx;
}
hyp_cands_struct;

static void
_hyp_cands_clear(hyp_cands_struct * C, gr_ctx_t ctx)
{
    slong i, n = C->p + C->q;
    for (i = 0; i < C->len; i++)
    {
        flint_free(C->c[i].off);
        if (C->c[i].type == 1)
        {
            gr_heap_clear_vec(C->c[i].vars, HYP_NVARS, ctx);
            gr_heap_clear_vec(C->c[i].Qp, FLINT_MAX(n, 1), ctx);
            gr_heap_clear(C->c[i].w, ctx);
        }
    }
    flint_free(C->c);
}

static hyp_cand_struct *
_hyp_cands_push(hyp_cands_struct * C)
{
    hyp_cand_struct * c;
    if (C->len == C->alloc)
    {
        C->alloc = FLINT_MAX(8, 2 * C->alloc);
        C->c = flint_realloc(C->c, sizeof(hyp_cand_struct) * C->alloc);
    }
    c = C->c + C->len;
    C->len++;
    c->off = flint_calloc(FLINT_MAX(C->p + C->q, 1), sizeof(slong));
    c->anchored = 0;
    return c;
}

/* found(): the point Q' = g(P0) + off' of the transformed coset; its
   preimage Q = g(Q') must be P0 + off with integer off */
static int
_hyp_cands_found(const slong * offt, gr_srcptr * vars, void * data)
{
    hyp_cands_struct * C = data;
    gr_ctx_struct * ctx = C->ctx;
    const hyp_params_struct * P0 = C->P0;
    slong p = C->p, q = C->q, n = p + q, sz = ctx->sizeof_elem, i, k;
    hyp_params_t Qp, Q;
    gr_ptr ww, f, dif;
    slong * oq;
    int ok = 1, status = GR_SUCCESS;

    if (C->len >= C->limit)
        return GR_SUCCESS;

    _hyp_params_init(Qp, p, q, ctx);
    _hyp_params_init(Q, p, q, ctx);
    GR_TMP_INIT3(ww, f, dif, ctx);
    oq = flint_malloc(sizeof(slong) * n);

    for (i = 0; i < n && status == GR_SUCCESS; i++)
        status = gr_add_si(ENTRY(Qp->v, i, sz), ENTRY(C->Pt->v, i, sz), offt[i], ctx);

    if (status == GR_SUCCESS)
        status = _hyp_transform(Q, ww, NULL, C->g, Qp, C->w, ctx);

    /* the offsets of Q from P0 (for 2F1 the transformation may exchange
       a and b) */
    ok = 0;
    if (status == GR_SUCCESS)
    {
        int perm;
        for (perm = 0; perm < ((p == 2) ? 2 : 1) && !ok; perm++)
        {
            ok = 1;
            for (i = 0; i < n && ok; i++)
            {
                slong src = (perm == 1 && i < 2) ? 1 - i : i;
                slong m;
                if (gr_sub(dif, ENTRY(Q->v, src, sz), ENTRY(P0->v, i, sz), ctx) != GR_SUCCESS ||
                    _hyp_integer(&m, dif, ctx) != 1)
                    ok = 0;
                else
                    oq[i] = m;
            }
        }
    }

    if (ok && _hyp_point_ok(P0, oq, C->red, ctx))
    {
        /* (the same point from the same entry: once) */
        for (k = 0; k < C->len && ok; k++)
        {
            if (C->c[k].type == 1 && C->c[k].entry == C->entry)
            {
                for (i = 0; i < n; i++)
                    if (C->c[k].off[i] != oq[i])
                        break;
                if (i == n)
                    ok = 0;
            }
        }
    }
    else
        ok = 0;

    if (ok)
    {
        hyp_cand_struct * c = _hyp_cands_push(C);
        c->type = 1;
        c->dist = 0;
        for (i = 0; i < n; i++)
        {
            c->off[i] = oq[i];
            c->dist += FLINT_ABS(oq[i]);
        }
        c->zconst = C->zconst;
        c->entry = C->entry;
        c->anchored = C->anchored;
        c->g = C->g;
        c->vars = gr_heap_init_vec(HYP_NVARS, ctx);
        for (i = 0; i < HYP_NVARS; i++)
        {
            c->bound[i] = (vars[i] != NULL);
            if (vars[i] != NULL)
                status |= gr_set(ENTRY(c->vars, i, sz), vars[i], ctx);
        }
        c->Qp = gr_heap_init_vec(FLINT_MAX(n, 1), ctx);
        status |= _gr_vec_set(c->Qp, Qp->v, n, ctx);
        c->w = gr_heap_init(ctx);
        status |= gr_set(c->w, C->w, ctx);
    }

    flint_free(oq);
    GR_TMP_CLEAR3(ww, f, dif, ctx);
    _hyp_params_clear(Qp, ctx);
    _hyp_params_clear(Q, ctx);
    return status;
}

/* the value of a table candidate: F(P0 + off; z) */
static int
_hyp_cand_value(gr_ptr res, const hyp_cand_struct * c, slong p, slong q, int depth, gr_ctx_t ctx)
{
    hyp_entry_struct R;
    gr_srcptr vars[HYP_NVARS];
    slong i, sz = ctx->sizeof_elem;
    int status;

    for (i = 0; i < HYP_NVARS; i++)
        vars[i] = c->bound[i] ? ENTRY(c->vars, i, sz) : NULL;

    if (!_hyp_entry_parse(&R, _hyp_table[c->entry]))
    {
        _hyp_entry_clear(&R);
        return GR_UNABLE;
    }

    if (c->anchored && _gr_tower_lazy_hyp_anchored(ctx, 0) > 0)
        status = GR_UNABLE;     /* (one linking step at a time) */
    else if (_hyp_entry_conditions(&R, vars, depth, ctx) != 1)
        status = GR_UNABLE;
    else
    {
        if (c->anchored)
            _gr_tower_lazy_hyp_anchored(ctx, 1);
        status = _hyp_entry_value(res, &R, vars, depth + 1, ctx);
        if (c->anchored)
            _gr_tower_lazy_hyp_anchored(ctx, -1);
    }
    _hyp_entry_clear(&R);

    if (status == GR_SUCCESS && c->g != HYP_T_ID)
    {
        /* F(Q) = F(Q'; w) / f(Q', w) */
        hyp_params_t Qp, Q;
        gr_ptr ww, f;
        _hyp_params_init(Qp, p, q, ctx);
        _hyp_params_init(Q, p, q, ctx);
        GR_TMP_INIT2(ww, f, ctx);
        status = _gr_vec_set(Qp->v, c->Qp, p + q, ctx);
        status |= _hyp_transform(Q, ww, f, c->g, Qp, c->w, ctx);
        if (status == GR_SUCCESS)
            status = gr_div(res, res, f, ctx);
        GR_TMP_CLEAR2(ww, f, ctx);
        _hyp_params_clear(Qp, ctx);
        _hyp_params_clear(Q, ctx);
    }

    return status;
}

/* for an entry with an anchor (a linking step, to the values of another
   hypergeometric function at another argument): whether a generator of
   the context is at that argument, numerically, and no linking step is
   in progress */
static int
_hyp_anchor_ok(const hyp_entry_struct * R, gr_srcptr w, gr_ctx_t ctx)
{
    gr_srcptr vars[HYP_NVARS];
    acb_t t;
    slong i;
    int ok;

    if (R->anchor == NULL)
        return 1;
    if (_gr_tower_lazy_hyp_anchored(ctx, 0) > 0)
        return 0;

    for (i = 0; i < HYP_NVARS; i++)
        vars[i] = NULL;
    vars['z' - 'a'] = w;
    acb_init(t);
    ok = (_hyp_eval_expr_acb(t, R->anchor, R->anchor_len, vars, ctx) == GR_SUCCESS) &&
        _gr_tower_lazy_hyp_anchor_present(t, 0, ctx);
    acb_clear(t);
    return ok;
}

/* collects the table candidates for the coset of P0 at z */
static int
_hyp_collect(hyp_cands_struct * C, const hyp_params_t P0, const slong * red, gr_srcptr z, int exact_only, gr_ctx_t ctx)
{
    int gs[8];
    slong ti, e, pass;
    int status = GR_SUCCESS;

    if (exact_only)
        gs[0] = HYP_T_ID, gs[1] = -1;
    else
        _hyp_transforms(gs, P0, z, ctx);

    /* the entries with a constant argument (summation theorems) first */
    for (pass = 0; pass < 2; pass++)
    {
        for (ti = 0; gs[ti] >= 0 && status == GR_SUCCESS; ti++)
        {
            hyp_params_t Pt;
            gr_ptr w;

            _hyp_params_init(Pt, P0->p, P0->q, ctx);
            GR_TMP_INIT(w, ctx);

            status = _hyp_transform(Pt, w, NULL, gs[ti], P0, z, ctx);

            C->P0 = P0;
            C->red = red;
            C->g = gs[ti];
            C->Pt = Pt;
            C->w = w;

            for (e = 0; _hyp_table[e] != NULL && status == GR_SUCCESS; e++)
            {
                hyp_entry_struct R;
                if (_hyp_entry_parse(&R, _hyp_table[e]) && R.p == P0->p && R.q == P0->q &&
                    (pass == 0) == (R.z_free == 0) && _hyp_anchor_ok(&R, w, ctx))
                {
                    C->entry = e;
                    C->zconst = !R.z_free;
                    C->anchored = (R.anchor != NULL);
                    status = _hyp_match_entry(&R, Pt, w, _hyp_cands_found, C, exact_only, ctx);
                }
                _hyp_entry_clear(&R);
            }

            GR_TMP_CLEAR(w, ctx);
            _hyp_params_clear(Pt, ctx);
        }
    }

    /* (a failure to match is not an error) */
    if (status == GR_UNABLE)
        status = GR_SUCCESS;
    return status;
}

/*
    The relations found so far, in reduced row echelon form: rows[k] has
    a 1 in the column pcol[k] and zeros in the other pivot columns;
    rows[k] . G = vals[k].
*/
typedef struct
{
    slong r, len;
    gr_ptr * rows;
    gr_ptr * vals;
    slong * pcol;
}
hyp_rref_struct;

/* reduces row (and val, if not NULL) by the relations */
static int
_hyp_rref_reduce(gr_ptr row, gr_ptr val, const hyp_rref_struct * E, gr_ctx_t ctx)
{
    slong k, j, sz = ctx->sizeof_elem, r = E->r;
    gr_ptr t, u;
    int status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);
    for (k = 0; k < E->len && status == GR_SUCCESS; k++)
    {
        slong c = E->pcol[k];
        if (gr_is_zero(ENTRY(row, c, sz), ctx) == T_TRUE)
            continue;
        status |= gr_set(t, ENTRY(row, c, sz), ctx);
        for (j = 0; j < r; j++)
        {
            status |= gr_mul(u, ENTRY(E->rows[k], j, sz), t, ctx);
            status |= gr_sub(ENTRY(row, j, sz), ENTRY(row, j, sz), u, ctx);
        }
        if (val != NULL)
        {
            status |= gr_mul(u, E->vals[k], t, ctx);
            status |= gr_sub(val, val, u, ctx);
        }
    }
    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

/* the pivot column of a reduced row: the nonzero entry of the highest
   index; -1 if the row is zero, -2 if unknown */
static slong
_hyp_rref_pivot(gr_srcptr row, slong r, gr_ctx_t ctx)
{
    slong c, sz = ctx->sizeof_elem;
    for (c = r - 1; c >= 0; c--)
    {
        truth_t z = gr_is_zero(ENTRY(row, c, sz), ctx);
        if (z == T_FALSE)
            return c;
        if (z == T_UNKNOWN)
            return -2;
    }
    return -1;
}

/* adds the reduced row with pivot c and value val */
static int
_hyp_rref_add(hyp_rref_struct * E, gr_ptr row, gr_ptr val, slong c, gr_ctx_t ctx)
{
    slong k, j, sz = ctx->sizeof_elem, r = E->r;
    gr_ptr t, u;
    int status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);
    status |= gr_inv(t, ENTRY(row, c, sz), ctx);
    for (j = 0; j < r; j++)
        status |= gr_mul(ENTRY(row, j, sz), ENTRY(row, j, sz), t, ctx);
    status |= gr_mul(val, val, t, ctx);

    /* eliminate the new pivot from the other rows */
    for (k = 0; k < E->len && status == GR_SUCCESS; k++)
    {
        if (gr_is_zero(ENTRY(E->rows[k], c, sz), ctx) == T_TRUE)
            continue;
        status |= gr_set(t, ENTRY(E->rows[k], c, sz), ctx);
        for (j = 0; j < r; j++)
        {
            status |= gr_mul(u, ENTRY(row, j, sz), t, ctx);
            status |= gr_sub(ENTRY(E->rows[k], j, sz), ENTRY(E->rows[k], j, sz), u, ctx);
        }
        status |= gr_mul(u, val, t, ctx);
        status |= gr_sub(E->vals[k], E->vals[k], u, ctx);
    }

    E->rows[E->len] = gr_heap_init_vec(r, ctx);
    E->vals[E->len] = gr_heap_init(ctx);
    status |= _gr_vec_set(E->rows[E->len], row, r, ctx);
    status |= gr_set(E->vals[E->len], val, ctx);
    E->pcol[E->len] = c;
    E->len++;

    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

/* sorts the candidates: reducible pairs, then constant arguments, then by
   distance from the cell */
static int
_hyp_cand_cmp(const void * x, const void * y)
{
    const hyp_cand_struct * a = x, * b = y;
    int ka = (a->type != 1) ? 0 : (a->zconst ? 1 : a->anchored ? 3 : 2);
    int kb = (b->type != 1) ? 0 : (b->zconst ? 1 : b->anchored ? 3 : 2);
    if (ka != kb)
        return ka - kb;
    if (a->dist != b->dist)
        return (a->dist < b->dist) ? -1 : 1;
    return (a->entry < b->entry) ? -1 : (a->entry > b->entry);
}

/*
    For a numerator parameter a_i = 1 (a positive integer, in its cell), the
    operator factors as theta L1 with L1 = prod_j (theta + b_j - 1) -
    z prod_{k != i} (theta + a_k), so that L1 F = L1 F (0) = prod_j (b_j - 1):
    the row of L1 in the basis v(P0) (degree < r), and that value.
*/
static int
_hyp_int_a_row(gr_ptr row, gr_ptr val, const hyp_params_t P0, slong ia, gr_srcptr z, gr_ctx_t ctx)
{
    slong p = P0->p, q = P0->q, r = _hyp_order(p, q), sz = ctx->sizeof_elem, i, j, k, deg;
    gr_ptr A, B, t;
    int status = GR_SUCCESS;

    A = gr_heap_init_vec(r + 1, ctx);
    B = gr_heap_init_vec(r + 1, ctx);
    GR_TMP_INIT(t, ctx);

    status |= gr_one(A, ctx);
    status |= gr_one(val, ctx);
    for (j = 0, deg = 0; j < q; j++, deg++)
    {
        status |= gr_sub_ui(t, HB(P0, j, sz), 1, ctx);
        status |= gr_mul(val, val, t, ctx);
        for (k = deg + 1; k >= 1; k--)
        {
            gr_ptr u;
            GR_TMP_INIT(u, ctx);
            status |= gr_mul(u, ENTRY(A, k, sz), t, ctx);
            status |= gr_add(ENTRY(A, k, sz), u, ENTRY(A, k - 1, sz), ctx);
            GR_TMP_CLEAR(u, ctx);
        }
        status |= gr_mul(A, A, t, ctx);
    }

    status |= gr_one(B, ctx);
    for (i = 0, deg = 0; i < p; i++)
    {
        if (i == ia)
            continue;
        for (k = deg + 1; k >= 1; k--)
        {
            gr_ptr u;
            GR_TMP_INIT(u, ctx);
            status |= gr_mul(u, ENTRY(B, k, sz), HA(P0, i, sz), ctx);
            status |= gr_add(ENTRY(B, k, sz), u, ENTRY(B, k - 1, sz), ctx);
            GR_TMP_CLEAR(u, ctx);
        }
        status |= gr_mul(B, B, HA(P0, i, sz), ctx);
        deg++;
    }

    for (k = 0; k < r; k++)
    {
        status |= gr_mul(t, ENTRY(B, k, sz), z, ctx);
        status |= gr_sub(ENTRY(row, k, sz), ENTRY(A, k, sz), t, ctx);
    }

    gr_heap_clear_vec(A, r + 1, ctx);
    gr_heap_clear_vec(B, r + 1, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/*
    The row (in the basis of S) of the relation of the candidate cc:
    (theta + a_i) F(P0) = a_i F_red (type 0), the integration in a_i
    (type 2: also its value val), or the path to an entry.
*/
static int
_hyp_cand_row(gr_ptr row, gr_ptr val, const hyp_cand_struct * cc, const hyp_params_t P0, const gr_mat_t S, slong r, gr_srcptr z, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, j;
    int status = GR_SUCCESS;

    if (cc->type == 0)
    {
        status |= _gr_vec_zero(row, r, ctx);
        status |= gr_set(GR_ENTRY(row, 0, sz), HA(P0, cc->i, sz), ctx);
        if (r >= 2)
            status |= gr_one(GR_ENTRY(row, 1, sz), ctx);
    }
    else if (cc->type == 2)
    {
        status |= _hyp_int_a_row(row, val, P0, cc->i, z, ctx);
    }
    else
    {
        gr_mat_t Tk;
        gr_mat_init(Tk, r, r, ctx);
        status = _hyp_path(Tk, P0, cc->off, z, ctx);
        for (j = 0; j < r && status == GR_SUCCESS; j++)
            status |= gr_set(GR_ENTRY(row, j, sz), GR_MAT_ENTRY(Tk, 0, j, sz), ctx);
        gr_mat_clear(Tk, ctx);
    }
    if (status == GR_SUCCESS)
        status = _hyp_row_mul(row, S, ctx);
    return status;
}

/*
    F(P; z) through the contiguity module of its coset: P has no
    cancelling pairs, no nonpositive integer parameters and no pairs with
    a_i - b_j a nonnegative integer; z != 0 (and z != 1 if p = q + 1).
*/
static int
_hyp_coset(gr_ptr res, const hyp_params_t P, gr_srcptr z, int depth, gr_ctx_t ctx)
{
    slong p = P->p, q = P->q, r = _hyp_order(p, q), sz = ctx->sizeof_elem, i, j, c;
    hyp_params_t P0;
    slong * n, * red;
    gr_mat_t T, S;
    gr_ptr rho, row, val;
    hyp_cands_struct C;
    hyp_rref_struct E;
    int status = GR_SUCCESS, done = 0;

    if (r > HYP_ORDER_LIMIT)
        return GR_UNABLE;

    _hyp_params_init(P0, p, q, ctx);
    n = flint_calloc(p + q, sizeof(slong));
    red = flint_malloc(sizeof(slong) * FLINT_MAX(q, 1));
    gr_mat_init(T, r, r, ctx);
    gr_mat_init(S, r, r, ctx);
    rho = gr_heap_init_vec(r, ctx);
    row = gr_heap_init_vec(r, ctx);
    GR_TMP_INIT(val, ctx);
    C.c = NULL;
    C.len = C.alloc = 0;
    C.p = p;
    C.q = q;
    C.limit = 64;
    C.ctx = ctx;
    E.r = r;
    E.len = 0;
    E.rows = flint_malloc(sizeof(gr_ptr) * r);
    E.vals = flint_malloc(sizeof(gr_ptr) * r);
    E.pcol = flint_malloc(sizeof(slong) * r);

    status = _hyp_cell(P0, n, red, P, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    {
        slong tot = 0;
        for (i = 0; i < p + q; i++)
            tot += FLINT_ABS(n[i]);
        if (tot > HYP_SHIFT_LIMIT)
        {
            status = GR_UNABLE;
            goto cleanup;
        }
    }

    /* the row of F(P) in the basis v(P0) = (theta^k F(P0)), then in the
       basis G_j = F(P0 + j e) */
    status = _hyp_path(T, P0, n, z, ctx);
    status |= _hyp_basis_matrix(S, P0, z, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;
    for (j = 0; j < r; j++)
        status |= gr_set(ENTRY(rho, j, sz), GR_MAT_ENTRY(T, 0, j, sz), ctx);
    status |= _hyp_row_mul(rho, S, ctx);

    /* the candidate relations */
    for (i = 0; i < p; i++)
    {
        /* a_i = 1 in the cell: a positive integer (once for equal ones) */
        slong m, k;
        int dup = 0;
        if (_hyp_integer(&m, HA(P0, i, sz), ctx) != 1)
            continue;
        for (k = 0; k < C.len; k++)
            if (C.c[k].type == 2)
                dup = 1;
        if (!dup)
        {
            hyp_cand_struct * cc = _hyp_cands_push(&C);
            cc->type = 2;
            cc->i = i;
            cc->j = -1;
            cc->dist = 0;
            cc->zconst = 0;
            cc->entry = -1;
        }
    }
    for (j = 0; j < q; j++)
    {
        if (red[j] >= 0)
        {
            hyp_cand_struct * cc = _hyp_cands_push(&C);
            cc->type = 0;
            cc->i = red[j];
            cc->j = j;
            cc->dist = 0;
            cc->zconst = 0;
            cc->entry = -1;
        }
    }
    if (status == GR_SUCCESS)
        status = _hyp_collect(&C, P0, red, z, 0, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;
    if (C.len > 1)
        qsort(C.c, C.len, sizeof(hyp_cand_struct), _hyp_cand_cmp);

    /* the relations, as long as they are independent and F(P) is not
       determined */
    for (c = 0; c < C.len && status == GR_SUCCESS && !done && E.len < r; c++)
    {
        hyp_cand_struct * cc = C.c + c;
        slong piv;

        status = _hyp_cand_row(row, val, cc, P0, S, r, z, ctx);
        if (status == GR_UNABLE && cc->type != 0 && cc->type != 2)
        {
            status = GR_SUCCESS;
            continue;
        }
        if (status != GR_SUCCESS)
            break;

        status = _hyp_rref_reduce(row, NULL, &E, ctx);
        piv = _hyp_rref_pivot(row, r, ctx);
        if (status != GR_SUCCESS || piv < 0)
            continue;

        /* an independent relation: its value */
        if (cc->type == 0)
        {
            hyp_params_t Pr;
            _hyp_params_init(Pr, p, q, ctx);
            status = _hyp_params_drop(Pr, P0, cc->i, cc->j, ctx);
            if (status == GR_SUCCESS)
                status = _hyp_eval(val, Pr, z, depth + 1, 0, ctx);
            if (status == GR_SUCCESS)
                status = gr_mul(val, val, HA(P0, cc->i, sz), ctx);
            _hyp_params_clear(Pr, ctx);
        }
        else if (cc->type == 2)
            status = _hyp_int_a_row(row, val, P0, cc->i, z, ctx);
        else
            status = _hyp_cand_value(val, cc, p, q, depth, ctx);

        if (status != GR_SUCCESS)
        {
            /* (an entry whose value cannot be computed: skipped) */
            status = GR_SUCCESS;
            continue;
        }

        /* (the row was reduced without the value: redo with it) */
        status = _hyp_cand_row(row, val, cc, P0, S, r, z, ctx);
        status |= _hyp_rref_reduce(row, val, &E, ctx);
        if (status == GR_SUCCESS)
            status = _hyp_rref_add(&E, row, val, piv, ctx);

        /* F(P) determined? */
        if (status == GR_SUCCESS)
        {
            status = _gr_vec_set(row, rho, r, ctx);
            status |= _hyp_rref_reduce(row, NULL, &E, ctx);
            if (status == GR_SUCCESS && _hyp_rref_pivot(row, r, ctx) == -1)
                done = 1;
        }
    }

    if (status != GR_SUCCESS)
        goto cleanup;

    /* F(P) = rho . G: the pivot columns through their relations, the
       others as generators */
    status = gr_zero(res, ctx);
    status |= gr_zero(val, ctx);
    status |= _gr_vec_set(row, rho, r, ctx);
    status |= _hyp_rref_reduce(row, val, &E, ctx);
    /* (now rho . G = -val + row . G with row zero in the pivot columns) */
    status |= gr_neg(res, val, ctx);

    for (j = 0; j < r && status == GR_SUCCESS; j++)
    {
        truth_t zero = gr_is_zero(ENTRY(row, j, sz), ctx);
        gr_ptr args, g;
        slong k;

        if (zero == T_TRUE)
            continue;
        if (zero == T_UNKNOWN)
        {
            status = GR_UNABLE;
            break;
        }

        args = gr_heap_init_vec(1 + p + q, ctx);
        GR_TMP_INIT(g, ctx);
        status |= gr_set(args, z, ctx);
        for (k = 0; k < p + q; k++)
            status |= gr_add_si(ENTRY(args, 1 + k, sz), ENTRY(P0->v, k, sz), j, ctx);
        /* (the parameters of the generator in the canonical order: the
           function is symmetric in the a_i and in the b_j) */
        _hyp_sort(ENTRY(args, 1, sz), p, ctx);
        _hyp_sort(ENTRY(args, 1 + p, sz), q, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_special_gen_multi_locked(g, args, 1 + p + q,
                GR_TOWER_HYPGEOM, GR_TOWER_HYPGEOM_PARAM(p, q), ctx);
        if (status == GR_SUCCESS)
        {
            status |= gr_mul(g, g, ENTRY(row, j, sz), ctx);
            status |= gr_add(res, res, g, ctx);
        }
        GR_TMP_CLEAR(g, ctx);
        gr_heap_clear_vec(args, 1 + p + q, ctx);
    }

cleanup:
    for (i = 0; i < E.len; i++)
    {
        gr_heap_clear_vec(E.rows[i], r, ctx);
        gr_heap_clear(E.vals[i], ctx);
    }
    flint_free(E.rows);
    flint_free(E.vals);
    flint_free(E.pcol);
    _hyp_cands_clear(&C, ctx);
    gr_heap_clear_vec(rho, r, ctx);
    gr_heap_clear_vec(row, r, ctx);
    GR_TMP_CLEAR(val, ctx);
    gr_mat_clear(T, ctx);
    gr_mat_clear(S, ctx);
    flint_free(n);
    flint_free(red);
    _hyp_params_clear(P0, ctx);
    return status;
}

/* the value at the point P itself from a table entry with a constant
   argument (summation theorems at z = 1, used without contiguity):
   1 if found */
static int
_hyp_exact(gr_ptr res, int * found, const hyp_params_t P, gr_srcptr z, int depth, gr_ctx_t ctx)
{
    hyp_cands_struct C;
    slong * red, i, j, n = P->p + P->q;
    int status;

    *found = 0;
    red = flint_malloc(sizeof(slong) * FLINT_MAX(P->q, 1));
    for (j = 0; j < P->q; j++)
        red[j] = -1;
    C.c = NULL;
    C.len = C.alloc = 0;
    C.p = P->p;
    C.q = P->q;
    C.limit = 16;
    C.ctx = ctx;

    status = _hyp_collect(&C, P, red, z, 1, ctx);
    for (i = 0; i < C.len && status == GR_SUCCESS && !*found; i++)
    {
        for (j = 0; j < n; j++)
            if (C.c[i].off[j] != 0)
                break;
        if (j == n && C.c[i].zconst && _hyp_cand_value(res, C.c + i, P->p, P->q, depth, ctx) == GR_SUCCESS)
            *found = 1;
    }

    _hyp_cands_clear(&C, ctx);
    flint_free(red);
    return status;
}

/* -------------------------------------------------------------------- */
/* layers 1 and 2                                                        */
/* -------------------------------------------------------------------- */

/* the terminating series sum_{k=0}^{N} (a)_k / (b)_k z^k / k! */
static int
_hyp_terminating(gr_ptr res, const hyp_params_t P, gr_srcptr z, slong N, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, i, k;
    gr_ptr t, u;
    int status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);
    status |= gr_one(t, ctx);
    status |= gr_one(res, ctx);
    for (k = 0; k < N && status == GR_SUCCESS; k++)
    {
        for (i = 0; i < P->p; i++)
        {
            status |= gr_add_si(u, HA(P, i, sz), k, ctx);
            status |= gr_mul(t, t, u, ctx);
        }
        for (i = 0; i < P->q; i++)
        {
            status |= gr_add_si(u, HB(P, i, sz), k, ctx);
            status |= gr_div(t, t, u, ctx);
        }
        status |= gr_mul(t, t, z, ctx);
        status |= gr_div_ui(t, t, k + 1, ctx);
        status |= gr_add(res, res, t, ctx);
    }
    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

/*
    The reduction of the order: a_i = b_j + m, m >= 1, F = (theta + b_j)_m
    F_red / (b_j)_m with F_red = F without a_i, b_j. With theta^k =
    sum_l S2(k, l) z^l D^l and D^l F_red = (a)_l / (b)_l F_red(P + l e)
    (exact for every k: iterating the companion matrix would need the
    derivatives of its coefficients), F = sum_l c_l F_red(P + l e) / (b_j)_m
    with c_l = z^l (a)_l / (b)_l sum_{k >= l} poly_k S2(k, l).
*/
static int
_hyp_reduce_order(gr_ptr res, const hyp_params_t P, slong i, slong j, slong m, gr_srcptr z, int depth, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, k, l, h;
    hyp_params_t R, Q;
    gr_ptr poly, f, t, b, c, g;
    fmpz_t s2;
    int status = GR_SUCCESS;

    _hyp_params_init(R, P->p, P->q, ctx);
    status |= _hyp_params_drop(R, P, i, j, ctx);
    _hyp_params_init(Q, R->p, R->q, ctx);

    poly = gr_heap_init_vec(m + 1, ctx);   /* (x + b)_m */
    f = gr_heap_init_vec(m + 1, ctx);      /* z^l (a)_l / (b)_l */
    GR_TMP_INIT4(t, b, c, g, ctx);
    fmpz_init(s2);

    status |= gr_one(poly, ctx);
    for (k = 0; k < m; k++)
    {
        status |= gr_add_si(b, HB(P, j, sz), k, ctx);
        for (l = k + 1; l >= 1; l--)
        {
            status |= gr_mul(t, ENTRY(poly, l, sz), b, ctx);
            status |= gr_add(ENTRY(poly, l, sz), t, ENTRY(poly, l - 1, sz), ctx);
        }
        status |= gr_mul(ENTRY(poly, 0, sz), ENTRY(poly, 0, sz), b, ctx);
    }

    status |= gr_one(f, ctx);
    for (l = 1; l <= m; l++)
    {
        status |= gr_mul(ENTRY(f, l, sz), ENTRY(f, l - 1, sz), z, ctx);
        for (h = 0; h < R->p; h++)
        {
            status |= gr_add_si(t, HA(R, h, sz), l - 1, ctx);
            status |= gr_mul(ENTRY(f, l, sz), ENTRY(f, l, sz), t, ctx);
        }
        for (h = 0; h < R->q; h++)
        {
            status |= gr_add_si(t, HB(R, h, sz), l - 1, ctx);
            status |= gr_div(ENTRY(f, l, sz), ENTRY(f, l, sz), t, ctx);
        }
    }

    status |= gr_zero(res, ctx);
    for (l = 0; l <= m && status == GR_SUCCESS; l++)
    {
        /* c_l = f_l sum_{k >= l} poly_k S2(k, l) */
        status |= gr_zero(c, ctx);
        for (k = l; k <= m; k++)
        {
            arith_stirling_number_2(s2, k, l);
            status |= gr_mul_fmpz(t, ENTRY(poly, k, sz), s2, ctx);
            status |= gr_add(c, c, t, ctx);
        }
        status |= gr_mul(c, c, ENTRY(f, l, sz), ctx);
        if (status != GR_SUCCESS || gr_is_zero(c, ctx) == T_TRUE)
            continue;

        for (h = 0; h < R->p + R->q; h++)
            status |= gr_add_si(ENTRY(Q->v, h, sz), ENTRY(R->v, h, sz), l, ctx);
        if (status == GR_SUCCESS)
            status = _hyp_eval(g, Q, z, depth + 1, 0, ctx);
        if (status == GR_SUCCESS)
        {
            status |= gr_mul(g, g, c, ctx);
            status |= gr_add(res, res, g, ctx);
        }
    }

    /* / (b)_m = poly(0) */
    if (status == GR_SUCCESS)
        status = gr_div(res, res, poly, ctx);

    gr_heap_clear_vec(poly, m + 1, ctx);
    gr_heap_clear_vec(f, m + 1, ctx);
    GR_TMP_CLEAR4(t, b, c, g, ctx);
    fmpz_clear(s2);
    _hyp_params_clear(R, ctx);
    _hyp_params_clear(Q, ctx);
    return status;
}

/* 1F0(a;; z) = (1 - z)^(-a) */
static int
_hyp_1f0(gr_ptr res, gr_srcptr a, gr_srcptr z, gr_ctx_t ctx)
{
    gr_ptr t, e;
    int status;

    GR_TMP_INIT2(t, e, ctx);
    status = gr_sub_ui(t, z, 1, ctx);
    status |= gr_neg(t, t, ctx);
    status |= gr_neg(e, a, ctx);
    if (status == GR_SUCCESS)
        status = _hyp_pow(res, t, e, ctx);
    GR_TMP_CLEAR2(t, e, ctx);
    return status;
}

/*
    The canonical member of the orbit of (P, z) under the transformations:
    for 2F1, Pfaff when |z - 1| > 1 (or |z - 1| = 1, Im(z) < 0), off the
    cut, then the smaller of P and its Euler transform; for 1F1, Kummer
    when Re(z) > 0 (or Re(z) = 0, Im(z) < 0). Sets *g to the
    transformation to apply first (HYP_T_ID if none), and *g2 to a second
    one (Euler after Pfaff).
*/
static int
_hyp_choose_transform(int * g, int * g2, const hyp_params_t P, gr_srcptr z, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    *g = *g2 = HYP_T_ID;

    if (P->p == 1 && P->q == 1)
    {
        int s;
        status = _gr_tower_lazy_re_cmp(&s, z, NULL, ctx);
        if (status == GR_SUCCESS && s == 0)
        {
            status = _gr_tower_lazy_im_sign(&s, z, ctx);
            s = -s;
        }
        if (status == GR_SUCCESS && s > 0)
            *g = HYP_T_KUMMER;
        return status;
    }

    if (P->p == 2 && P->q == 1 && _hyp_off_cut(z, ctx) == 1)
    {
        gr_ptr u, v;
        int sg = 0, im = 0;
        GR_TMP_INIT2(u, v, ctx);
        /* |z - 1|^2 - 1 */
        status = gr_sub_ui(u, z, 1, ctx);
        status |= gr_conj(v, u, ctx);
        status |= gr_mul(u, u, v, ctx);
        status |= gr_sub_ui(u, u, 1, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_re_cmp(&sg, u, NULL, ctx);
        if (status == GR_SUCCESS && sg == 0)
            status = _gr_tower_lazy_im_sign(&im, z, ctx);
        if (status == GR_SUCCESS && (sg > 0 || (sg == 0 && im < 0)))
            *g = HYP_T_PFAFF_A;
        GR_TMP_CLEAR2(u, v, ctx);
    }

    return status;
}

/* lexicographic comparison of the cells of two parameter vectors (sorted) */
static int
_hyp_cmp_cells(const hyp_params_t P, const hyp_params_t Q, gr_ctx_t ctx)
{
    hyp_params_t A, B;
    slong * n, * red, i, sz = ctx->sizeof_elem;
    int res = 0;

    _hyp_params_init(A, P->p, P->q, ctx);
    _hyp_params_init(B, Q->p, Q->q, ctx);
    n = flint_malloc(sizeof(slong) * (P->p + P->q));
    red = flint_malloc(sizeof(slong) * FLINT_MAX(P->q, 1));

    if (_hyp_cell(A, n, red, P, ctx) == GR_SUCCESS && _hyp_cell(B, n, red, Q, ctx) == GR_SUCCESS)
    {
        _hyp_params_sort(A, ctx);
        _hyp_params_sort(B, ctx);
        for (i = 0; i < P->p + P->q && res == 0; i++)
            res = _hyp_cmp(ENTRY(A->v, i, sz), ENTRY(B->v, i, sz), ctx);
    }

    flint_free(n);
    flint_free(red);
    _hyp_params_clear(A, ctx);
    _hyp_params_clear(B, ctx);
    return res;
}

static int
_hyp_eval(gr_ptr res, const hyp_params_t P_in, gr_srcptr z, int depth, int flags, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, i, j;
    hyp_params_t P;
    truth_t t;
    int status = GR_SUCCESS;

    if (depth > HYP_DEPTH_LIMIT)
        return GR_UNABLE;

    /* z = 0 */
    t = gr_is_zero(z, ctx);
    if (t == T_UNKNOWN)
        return GR_UNABLE;
    if (t == T_TRUE)
        return gr_one(res, ctx);

    _hyp_params_init(P, P_in->p, P_in->q, ctx);
    status |= _gr_vec_set(P->v, P_in->v, P_in->p + P_in->q, ctx);

    /* cancellation a_i = b_j */
    for (i = 0; i < P->p && status == GR_SUCCESS; i++)
    {
        for (j = 0; j < P->q; j++)
        {
            t = gr_equal(HA(P, i, sz), HB(P, j, sz), ctx);
            if (t == T_UNKNOWN)
            {
                status = GR_UNABLE;
                break;
            }
            if (t == T_TRUE)
            {
                status = _hyp_params_drop(P, P, i, j, ctx);
                i = -1;   /* (restart) */
                break;
            }
        }
    }
    if (status != GR_SUCCESS)
        goto cleanup;

    /* 0F0 = exp(z), 1F0 = (1 - z)^(-a) */
    if (P->p == 0 && P->q == 0)
    {
        status = gr_exp(res, z, ctx);
        goto cleanup;
    }
    if (P->p == 1 && P->q == 0)
    {
        status = _hyp_1f0(res, HA(P, 0, sz), z, ctx);
        goto cleanup;
    }

    /* terminating series and poles */
    {
        slong N = WORD_MAX, poleb = WORD_MAX, m;
        for (i = 0; i < P->p; i++)
        {
            int r = _hyp_integer(&m, HA(P, i, sz), ctx);
            if (r == -1)
            {
                status = GR_UNABLE;
                goto cleanup;
            }
            if (r == 1 && m <= 0)
                N = FLINT_MIN(N, -m);
        }
        for (j = 0; j < P->q; j++)
        {
            int r = _hyp_integer(&m, HB(P, j, sz), ctx);
            if (r == -1)
            {
                status = GR_UNABLE;
                goto cleanup;
            }
            if (r == 1 && m <= 0)
                poleb = FLINT_MIN(poleb, -m);
        }
        if (N != WORD_MAX && N <= poleb)
        {
            /* (a pole at the term N when b_j = -N is reached at the same
               time: (a)_N / (b)_N is then 0/0, taken as the limit along
               the terminating direction, the usual convention, only when
               a_i is closer to zero) */
            if (N == poleb)
            {
                status = GR_DOMAIN;
                goto cleanup;
            }
            if (N > HYP_TERM_LIMIT)
            {
                status = GR_UNABLE;
                goto cleanup;
            }
            status = _hyp_terminating(res, P, z, N, ctx);
            goto cleanup;
        }
        if (poleb != WORD_MAX)
        {
            status = GR_DOMAIN;
            goto cleanup;
        }
    }

    /* divergent series */
    if (P->p > P->q + 1)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* z = 1 for p = q + 1 */
    if (P->p == P->q + 1)
    {
        t = gr_is_one(z, ctx);
        if (t == T_UNKNOWN)
        {
            status = GR_UNABLE;
            goto cleanup;
        }
        if (t == T_TRUE)
        {
            gr_ptr s, u;
            int sgn;
            GR_TMP_INIT2(s, u, ctx);
            status = gr_zero(s, ctx);
            for (j = 0; j < P->q; j++)
                status |= gr_add(s, s, HB(P, j, sz), ctx);
            for (i = 0; i < P->p; i++)
                status |= gr_sub(s, s, HA(P, i, sz), ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_re_cmp(&sgn, s, NULL, ctx);
            if (status == GR_SUCCESS && sgn <= 0)
                status = GR_DOMAIN;    /* (divergent at z = 1) */
            if (status == GR_SUCCESS && P->p == 2)
            {
                /* Gauss: Gamma(c) Gamma(c - a - b) / (Gamma(c - a) Gamma(c - b)) */
                status = gr_gamma(res, HB(P, 0, sz), ctx);
                status |= gr_gamma(u, s, ctx);
                status |= gr_mul(res, res, u, ctx);
                status |= gr_sub(u, HB(P, 0, sz), HA(P, 0, sz), ctx);
                status |= gr_rgamma(u, u, ctx);
                status |= gr_mul(res, res, u, ctx);
                status |= gr_sub(u, HB(P, 0, sz), HA(P, 1, sz), ctx);
                status |= gr_rgamma(u, u, ctx);
                status |= gr_mul(res, res, u, ctx);
            }
            else if (status == GR_SUCCESS)
            {
                /* the table at the point itself (no contiguity at z = 1) */
                int found;
                status = _hyp_exact(res, &found, P, z, depth, ctx);
                if (status == GR_SUCCESS && !found)
                {
                    /* a generator (p = q + 1 >= 3: not yet, the
                       enclosures of _gr_tower_special_eval_multi being
                       for |z| < 1 only, so that the adjunction fails with
                       GR_UNABLE) */
                    gr_ptr args;
                    slong k;
                    _hyp_params_sort(P, ctx);
                    args = gr_heap_init_vec(1 + P->p + P->q, ctx);
                    status |= gr_set(args, z, ctx);
                    for (k = 0; k < P->p + P->q; k++)
                        status |= gr_set(ENTRY(args, 1 + k, sz), ENTRY(P->v, k, sz), ctx);
                    if (status == GR_SUCCESS)
                        status = _gr_tower_lazy_special_gen_multi_locked(res, args, 1 + P->p + P->q,
                            GR_TOWER_HYPGEOM, GR_TOWER_HYPGEOM_PARAM(P->p, P->q), ctx);
                    gr_heap_clear_vec(args, 1 + P->p + P->q, ctx);
                }
            }
            GR_TMP_CLEAR2(s, u, ctx);
            goto cleanup;
        }
    }

    /* the reduction of the order: a_i - b_j a positive integer */
    for (i = 0; i < P->p; i++)
    {
        for (j = 0; j < P->q; j++)
        {
            gr_ptr d;
            slong m;
            int r;
            GR_TMP_INIT(d, ctx);
            status = gr_sub(d, HA(P, i, sz), HB(P, j, sz), ctx);
            r = (status == GR_SUCCESS) ? _hyp_integer(&m, d, ctx) : -1;
            GR_TMP_CLEAR(d, ctx);
            if (r == -1)
            {
                status = GR_UNABLE;
                goto cleanup;
            }
            if (r == 1 && m >= 1)
            {
                if (m > HYP_REDUCE_LIMIT)
                    status = GR_UNABLE;
                else
                    status = _hyp_reduce_order(res, P, i, j, m, z, depth, ctx);
                goto cleanup;
            }
        }
    }

    /* the transformations of one term */
    if (!(flags & HYP_NO_ORBIT))
    {
        int g, g2;
        status = _hyp_choose_transform(&g, &g2, P, z, ctx);
        if (status != GR_SUCCESS)
            goto cleanup;

        if (g != HYP_T_ID || (P->p == 2 && P->q == 1 && _hyp_off_cut(z, ctx) == 1))
        {
            hyp_params_t Q;
            gr_ptr w, f;
            _hyp_params_init(Q, P->p, P->q, ctx);
            GR_TMP_INIT2(w, f, ctx);
            status = _hyp_transform(Q, w, f, g, P, z, ctx);

            /* 2F1: the smaller of Q and its Euler transform */
            if (status == GR_SUCCESS && P->p == 2 && P->q == 1)
            {
                hyp_params_t E;
                gr_ptr w2, f2;
                _hyp_params_init(E, 2, 1, ctx);
                GR_TMP_INIT2(w2, f2, ctx);
                status = _hyp_transform(E, w2, NULL, HYP_T_EULER, Q, w, ctx);
                if (status == GR_SUCCESS && _hyp_cmp_cells(E, Q, ctx) < 0)
                {
                    status = _hyp_transform(E, w2, f2, HYP_T_EULER, Q, w, ctx);
                    status |= _hyp_params_set(Q, E, ctx);
                    status |= gr_mul(f, f, f2, ctx);
                }
                GR_TMP_CLEAR2(w2, f2, ctx);
                _hyp_params_clear(E, ctx);
            }

            if (status == GR_SUCCESS)
                status = _hyp_eval(res, Q, w, depth, HYP_NO_ORBIT, ctx);
            if (status == GR_SUCCESS)
                status = gr_mul(res, res, f, ctx);

            GR_TMP_CLEAR2(w, f, ctx);
            _hyp_params_clear(Q, ctx);
            goto cleanup;
        }
    }

    /* summation theorems at the point itself (no linear algebra) */
    {
        int found;
        status = _hyp_exact(res, &found, P, z, depth, ctx);
        if (status == GR_SUCCESS && found)
            goto cleanup;
    }

    if (status == GR_SUCCESS)
        status = _hyp_coset(res, P, z, depth, ctx);

cleanup:
    _hyp_params_clear(P, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* interface                                                             */
/* -------------------------------------------------------------------- */

/* regularized: F / prod Gamma(b_j), finite at b_j = -n (see the docs) */
static int
_hyp_eval_flags(gr_ptr res, const hyp_params_t P, gr_srcptr z, int regularized, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem, j;
    int status;

    if (!regularized)
        return _hyp_eval(res, P, z, 0, 0, ctx);

    /* b_j = -n: F~ = prod (a)_{n+1} z^(n+1) F~(a + n + 1; n + 2, b + n + 1) */
    for (j = 0; j < P->q; j++)
    {
        slong m;
        int r = _hyp_integer(&m, HB(P, j, sz), ctx);
        if (r == -1)
            return GR_UNABLE;
        if (r == 1 && m <= 0)
        {
            slong n = -m, k, i;
            hyp_params_t Q;
            gr_ptr f, t;

            if (n > HYP_TERM_LIMIT)
                return GR_UNABLE;

            _hyp_params_init(Q, P->p, P->q, ctx);
            GR_TMP_INIT2(f, t, ctx);
            status = gr_pow_ui(f, z, n + 1, ctx);
            for (i = 0; i < P->p; i++)
            {
                status |= gr_rising_ui(t, HA(P, i, sz), n + 1, ctx);
                status |= gr_mul(f, f, t, ctx);
                status |= gr_add_ui(HA(Q, i, sz), HA(P, i, sz), n + 1, ctx);
            }
            for (k = 0; k < P->q; k++)
            {
                if (k == j)
                    status |= gr_set_ui(HB(Q, k, sz), n + 2, ctx);
                else
                    status |= gr_add_ui(HB(Q, k, sz), HB(P, k, sz), n + 1, ctx);
            }
            if (status == GR_SUCCESS)
                status = _hyp_eval_flags(res, Q, z, 1, ctx);
            if (status == GR_SUCCESS)
                status = gr_mul(res, res, f, ctx);
            GR_TMP_CLEAR2(f, t, ctx);
            _hyp_params_clear(Q, ctx);
            return status;
        }
    }

    status = _hyp_eval(res, P, z, 0, 0, ctx);
    for (j = 0; j < P->q && status == GR_SUCCESS; j++)
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status = gr_rgamma(t, HB(P, j, sz), ctx);
        status |= gr_mul(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);
    }
    return status;
}

int
gr_tower_lazy_hypgeom_pfq_vec(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_struct * a, slong p,
    const gr_tower_lazy_elem_struct * b, slong q, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
{
    hyp_params_t P;
    slong sz = ctx->sizeof_elem, i;
    int status = GR_SUCCESS, real, alg;
    gr_ptr zz;

    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);

    _hyp_params_init(P, p, q, ctx);
    GR_TMP_INIT(zz, ctx);
    for (i = 0; i < p; i++)
        status |= gr_set(HA(P, i, sz), a + i, ctx);
    for (i = 0; i < q; i++)
        status |= gr_set(HB(P, i, sz), b + i, ctx);
    status |= gr_set(zz, z, ctx);

    if (status == GR_SUCCESS)
        status = _hyp_eval_flags(res, P, zz, flags & 1, ctx);

    GR_TMP_CLEAR(zz, ctx);
    _hyp_params_clear(P, ctx);

    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_hypgeom_pfq(gr_tower_lazy_elem_t res, const gr_vec_t a, const gr_vec_t b, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
{
    return gr_tower_lazy_hypgeom_pfq_vec(res, a->entries, a->length, b->entries, b->length, z, flags, ctx);
}

int
gr_tower_lazy_hypgeom_0f1(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t b, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
{
    return gr_tower_lazy_hypgeom_pfq_vec(res, NULL, 0, b, 1, z, flags, ctx);
}

int
gr_tower_lazy_hypgeom_1f1(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t a, const gr_tower_lazy_elem_t b, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
{
    return gr_tower_lazy_hypgeom_pfq_vec(res, a, 1, b, 1, z, flags, ctx);
}

int
gr_tower_lazy_hypgeom_2f1(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t a, const gr_tower_lazy_elem_t b,
    const gr_tower_lazy_elem_t c, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct ab[2];
    int status;
    /* (not under the lock: the view restrictions of the context apply in
       gr_tower_lazy_hypgeom_pfq_vec only at the outermost level) */
    gr_tower_lazy_init(ab, ctx);
    gr_tower_lazy_init(ab + 1, ctx);
    status = gr_set(ab, a, ctx);
    status |= gr_set(ab + 1, b, ctx);
    if (status == GR_SUCCESS)
        status = gr_tower_lazy_hypgeom_pfq_vec(res, ab, 2, c, 1, z, flags, ctx);
    gr_tower_lazy_clear(ab, ctx);
    gr_tower_lazy_clear(ab + 1, ctx);
    return status;
}

/* the value of the generator definition (z, a, b): non-regularized, with
   canonicalization (conjugation, transfer between contexts) */
int
_gr_tower_lazy_hypgeom_args(gr_tower_lazy_elem_t res, slong param, const gr_tower_lazy_elem_struct * args, slong nargs, gr_ctx_t ctx)
{
    slong p = GR_TOWER_HYPGEOM_P(param), q = GR_TOWER_HYPGEOM_Q(param);
    hyp_params_t P;
    slong sz = ctx->sizeof_elem, i;
    int status = GR_SUCCESS;

    if (nargs != 1 + p + q)
        return GR_UNABLE;

    _gr_tower_lazy_lock(ctx);
    _hyp_params_init(P, p, q, ctx);
    for (i = 0; i < p + q; i++)
        status |= gr_set(ENTRY(P->v, i, sz), args + 1 + i, ctx);
    if (status == GR_SUCCESS)
        status = _hyp_eval(res, P, args, 0, 0, ctx);
    _hyp_params_clear(P, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* K(m) = (pi/2) 2F1(1/2, 1/2; 1; m), E(m) = (pi/2) 2F1(-1/2, 1/2; 1; m) */
int
_gr_tower_lazy_elliptic_hypgeom(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t m, int kind, gr_ctx_t ctx)
{
    hyp_params_t P;
    slong sz = ctx->sizeof_elem;
    gr_ptr t;
    int status;

    _gr_tower_lazy_lock(ctx);
    _hyp_params_init(P, 2, 1, ctx);
    GR_TMP_INIT(t, ctx);
    status = gr_set_si(HA(P, 0, sz), (kind == GR_TOWER_ELLIPTIC_K) ? 1 : -1, ctx);
    status |= gr_div_ui(HA(P, 0, sz), HA(P, 0, sz), 2, ctx);
    status |= gr_one(HA(P, 1, sz), ctx);
    status |= gr_div_ui(HA(P, 1, sz), HA(P, 1, sz), 2, ctx);
    status |= gr_one(HB(P, 0, sz), ctx);
    status |= gr_set(t, m, ctx);
    if (status == GR_SUCCESS)
        status = _hyp_eval(res, P, t, 0, 0, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_pi(t, ctx);
        status |= gr_div_ui(t, t, 2, ctx);
        status |= gr_mul(res, res, t, ctx);
    }
    GR_TMP_CLEAR(t, ctx);
    _hyp_params_clear(P, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* testing                                                               */
/* -------------------------------------------------------------------- */

/*
    For the tests: the shape of the table entry i (*p, *q, the symbols
    used in used[0..24], *zfree, the constant argument z0), and, with
    vals (HYP_NVARS values: the symbols a, ..., y and z = vals[25]), the
    parameters at those symbols (params: p + q elements), whether the
    conditions hold (*cond: 1, 0 or -1) and the value. Returns GR_UNABLE
    for an invalid entry.
*/
int
_gr_tower_hypgeom_table_shape(slong * p, slong * q, int * used, int * zfree, fmpq_t z0, slong i)
{
    hyp_entry_struct R;
    slong k;
    int ok = _hyp_entry_parse(&R, _hyp_table[i]);
    if (ok)
    {
        *p = R.p;
        *q = R.q;
        for (k = 0; k < HYP_NSYM - 1; k++)
            used[k] = R.used[k];
        *zfree = R.z_free;
        fmpq_set(z0, R.z0);
    }
    _hyp_entry_clear(&R);
    return ok ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_tower_hypgeom_table_eval(gr_ptr params, int * cond, gr_ptr value, slong i, gr_srcptr * vals, gr_ctx_t ctx)
{
    hyp_entry_struct R;
    slong k, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (!_hyp_entry_parse(&R, _hyp_table[i]))
    {
        _hyp_entry_clear(&R);
        return GR_UNABLE;
    }

    _gr_tower_lazy_lock(ctx);
    for (k = 0; k < R.p + R.q && status == GR_SUCCESS; k++)
        status = _hyp_aff_eval(ENTRY(params, k, sz), R.par + k, vals, ctx);
    if (status == GR_SUCCESS)
    {
        *cond = _hyp_entry_conditions(&R, vals, 0, ctx);
        if (*cond == 1)
            status = _hyp_entry_value(value, &R, vals, 1, ctx);
    }
    _gr_tower_lazy_unlock(ctx);

    _hyp_entry_clear(&R);
    return status;
}

POP_OPTIONS
