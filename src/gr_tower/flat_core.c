/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Core of the flat representation: multivariate rational functions in
    the generators of a tower, reduced modulo the triangular set with
    denominators cleared. See gr_tower.h for the variable layout.
*/

#include <string.h>
#include <stdio.h>
#include "ulong_extras.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_q.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

static void
_flat_clear_ideal(gr_tower_flat_t F)
{
    slong i;

    for (i = 0; i < F->ideal_len; i++)
    {
        if (F->ideal[i] != NULL)
        {
            fmpz_mpoly_clear(F->ideal[i], F->mctx);
            flint_free(F->ideal[i]);
        }
        if (F->ideal_univar[i] != NULL)
        {
            fmpz_poly_clear(F->ideal_univar[i]);
            flint_free(F->ideal_univar[i]);
        }
    }
    flint_free(F->ideal);
    flint_free(F->ideal_univar);
    F->ideal = NULL;
    F->ideal_univar = NULL;
    F->ideal_len = 0;
}

/* var_gid for the current definition order and the given capacity */
static slong *
_make_var_gid(gr_tower_t T, slong cap)
{
    slong * var_gid = flint_malloc(sizeof(slong) * cap);
    slong d;
    for (d = 0; d < cap; d++)
        var_gid[d] = -1;
    for (d = 0; d < T->num_gens; d++)
        var_gid[cap - 1 - d] = T->gens[d].gid;
    return var_gid;
}

void
gr_tower_flat_init(gr_tower_flat_t F, gr_tower_t T, slong cap)
{
    cap = FLINT_MAX(cap, FLINT_MAX(T->num_gens, 1));

    F->T = T;
    F->cap = cap;
    F->mctx = flint_malloc(sizeof(fmpz_mpoly_ctx_struct));
    fmpz_mpoly_ctx_init(F->mctx, cap, ORD_LEX);
    F->var_gid = _make_var_gid(T, cap);
    F->ideal = NULL;
    F->ideal_univar = NULL;
    F->ideal_len = 0;
    F->ideal_version = 0;
    F->ideal_moduli_version = 0;
    F->ideal_lc_const = 1;
    F->old_mctx = NULL;
    F->old_var_gid = NULL;
    F->num_old_mctx = 0;
    F->layout_version = T->structure_version;
    F->stale = NULL;
    F->stale_mctx = NULL;
    F->num_stale = 0;
    F->inv_cache = NULL;
    F->inv_cache_mctx = NULL;
    F->modp = NULL;
    F->refs = 0;
    F->gc = 0;
    F->dense_fields = NULL;
    F->primitive_tried[0] = F->primitive_tried[1] = 0;
    F->rebase_key[0] = F->rebase_key[1] = 0;
    F->rebase_pending = 0;
}

void
_gr_tower_dense_fields_clear(gr_tower_flat_t F)
{
    _gr_tower_dense_field_struct * D = F->dense_fields, * next;
    while (D != NULL)
    {
        next = D->next;
        fmpz_poly_clear(&D->modulus);
        flint_free(D);
        D = next;
    }
    F->dense_fields = NULL;
}

/* Allocates an element (in the current context) owned by F. */
fmpz_mpoly_q_struct *
_gr_tower_flat_stale_alloc(gr_tower_flat_t F)
{
    fmpz_mpoly_q_struct * p = flint_malloc(sizeof(fmpz_mpoly_q_struct));
    fmpz_mpoly_q_init(p, F->mctx);
    F->stale = flint_realloc(F->stale, sizeof(fmpz_mpoly_q_struct *) * (F->num_stale + 1));
    F->stale_mctx = flint_realloc(F->stale_mctx, sizeof(fmpz_mpoly_ctx_struct *) * (F->num_stale + 1));
    F->stale[F->num_stale] = p;
    F->stale_mctx[F->num_stale] = F->mctx;
    F->num_stale++;
    return p;
}

void
gr_tower_flat_clear(gr_tower_flat_t F)
{
    slong i;

    for (i = 0; i < F->num_stale; i++)
    {
        fmpz_mpoly_q_clear(F->stale[i], F->stale_mctx[i]);
        flint_free(F->stale[i]);
    }
    flint_free(F->stale);
    flint_free(F->stale_mctx);

    _gr_tower_flat_modp_clear(F);
    _gr_tower_dense_fields_clear(F);
    if (F->inv_cache != NULL)
    {
        fmpz_mpoly_q_clear(F->inv_cache, F->inv_cache_mctx);
        fmpz_mpoly_q_clear(F->inv_cache + 1, F->inv_cache_mctx);
        flint_free(F->inv_cache);
    }

    _flat_clear_ideal(F);
    fmpz_mpoly_ctx_clear(F->mctx);
    flint_free(F->mctx);

    flint_free(F->var_gid);
    for (i = 0; i < F->num_old_mctx; i++)
    {
        fmpz_mpoly_ctx_clear(F->old_mctx[i]);
        flint_free(F->old_mctx[i]);
        flint_free(F->old_var_gid[i]);
    }
    flint_free(F->old_mctx);
    flint_free(F->old_var_gid);
}

/* Layout of the current or a superseded context of F. */
int
_gr_tower_flat_find_layout(slong * cap, slong ** var_gid, const fmpz_mpoly_ctx_struct * mctx, gr_tower_flat_t F)
{
    slong i;

    if (mctx == F->mctx)
    {
        *cap = F->cap;
        *var_gid = F->var_gid;
        return 1;
    }

    for (i = 0; i < F->num_old_mctx; i++)
    {
        if (F->old_mctx[i] == mctx)
        {
            *cap = mctx->minfo->nvars;
            *var_gid = F->old_var_gid[i];
            return 1;
        }
    }
    return 0;
}

/* Renames variables: c[v] = target variable of the generator at variable v
   of a layout of capacity old_cap; the target is the current context of
   dst, by definition order in dst (which must contain the generators:
   either dst->T is the same tower (by gid) or a tower sharing the
   definition order of the generators occurring (by definition order in
   src, when src != NULL)). */
static void
_rename(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t old_mctx, slong old_cap, const slong * old_var_gid, gr_tower_flat_t src, gr_tower_flat_t dst, const slong * order_map)
{
    slong * c = flint_malloc(sizeof(slong) * old_cap);
    fmpz_mpoly_t t;
    slong v;

    for (v = 0; v < old_cap; v++)
    {
        slong g = old_var_gid[v], d;
        c[v] = 0;
        if (g < 0)
            continue;
        d = gr_tower_gid_order(src == NULL ? dst->T : src->T, g);
        if (d >= 0 && order_map != NULL)
            d = order_map[d];
        if (d >= 0 && d < dst->T->num_gens)
            c[v] = GR_TOWER_FLAT_VAR_D(dst, d);
    }

    fmpz_mpoly_init(t, dst->mctx);
    fmpz_mpoly_compose_fmpz_mpoly_gen(t, fmpz_mpoly_q_numref(x), c, old_mctx, dst->mctx);
    fmpz_mpoly_swap(fmpz_mpoly_q_numref(res), t, dst->mctx);
    fmpz_mpoly_compose_fmpz_mpoly_gen(t, fmpz_mpoly_q_denref(x), c, old_mctx, dst->mctx);
    fmpz_mpoly_swap(fmpz_mpoly_q_denref(res), t, dst->mctx);
    fmpz_mpoly_clear(t, dst->mctx);

    flint_free(c);
}

void
gr_tower_flat_convert(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t old_mctx, gr_tower_flat_t F)
{
    slong old_cap;
    slong * old_var_gid;

    _gr_tower_flat_grow(F);

    if (old_mctx == F->mctx)
    {
        fmpz_mpoly_q_set(res, x, F->mctx);
        return;
    }

    if (!_gr_tower_flat_find_layout(&old_cap, &old_var_gid, old_mctx, F))
        flint_throw(FLINT_ERROR, "(%s): unknown polynomial context\n", __func__);

    _rename(res, x, old_mctx, old_cap, old_var_gid, NULL, F, NULL);
}

/* Transports a flat element between the (current) contexts of two flat
   machineries whose towers list the generators occurring in x at the
   same definition orders. */
void
_gr_tower_flat_transport(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t src, gr_tower_flat_t dst)
{
    _gr_tower_flat_grow(src);
    _gr_tower_flat_grow(dst);
    _rename(res, x, src->mctx, src->cap, src->var_gid, src, dst, NULL);
}

/* The same, the generator of definition order d in src being the
   generator of definition order order_map[d] in dst. */
void
_gr_tower_flat_transport_map(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t src, gr_tower_flat_t dst, const slong * order_map)
{
    _gr_tower_flat_grow(src);
    _gr_tower_flat_grow(dst);
    _rename(res, x, src->mctx, src->cap, src->var_gid, src, dst, order_map);
}

/* Number of variables of the base field context (0 for QQ). */
static slong
_base_nvars(gr_tower_flat_t F)
{
    slong n = 0;
    if (F->T->base->which_ring == GR_CTX_FMPZ_MPOLY_Q)
        GR_MUST_SUCCEED(gr_ctx_ngens(&n, F->T->base));
    return n;
}

/* base element (constant or rational function in t's) -> flat element */
static int
_base_to_flat(fmpz_mpoly_q_t res, gr_srcptr x, gr_tower_flat_t F)
{
    gr_ctx_struct * base = F->T->base;

    if (base->which_ring == GR_CTX_FMPQ)
    {
        fmpz_mpoly_q_set_fmpq(res, x, F->mctx);
        return GR_SUCCESS;
    }
    else if (base->which_ring == GR_CTX_FMPZ_MPOLY_Q)
    {
        const fmpz_mpoly_q_struct * q = x;
        slong r = _base_nvars(F), j;
        slong * c = flint_malloc(sizeof(slong) * FLINT_MAX(r, 1));
        fmpz_mpoly_t t;
        const fmpz_mpoly_ctx_struct * bctx = gr_ctx_fmpz_mpoly_q_mctx(base);

        for (j = 0; j < r; j++)
            c[j] = GR_TOWER_FLAT_TVAR(F, j + 1);

        fmpz_mpoly_init(t, F->mctx);
        fmpz_mpoly_compose_fmpz_mpoly_gen(t, fmpz_mpoly_q_numref(q), c, bctx, F->mctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_numref(res), t, F->mctx);
        fmpz_mpoly_compose_fmpz_mpoly_gen(t, fmpz_mpoly_q_denref(q), c, bctx, F->mctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_denref(res), t, F->mctx);
        fmpz_mpoly_clear(t, F->mctx);
        flint_free(c);
        return GR_SUCCESS;
    }
    else
    {
        fmpq_t c;
        int status;
        fmpq_init(c);
        status = gr_get_fmpq(c, x, base);
        if (status == GR_SUCCESS)
            fmpz_mpoly_q_set_fmpq(res, c, F->mctx);
        fmpq_clear(c);
        return status;
    }
}

/* polynomial in the transcendental variables (in F's context) -> base element */
static int
_flat_to_base(gr_ptr res, const fmpz_mpoly_t f, gr_tower_flat_t F)
{
    gr_ctx_struct * base = F->T->base;

    if (base->which_ring == GR_CTX_FMPQ)
    {
        fmpz_t c;
        int status;
        if (!fmpz_mpoly_is_fmpz(f, F->mctx))
            return GR_DOMAIN;
        fmpz_init(c);
        fmpz_mpoly_get_fmpz(c, f, F->mctx);
        status = gr_set_fmpz(res, c, base);
        fmpz_clear(c);
        return status;
    }
    else
    {
        fmpz_mpoly_q_struct * q = res;
        slong r = _base_nvars(F), j, tot = F->cap;
        slong * c = flint_malloc(sizeof(slong) * tot);
        int * is_trans = flint_calloc(tot, sizeof(int));
        const fmpz_mpoly_ctx_struct * bctx = gr_ctx_fmpz_mpoly_q_mctx(base);

        /* only transcendental variables may occur */
        for (j = 0; j < tot; j++)
            c[j] = 0;
        for (j = 1; j <= r; j++)
        {
            c[GR_TOWER_FLAT_TVAR(F, j)] = j - 1;
            is_trans[GR_TOWER_FLAT_TVAR(F, j)] = 1;
        }

        {
            int * used = flint_malloc(sizeof(int) * tot);
            fmpz_mpoly_used_vars(used, f, F->mctx);
            for (j = 0; j < tot; j++)
            {
                if (used[j] && !is_trans[j])
                {
                    flint_free(used);
                    flint_free(is_trans);
                    flint_free(c);
                    return GR_DOMAIN;
                }
            }
            flint_free(used);
        }
        flint_free(is_trans);

        fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_numref(q), f, c, F->mctx, bctx);
        fmpz_mpoly_one(fmpz_mpoly_q_denref(q), bctx);
        flint_free(c);
        return GR_SUCCESS;
    }
}

/* nested element x of F_k -> flat element */
/* lcm of the denominators of the rational coordinates of x in F_k */
static void
_nested_denominator_lcm(fmpz_t den, gr_srcptr x, slong k, gr_tower_t T)
{
    if (k == 0)
    {
        const fmpz * d = fmpq_denref((const fmpq *) x);
        /* (the denominators repeat a lot: divisibility is cheaper than
           the lcm, and skipping the units cheaper still) */
        if (!fmpz_is_one(d) && !fmpz_divisible(den, d))
            fmpz_lcm(den, den, d);
    }
    else
    {
        const gr_poly_struct * poly = x;
        slong i;
        for (i = 0; i < poly->length; i++)
            _nested_denominator_lcm(den, gr_poly_coeff_srcptr(poly, i, gr_tower_field_at(T, k - 1)), k - 1, T);
    }
}

/* Pushes the terms of den * x (x in F_k, packed holding the packed
   exponents of the generators above k) onto res, whose exponent field
   width bits (at most FLINT_BITS) accommodates the degrees. The terms
   come out in decreasing lex order (the variable of a higher generator
   is more significant), distinct. */
#define COFACTOR_CACHE 8

/* den / d for the denominators d met recently (the leaves of a nested
   element share few distinct denominators, and den is large) */
typedef struct
{
    fmpz d[COFACTOR_CACHE];
    fmpz q[COFACTOR_CACHE];
    slong len, next;
}
cofactor_cache_struct;

/* returns a pointer to den / d (t being scratch) */
static const fmpz *
_cofactor(fmpz_t t, cofactor_cache_struct * cache, const fmpz_t den, const fmpz_t d)
{
    slong i;

    if (fmpz_is_one(d))
        return den;

    for (i = 0; i < cache->len; i++)
    {
        if (fmpz_equal(cache->d + i, d))
            return cache->q + i;
    }

    fmpz_divexact(t, den, d);
    i = cache->next;
    fmpz_set(cache->d + i, d);
    fmpz_swap(cache->q + i, t);
    cache->next = (i + 1) % COFACTOR_CACHE;
    cache->len = FLINT_MIN(cache->len + 1, COFACTOR_CACHE);
    return cache->q + i;
}

static void
_nested_push_terms(fmpz_mpoly_t res, gr_srcptr x, slong k, const ulong * packed, const ulong * gen_packed, slong N, const fmpz_t den, cofactor_cache_struct * cache, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;

    if (k == 0)
    {
        const fmpq * c = x;
        if (!fmpq_is_zero(c))
        {
            fmpz_t t;
            slong len = res->length;
            fmpz_init(t);
            fmpz_mul(t, _cofactor(t, cache, den, fmpq_denref(c)), fmpq_numref(c));
            fmpz_mpoly_fit_length(res, len + 1, F->mctx);
            mpoly_monomial_set(res->exps + N * len, packed, N);
            fmpz_swap(res->coeffs + len, t);
            res->length = len + 1;
            fmpz_clear(t);
        }
    }
    else
    {
        const gr_poly_struct * poly = x;
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        slong i, v = GR_TOWER_FLAT_VAR(F, k);
        ulong * cur = flint_malloc(sizeof(ulong) * N);
        for (i = poly->length - 1; i >= 0; i--)
        {
            mpoly_monomial_madd(cur, packed, i, gen_packed + N * v, N);
            _nested_push_terms(res, gr_poly_coeff_srcptr(poly, i, below), k - 1, cur, gen_packed, N, den, cache, F);
        }
        flint_free(cur);
    }
}

/* Largest length of a polynomial in the nested element x of F_k, per level. */
static void
_nested_max_lengths(slong * maxlen, gr_srcptr x, slong k, gr_tower_t T)
{
    if (k > 0)
    {
        const gr_poly_struct * poly = x;
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        slong i;
        maxlen[k] = FLINT_MAX(maxlen[k], poly->length);
        for (i = 0; i < poly->length; i++)
            _nested_max_lengths(maxlen, gr_poly_coeff_srcptr(poly, i, below), k - 1, T);
    }
}

/*
    Nested element of F_k -> flat element, term by term (the nested
    representation is sparse in the flat variables). Over a base of
    rational functions, the coordinates are combined as rational
    functions.
*/
static int
_nested_to_flat(fmpz_mpoly_q_t res, gr_srcptr x, slong k, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    gr_ctx_struct * F0 = T->base;
    int status = GR_SUCCESS;

    if (F0->which_ring == GR_CTX_FMPQ)
    {
        fmpz_t den;
        ulong * exp;

        fmpz_init(den);
        fmpz_one(den);
        _nested_denominator_lcm(den, x, k, T);

        /* exponent field width from the lengths of the polynomials */
        {
            slong * maxlen = flint_calloc(k + 1, sizeof(slong));
            fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(res);
            flint_bitcnt_t bits;
            slong N, j, v;
            ulong * gen_packed, * zero;

            _nested_max_lengths(maxlen, x, k, T);
            exp = flint_calloc(F->cap, sizeof(ulong));
            for (j = 1; j <= k; j++)
                exp[GR_TOWER_FLAT_VAR(F, j)] = FLINT_MAX(maxlen[j] - 1, 0);
            bits = mpoly_exp_bits_required_ui(exp, F->mctx->minfo);
            bits = mpoly_fix_bits(bits, F->mctx->minfo);
            if (bits > FLINT_BITS)
                bits = FLINT_BITS;   /* not reached: degrees are small */

            fmpz_mpoly_zero(num, F->mctx);
            fmpz_mpoly_fit_bits(num, bits, F->mctx);
            num->bits = bits;
            N = mpoly_words_per_exp(bits, F->mctx->minfo);
            gen_packed = flint_calloc(N * F->cap, sizeof(ulong));
            zero = flint_calloc(N, sizeof(ulong));
            for (j = 1; j <= k; j++)
            {
                v = GR_TOWER_FLAT_VAR(F, j);
                mpoly_gen_monomial_sp(gen_packed + N * v, v, bits, F->mctx->minfo);
            }

            {
                cofactor_cache_struct cache;
                cache.len = cache.next = 0;
                for (j = 0; j < COFACTOR_CACHE; j++)
                {
                    fmpz_init(cache.d + j);
                    fmpz_init(cache.q + j);
                }
                _nested_push_terms(num, x, k, zero, gen_packed, N, den, &cache, F);
                for (j = 0; j < COFACTOR_CACHE; j++)
                {
                    fmpz_clear(cache.d + j);
                    fmpz_clear(cache.q + j);
                }
            }
            /* the terms were pushed in decreasing order, distinct */

            flint_free(gen_packed);
            flint_free(zero);
            flint_free(maxlen);
        }
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), den, F->mctx);
        fmpz_mpoly_q_canonicalise(res, F->mctx);

        flint_free(exp);
        fmpz_clear(den);
        return GR_SUCCESS;
    }

    if (k == 0)
        return _base_to_flat(res, x, F);

    {
        /* Horner in the top generator */
        const gr_poly_struct * poly = x;
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        fmpz_mpoly_q_t t;
        fmpz_mpoly_t xk;
        slong i;

        fmpz_mpoly_q_init(t, F->mctx);
        fmpz_mpoly_init(xk, F->mctx);
        fmpz_mpoly_gen(xk, GR_TOWER_FLAT_VAR(F, k), F->mctx);
        fmpz_mpoly_q_zero(res, F->mctx);

        for (i = poly->length - 1; i >= 0 && status == GR_SUCCESS; i--)
        {
            fmpz_mpoly_mul(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), xk, F->mctx);
            status |= _nested_to_flat(t, gr_poly_coeff_srcptr(poly, i, below), k - 1, F);
            fmpz_mpoly_q_add(res, res, t, F->mctx);
        }

        fmpz_mpoly_q_clear(t, F->mctx);
        fmpz_mpoly_clear(xk, F->mctx);
    }

    return status;
}

/* Number of base coefficients of a nested element, counting up to the
   limit. */
slong
_gr_tower_nested_num_coeffs(gr_srcptr x, gr_ctx_t ctx, slong limit)
{
    if (ctx->which_ring == GR_CTX_GR_POLY_QUOTIENT)
    {
        const gr_poly_struct * p = (const gr_poly_struct *) x;
        gr_ctx_struct * cctx = gr_poly_quotient_ctx_base(ctx);
        slong j, s = 0;
        for (j = 0; j < p->length && s <= limit; j++)
            s += _gr_tower_nested_num_coeffs(gr_poly_coeff_srcptr(p, j, cctx), cctx, limit - s);
        return s;
    }
    return 1;
}

/*
    Whether the conversion of the modulus of step k to the flat
    representation is deferred until an element actually has to be
    reduced by it. Over QQ the flat element of a modulus is a primitive
    integer polynomial with an integer leading coefficient, so the
    properties of the ideal recorded at build time are known without it.
    Splitting towers of large polynomials have moduli of binomial size,
    of which elements of moderate degree never need the larger ones
    (symmetric functions of the roots reduce through the last, linear
    steps); keeping only the nested moduli halves the memory of such a
    tower.
*/
static int
_build_ideal_step_deferred(gr_tower_flat_t F, slong k)
{
    gr_tower_struct * T = F->T;
    const gr_poly_struct * m = gr_tower_step_minpoly(T, k);
    gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
    slong i;
    int constant = 1;

    if (!GR_TOWER_BASE_IS_CONSTS(T))
        return 0;

    /* a modulus with constant coefficients is a univariate integer
       polynomial: as compact flat as nested, and wanted for the fast
       univariate reductions */
    for (i = 0; i < m->length && constant; i++)
        constant = (_gr_tower_nested_num_coeffs(gr_poly_coeff_srcptr(m, i, below), below, 1) <= 1);
    if (constant)
        return 0;

    /* (moduli with more base coefficients than the option are converted
       to the flat representation only when needed) */
    {
        slong limit = GR_TOWER_OPTION(T, GR_TOWER_OPT_DEFERRED_IDEAL_SIZE);
        return _gr_tower_nested_num_coeffs(m, gr_tower_field_at(T, k), limit) > limit;
    }
}

static int _materialize_ideal_step(gr_tower_flat_t F, slong k);

/* The flat ideal element of step k, converted from the nested modulus
   on first use when the conversion was deferred. */
const fmpz_mpoly_struct *
_gr_tower_flat_ideal_elem(gr_tower_flat_t F, slong k)
{
    if (F->ideal[k - 1] == NULL)
        GR_MUST_SUCCEED(_materialize_ideal_step(F, k));
    return F->ideal[k - 1];
}

/* Converts the modulus of step k to the flat ideal element k - 1 and
   records its properties (integer leading coefficient, monic univariate
   integer polynomial); or defers the conversion. */
static int
_build_ideal_step(gr_tower_flat_t F, slong k)
{
    F->ideal[k - 1] = NULL;
    F->ideal_univar[k - 1] = NULL;

    if (_build_ideal_step_deferred(F, k))
        return GR_SUCCESS;

    return _materialize_ideal_step(F, k);
}

static int
_materialize_ideal_step(gr_tower_flat_t F, slong k)
{
    gr_tower_struct * T = F->T;
    const gr_poly_struct * m = gr_tower_step_minpoly(T, k);
    gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
    fmpz_mpoly_q_t num, t;
    fmpz_mpoly_t xk;
    slong i;
    int status = GR_SUCCESS;

    F->ideal[k - 1] = flint_malloc(sizeof(fmpz_mpoly_struct));
    fmpz_mpoly_init(F->ideal[k - 1], F->mctx);
    F->ideal_univar[k - 1] = NULL;
    fmpz_mpoly_q_init(num, F->mctx);
    fmpz_mpoly_q_init(t, F->mctx);
    fmpz_mpoly_init(xk, F->mctx);

    fmpz_mpoly_gen(xk, GR_TOWER_FLAT_VAR(F, k), F->mctx);

    if (GR_TOWER_BASE_IS_CONSTS(T))
    {
        /* the polynomial m over F_{k-1} has the shape of an element of
           F_k (of larger length): one term-by-term conversion */
        status |= _nested_to_flat(num, m, k, F);
    }
    else
    {
        /* Horner: num = sum c_i x_k^i as a rational function, then clear
           the denominator */
        for (i = m->length - 1; i >= 0 && status == GR_SUCCESS; i--)
        {
            status |= _nested_to_flat(t, gr_poly_coeff_srcptr(m, i, below), k - 1, F);
            fmpz_mpoly_mul(fmpz_mpoly_q_numref(num), fmpz_mpoly_q_numref(num), xk, F->mctx);
            fmpz_mpoly_q_add(num, num, t, F->mctx);
        }
    }

    /* the denominator is a polynomial in the transcendental variables
       only (coordinates have no algebraic denominators); the ideal
       element is the numerator, made primitive */
    {
        fmpz_t g;
        fmpz_init(g);
        _fmpz_vec_content(g, fmpz_mpoly_q_numref(num)->coeffs, fmpz_mpoly_q_numref(num)->length);
        if (!fmpz_is_one(g) && !fmpz_is_zero(g))
            fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_q_numref(num), fmpz_mpoly_q_numref(num), g, F->mctx);
        fmpz_clear(g);
    }

    fmpz_mpoly_swap(F->ideal[k - 1], fmpz_mpoly_q_numref(num), F->mctx);

    fmpz_mpoly_q_clear(num, F->mctx);
    fmpz_mpoly_q_clear(t, F->mctx);
    fmpz_mpoly_clear(xk, F->mctx);

    /* with integer leading coefficients, plain division applies */
    {
        slong v = GR_TOWER_FLAT_VAR(F, k);
        slong d = gr_tower_step_degree(T, k);
        fmpz_mpoly_t lc;
        ulong dd = d;
        fmpz_mpoly_init(lc, F->mctx);
        /* the leading coefficient in the variable of the generator */
        if (fmpz_mpoly_degree_si(F->ideal[k - 1], v, F->mctx) == d)
            fmpz_mpoly_get_coeff_vars_ui(lc, F->ideal[k - 1], &v, &dd, 1, F->mctx);
        if (fmpz_mpoly_is_zero(lc, F->mctx) || !fmpz_mpoly_is_fmpz(lc, F->mctx))
            F->ideal_lc_const = 0;
        else if (fmpz_mpoly_is_one(lc, F->mctx))
        {
            /* a monic univariate integer modulus: reductions of elements
               involving that variable only can use fast univariate
               division */
            int * used = flint_malloc(sizeof(int) * F->cap);
            slong nused = 0;

            fmpz_mpoly_used_vars(used, F->ideal[k - 1], F->mctx);
            for (i = 0; i < F->cap; i++)
                nused += (used[i] != 0);
            flint_free(used);

            if (nused == 1)
            {
                fmpz_poly_t P;
                fmpz_poly_init(P);
                fmpz_mpoly_get_fmpz_poly(P, F->ideal[k - 1], v, F->mctx);
                F->ideal_univar[k - 1] = flint_malloc(sizeof(fmpz_poly_struct));
                fmpz_poly_init(F->ideal_univar[k - 1]);
                fmpz_poly_swap(F->ideal_univar[k - 1], P);
                fmpz_poly_clear(P);
            }
        }
        fmpz_mpoly_clear(lc, F->mctx);
    }

    return status;
}

/*
    Builds the flat ideal from the moduli. When steps were only appended
    since the last build (the moduli of the existing steps unchanged, as
    recorded by the moduli version of the tower), only the new steps are
    converted.
*/
static int
_build_ideal(gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong k, n = T->length, first;
    int status = GR_SUCCESS;

    if (F->ideal_len > 0 && F->ideal_len <= n && F->ideal_moduli_version == T->moduli_version)
    {
        first = F->ideal_len + 1;
        F->ideal = flint_realloc(F->ideal, sizeof(fmpz_mpoly_struct *) * FLINT_MAX(n, 1));
        F->ideal_univar = flint_realloc(F->ideal_univar, sizeof(fmpz_poly_struct *) * FLINT_MAX(n, 1));
    }
    else
    {
        _flat_clear_ideal(F);
        F->ideal = flint_malloc(sizeof(fmpz_mpoly_struct *) * FLINT_MAX(n, 1));
        F->ideal_univar = flint_calloc(FLINT_MAX(n, 1), sizeof(fmpz_poly_struct *));
        F->ideal_lc_const = 1;
        first = 1;
    }

    for (k = first; k <= n; k++)
    {
        status |= _build_ideal_step(F, k);
        F->ideal_len = k;
    }

    F->ideal_version = T->version;
    F->ideal_moduli_version = T->moduli_version;

    return status;
}

/*
    Reduces the polynomial A if it involves at most one variable and that
    variable has a monic univariate integer modulus (or is transcendental),
    using univariate division. Returns 1 if A is reduced on return.
*/
static int
_univar_reduce(fmpz_mpoly_t A, gr_tower_flat_t F)
{
    int * used;
    slong i, v = -1, nused = 0, k;

    if (fmpz_mpoly_is_fmpz(A, F->mctx))
        return 1;

    used = flint_malloc(sizeof(int) * F->cap);
    fmpz_mpoly_used_vars(used, A, F->mctx);
    for (i = 0; i < F->cap; i++)
    {
        if (used[i])
        {
            nused++;
            v = i;
        }
    }
    flint_free(used);

    if (nused != 1)
        return 0;

    for (k = 1; k <= F->ideal_len; k++)
    {
        if (GR_TOWER_FLAT_VAR(F, k) == v)
        {
            fmpz_poly_t P;

            if (F->ideal_univar[k - 1] == NULL)
                return 0;

            fmpz_poly_init(P);
            fmpz_mpoly_get_fmpz_poly(P, A, v, F->mctx);
            if (fmpz_poly_degree(P) >= fmpz_poly_degree(F->ideal_univar[k - 1]))
            {
                fmpz_poly_rem(P, P, F->ideal_univar[k - 1]);
                fmpz_mpoly_set_fmpz_poly(A, P, v, F->mctx);
            }
            fmpz_poly_clear(P);
            return 1;
        }
    }

    /* a transcendental variable: nothing to reduce */
    return 1;
}

/*
    Replaces the polynomial context (keeping the old one alive) when the
    tower has grown beyond the capacity or its definition order has
    changed; the new context has the given capacity. Returns 1 if the
    context was replaced.
*/
static int
_gr_tower_flat_relayout(gr_tower_flat_t F, slong new_cap)
{
    gr_tower_struct * T = F->T;

    _flat_clear_ideal(F);

    F->old_mctx = flint_realloc(F->old_mctx, (F->num_old_mctx + 1) * sizeof(fmpz_mpoly_ctx_struct *));
    F->old_var_gid = flint_realloc(F->old_var_gid, (F->num_old_mctx + 1) * sizeof(slong *));
    F->old_mctx[F->num_old_mctx] = F->mctx;
    F->old_var_gid[F->num_old_mctx] = F->var_gid;
    F->num_old_mctx++;

    F->mctx = flint_malloc(sizeof(fmpz_mpoly_ctx_struct));
    fmpz_mpoly_ctx_init(F->mctx, new_cap, ORD_LEX);
    F->cap = new_cap;
    F->var_gid = _make_var_gid(T, new_cap);
    F->layout_version = T->structure_version;
    return 1;
}

int
_gr_tower_flat_grow(gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;

    if (T->num_gens > F->cap)
        return _gr_tower_flat_relayout(F, FLINT_MAX(T->num_gens, 2 * F->cap));

    if (F->layout_version != T->structure_version)
        return _gr_tower_flat_relayout(F, F->cap);

    /* new generators (appended) get their variables (all of them are
       set whenever the last one is not: nothing to do in the common
       case) */
    if (T->num_gens > 0 && F->var_gid[F->cap - T->num_gens] != T->gens[T->num_gens - 1].gid)
    {
        slong d;
        for (d = 0; d < T->num_gens; d++)
            F->var_gid[F->cap - 1 - d] = T->gens[d].gid;
    }

    return 0;
}

int
gr_tower_flat_ensure(gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    int grown = _gr_tower_flat_grow(F);

    if (F->ideal_len != T->length || F->ideal_version != T->version)
        GR_MUST_SUCCEED(_build_ideal(F));

    return grown;
}

/*
    Reduces the polynomial P modulo the triangular set by pseudo-division,
    from the top generator down: whenever deg_{x_k}(P) = e >= d_k, P is
    replaced by its pseudo-remainder by m_k, lc_k^{e - d_k + 1} P mod m_k,
    where lc_k (the leading coefficient of m_k, a polynomial in the
    transcendental variables only) is accumulated in mult. Since lc_k
    involves no algebraic variable, multiplications by it never unreduce
    anything, so a single top-down pass gives exponents < d_k for every
    algebraic generator.
*/
static void
_flat_reduce_poly(fmpz_mpoly_t P, fmpz_mpoly_t mult, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong k;

    fmpz_mpoly_one(mult, F->mctx);

    for (k = T->length; k >= 1; k--)
    {
        slong v = GR_TOWER_FLAT_VAR(F, k);
        slong d = gr_tower_step_degree(T, k);
        slong e = fmpz_mpoly_degree_si(P, v, F->mctx);
        fmpz_mpoly_univar_t Pu, mu;
        fmpz_mpoly_t lc;

        if (e < d)
            continue;

        fmpz_mpoly_univar_init(Pu, F->mctx);
        fmpz_mpoly_univar_init(mu, F->mctx);
        fmpz_mpoly_init(lc, F->mctx);

        fmpz_mpoly_to_univar(mu, _gr_tower_flat_ideal_elem(F, k), v, F->mctx);
        fmpz_mpoly_to_univar(Pu, P, v, F->mctx);
        fmpz_mpoly_univar_pseudo_rem(Pu, Pu, mu, F->mctx);
        fmpz_mpoly_from_univar(P, Pu, v, F->mctx);

        fmpz_mpoly_univar_get_term_coeff(lc, mu, 0, F->mctx);
        if (!fmpz_mpoly_is_one(lc, F->mctx))
        {
            fmpz_mpoly_pow_ui(lc, lc, e - d + 1, F->mctx);
            fmpz_mpoly_mul(mult, mult, lc, F->mctx);
        }

        fmpz_mpoly_univar_clear(Pu, F->mctx);
        fmpz_mpoly_univar_clear(mu, F->mctx);
        fmpz_mpoly_clear(lc, F->mctx);
    }
}

/* scale P = (multiple of the ideal) + R, R reduced, by one division per
   step from the top down (see gr_tower_flat_reduce) */
static void
_flat_reduce_sequential(fmpz_mpoly_t P, fmpz_t scale, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong k;

    fmpz_one(scale);

    for (k = T->length; k >= 1; k--)
    {
        slong v = GR_TOWER_FLAT_VAR(F, k);
        slong d = gr_tower_step_degree(T, k);

        if (fmpz_mpoly_degree_si(P, v, F->mctx) >= d)
        {
            const fmpz_mpoly_struct * m = _gr_tower_flat_ideal_elem(F, k);
            fmpz_mpoly_struct * ms = (fmpz_mpoly_struct *) m;
            fmpz_mpoly_t Q, R;
            fmpz_mpoly_struct * Qs = Q;
            fmpz_t s;

            fmpz_mpoly_init(Q, F->mctx);
            fmpz_mpoly_init(R, F->mctx);
            fmpz_init(s);
            fmpz_mpoly_quasidivrem_ideal(s, &Qs, R, P, &ms, 1, F->mctx);
            fmpz_mpoly_swap(P, R, F->mctx);
            fmpz_mul(scale, scale, s);
            fmpz_mpoly_clear(Q, F->mctx);
            fmpz_mpoly_clear(R, F->mctx);
            fmpz_clear(s);
        }
    }
}

/* whether the degree of P in every algebraic variable is below the
   degree of its modulus */
static int
_flat_is_reduced(const fmpz_mpoly_t P, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong * degs;
    slong k;
    int ok = 1;

    if (P->length <= 1 && fmpz_mpoly_is_fmpz(P, F->mctx))
        return 1;

    degs = flint_malloc(sizeof(slong) * F->cap);
    fmpz_mpoly_degrees_si(degs, P, F->mctx);
    for (k = 1; k <= T->length && ok; k++)
        if (degs[GR_TOWER_FLAT_VAR(F, k)] >= gr_tower_step_degree(T, k))
            ok = 0;
    flint_free(degs);
    return ok;
}

int
gr_tower_flat_reduce(fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    fmpz_mpoly_t mult;
    int changed = 0;

    gr_tower_flat_ensure(F);

    if (F->ideal_len == 0)
        return GR_SUCCESS;

    /* nothing to reduce (all degrees below the moduli's): the fraction,
       being canonical already, stays as it is -- the gcd of the
       canonicalisation would only be recomputed */
    if (_flat_is_reduced(fmpz_mpoly_q_numref(x), F) && _flat_is_reduced(fmpz_mpoly_q_denref(x), F))
        return GR_SUCCESS;

    if (_univar_reduce(fmpz_mpoly_q_numref(x), F) && _univar_reduce(fmpz_mpoly_q_denref(x), F))
    {
        fmpz_mpoly_q_canonicalise(x, F->mctx);
        return GR_SUCCESS;
    }

    if (F->ideal_lc_const)
    {
        /* integer leading coefficients: the triangular set is a Groebner
           basis with pure-power leading monomials; a single division pass */
        slong i, n = F->ideal_len;
        fmpz_mpoly_struct ** Q;
        fmpz_mpoly_t R;
        fmpq_t scale;

        for (i = 0; i < n; i++)
            if (F->ideal[i] == NULL)
                break;

        if (i < n)
        {
            /* some conversions deferred: one division per step from the
               top down, converting the moduli which are actually needed
               (a reduction by the modulus of step k raises only degrees
               in the variables of steps below k) */
            fmpq_init(scale);
            _flat_reduce_sequential(fmpz_mpoly_q_numref(x), fmpq_denref(scale), F);
            if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), F->mctx))
            {
                fmpz_mpoly_q_zero(x, F->mctx);
                fmpq_clear(scale);
                return GR_SUCCESS;
            }
            _flat_reduce_sequential(fmpz_mpoly_q_denref(x), fmpq_numref(scale), F);
            fmpz_mpoly_q_canonicalise(x, F->mctx);
            if (!fmpq_is_one(scale))
            {
                fmpq_canonicalise(scale);
                fmpz_mpoly_q_mul_fmpq(x, x, scale, F->mctx);
            }
            fmpq_clear(scale);
            return GR_SUCCESS;
        }

        Q = flint_malloc(sizeof(fmpz_mpoly_struct *) * n);
        for (i = 0; i < n; i++)
        {
            Q[i] = flint_malloc(sizeof(fmpz_mpoly_struct));
            fmpz_mpoly_init(Q[i], F->mctx);
        }
        fmpz_mpoly_init(R, F->mctx);
        fmpq_init(scale);

        fmpz_mpoly_quasidivrem_ideal(fmpq_denref(scale), Q, R, fmpz_mpoly_q_numref(x), F->ideal, n, F->mctx);
        fmpz_mpoly_swap(R, fmpz_mpoly_q_numref(x), F->mctx);

        if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), F->mctx))
        {
            fmpz_mpoly_q_zero(x, F->mctx);
            fmpq_one(scale);
        }
        else if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx))
        {
            fmpz_one(fmpq_numref(scale));
        }
        else
        {
            fmpz_mpoly_quasidivrem_ideal(fmpq_numref(scale), Q, R, fmpz_mpoly_q_denref(x), F->ideal, n, F->mctx);
            fmpz_mpoly_swap(R, fmpz_mpoly_q_denref(x), F->mctx);
        }

        fmpz_mpoly_q_canonicalise(x, F->mctx);

        if (!fmpq_is_one(scale))
        {
            fmpq_canonicalise(scale);
            fmpz_mpoly_q_mul_fmpq(x, x, scale, F->mctx);
        }

        for (i = 0; i < n; i++)
        {
            fmpz_mpoly_clear(Q[i], F->mctx);
            flint_free(Q[i]);
        }
        flint_free(Q);
        fmpz_mpoly_clear(R, F->mctx);
        fmpq_clear(scale);
        return GR_SUCCESS;
    }

    fmpz_mpoly_init(mult, F->mctx);

    _flat_reduce_poly(fmpz_mpoly_q_numref(x), mult, F);

    /* a zero (the common case of a zero test): the denominator, often
       the larger part (the product of those of the operands), is not
       needed */
    if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), F->mctx))
    {
        fmpz_mpoly_clear(mult, F->mctx);
        fmpz_mpoly_q_zero(x, F->mctx);
        return GR_SUCCESS;
    }

    if (!fmpz_mpoly_is_one(mult, F->mctx))
    {
        fmpz_mpoly_mul(fmpz_mpoly_q_denref(x), fmpz_mpoly_q_denref(x), mult, F->mctx);
        changed = 1;
    }

    _flat_reduce_poly(fmpz_mpoly_q_denref(x), mult, F);
    if (!fmpz_mpoly_is_one(mult, F->mctx))
    {
        fmpz_mpoly_mul(fmpz_mpoly_q_numref(x), fmpz_mpoly_q_numref(x), mult, F->mctx);
        changed = 1;
    }

    fmpz_mpoly_clear(mult, F->mctx);

    /* the multipliers are polynomials in the transcendental variables, so
       the fraction is still reduced with respect to the algebraic ones;
       canonicalise (gcd) in any case since reductions change the parts */
    fmpz_mpoly_q_canonicalise(x, F->mctx);
    (void) changed;

    return GR_SUCCESS;
}

int
gr_tower_flat_set_nested_at(fmpz_mpoly_q_t res, gr_srcptr x, slong k, gr_tower_flat_t F)
{
    gr_tower_flat_ensure(F);
    return _nested_to_flat(res, x, k, F);
}

/*
    Polynomial in the flat variables (reduced, involving algebraic
    generators <= k) -> nested element of F_k: the coefficients with
    respect to the generator a_k are converted recursively (sparse in
    the flat variables); at level 0 the polynomial is in the
    transcendental variables only.
*/
static int
_flat_poly_to_nested(gr_ptr res, const fmpz_mpoly_t f, slong k, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    int status = GR_SUCCESS;

    if (k == 0)
        return _flat_to_base(res, f, F);

    {
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        gr_poly_struct * poly = res;
        slong v = GR_TOWER_FLAT_VAR(F, k), d = gr_tower_step_degree(T, k), i, n;
        fmpz_mpoly_univar_t u;
        fmpz_mpoly_t c;

        fmpz_mpoly_univar_init(u, F->mctx);
        fmpz_mpoly_init(c, F->mctx);
        fmpz_mpoly_to_univar(u, f, v, F->mctx);
        n = fmpz_mpoly_univar_length(u, F->mctx);

        gr_poly_fit_length(poly, d, below);
        _gr_poly_set_length(poly, d, below);
        status |= _gr_vec_zero(poly->coeffs, d, below);

        for (i = 0; i < n && status == GR_SUCCESS; i++)
        {
            slong e = fmpz_mpoly_univar_get_term_exp_si(u, i, F->mctx);
            if (e >= d)
            {
                status = GR_DOMAIN;   /* not reduced */
                break;
            }
            fmpz_mpoly_univar_swap_term_coeff(c, u, i, F->mctx);
            status |= _flat_poly_to_nested(gr_poly_coeff_ptr(poly, e, below), c, k - 1, F);
        }

        _gr_poly_normalise(poly, below);

        fmpz_mpoly_univar_clear(u, F->mctx);
        fmpz_mpoly_clear(c, F->mctx);
    }

    return status;
}

int
gr_tower_flat_poly_get_nested_at(gr_ptr res, const fmpz_mpoly_t f_in, slong k, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    gr_ctx_struct * F0 = T->base;
    fmpz_mpoly_q_t q;
    int status = GR_SUCCESS;

    gr_tower_flat_ensure(F);

    fmpz_mpoly_q_init(q, F->mctx);
    fmpz_mpoly_set(fmpz_mpoly_q_numref(q), f_in, F->mctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(q), F->mctx);
    status |= gr_tower_flat_reduce(q, F);

    /* generators above level k must not occur */
    if (gr_tower_flat_max_step(fmpz_mpoly_q_numref(q), F) > k)
        status = GR_DOMAIN;

    if (status == GR_SUCCESS)
        status = _flat_poly_to_nested(res, fmpz_mpoly_q_numref(q), k, F);

    /* the reduction may have introduced a denominator (an integer, or a
       polynomial in the transcendental variables) */
    if (status == GR_SUCCESS && fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(q), F->mctx) && !fmpz_mpoly_is_one(fmpz_mpoly_q_denref(q), F->mctx))
    {
        fmpz_t c;
        fmpz_init(c);
        fmpz_mpoly_get_fmpz(c, fmpz_mpoly_q_denref(q), F->mctx);
        status |= gr_div_fmpz(res, res, c, gr_tower_field_at(T, k));
        fmpz_clear(c);
    }
    else if (status == GR_SUCCESS && !fmpz_mpoly_is_one(fmpz_mpoly_q_denref(q), F->mctx))
    {
        gr_ctx_struct * Fk = gr_tower_field_at(T, k);
        gr_ptr d;
        GR_TMP_INIT(d, F0);
        status |= _flat_to_base(d, fmpz_mpoly_q_denref(q), F);
        if (status == GR_SUCCESS)
        {
            gr_ptr dk;
            GR_TMP_INIT(dk, Fk);
            status |= gr_tower_promote(dk, d, 0, k, T);
            status |= gr_div(res, res, dk, Fk);
            GR_TMP_CLEAR(dk, Fk);
        }
        GR_TMP_CLEAR(d, F0);
    }

    fmpz_mpoly_q_clear(q, F->mctx);
    return status;
}

int
gr_tower_flat_get_nested_at(gr_ptr res, const fmpz_mpoly_q_t x, slong k, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    gr_ctx_struct * Fk = gr_tower_field_at(T, k);
    gr_ptr d;
    int status;

    status = gr_tower_flat_poly_get_nested_at(res, fmpz_mpoly_q_numref(x), k, F);

    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(x), F->mctx))
        return status;

    /* an integer denominator: scalar division (a field inversion would
       be far more expensive) */
    if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx))
    {
        fmpz_t c;
        fmpz_init(c);
        fmpz_mpoly_get_fmpz(c, fmpz_mpoly_q_denref(x), F->mctx);
        status |= gr_div_fmpz(res, res, c, Fk);
        fmpz_clear(c);
        return status;
    }

    /* a denominator which is a monomial in the algebraic generators is
       cleared through the moduli (cheaply); otherwise the division is
       an inversion in the nested field */
    if (gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(x), F))
    {
        fmpz_mpoly_q_t y;
        fmpz_mpoly_q_init(y, F->mctx);
        fmpz_mpoly_q_set(y, x, F->mctx);
        if (gr_tower_flat_rationalize(y, F) == GR_SUCCESS)
        {
            status = gr_tower_flat_get_nested_at(res, y, k, F);
            fmpz_mpoly_q_clear(y, F->mctx);
            return status;
        }
        fmpz_mpoly_q_clear(y, F->mctx);
    }

    GR_TMP_INIT(d, Fk);
    status |= gr_tower_flat_poly_get_nested_at(d, fmpz_mpoly_q_denref(x), k, F);
    if (status == GR_SUCCESS)
        status = gr_tower_div_at(res, res, d, k, T);
    GR_TMP_CLEAR(d, Fk);
    return status;
}

/* The largest index k of an algebraic generator a_k occurring in p
   (0 if none). */
slong
gr_tower_flat_max_step(const fmpz_mpoly_t p, gr_tower_flat_t F)
{
    int * used;
    slong k;

    if (fmpz_mpoly_is_fmpz(p, F->mctx))
        return 0;
    used = flint_malloc(sizeof(int) * F->cap);
    fmpz_mpoly_used_vars(used, p, F->mctx);
    for (k = F->T->length; k >= 1; k--)
        if (used[GR_TOWER_FLAT_VAR(F, k)])
            break;
    flint_free(used);
    return k;
}

/* Whether p involves an algebraic variable of F. */
int
gr_tower_flat_has_alg_var(const fmpz_mpoly_t p, gr_tower_flat_t F)
{
    return gr_tower_flat_max_step(p, F) > 0;
}

int gr_tower_flat_rationalize(fmpz_mpoly_q_t x, gr_tower_flat_t F);

/*
    The inverse of the algebraic generator a_k as a flat element with
    denominator free of algebraic variables: from the ideal element
    c_0 + c_1 a + ... + c_d a^d = 0, 1/a = -(c_1 + c_2 a + ... + c_d
    a^(d-1)) / c_0, where 1/c_0 is rationalized recursively (c_0 involves
    only lower generators).
*/
static int
_flat_gen_inverse(fmpz_mpoly_q_t res, slong k, gr_tower_flat_t F)
{
    slong v = GR_TOWER_FLAT_VAR(F, k), i, n;
    fmpz_mpoly_univar_t u;
    fmpz_mpoly_t c, num, xpow;
    fmpz_mpoly_t c0;
    int status = GR_SUCCESS;

    gr_tower_flat_ensure(F);

    fmpz_mpoly_univar_init(u, F->mctx);
    fmpz_mpoly_init(c, F->mctx);
    fmpz_mpoly_init(c0, F->mctx);
    fmpz_mpoly_init(num, F->mctx);
    fmpz_mpoly_init(xpow, F->mctx);

    fmpz_mpoly_to_univar(u, _gr_tower_flat_ideal_elem(F, k), v, F->mctx);
    n = fmpz_mpoly_univar_length(u, F->mctx);

    for (i = 0; i < n; i++)
    {
        slong e = fmpz_mpoly_univar_get_term_exp_si(u, i, F->mctx);
        fmpz_mpoly_univar_get_term_coeff(c, u, i, F->mctx);
        if (e == 0)
        {
            fmpz_mpoly_set(c0, c, F->mctx);
        }
        else
        {
            fmpz_mpoly_gen(xpow, v, F->mctx);
            fmpz_mpoly_pow_ui(xpow, xpow, e - 1, F->mctx);
            fmpz_mpoly_mul(xpow, xpow, c, F->mctx);
            fmpz_mpoly_add(num, num, xpow, F->mctx);
        }
    }

    if (fmpz_mpoly_is_zero(c0, F->mctx))
    {
        status = GR_DOMAIN;   /* the generator is zero (not a field) */
    }
    else
    {
        fmpz_mpoly_neg(num, num, F->mctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_numref(res), num, F->mctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_denref(res), c0, F->mctx);
        fmpz_mpoly_q_canonicalise(res, F->mctx);
        status = gr_tower_flat_rationalize(res, F);
        if (status == GR_SUCCESS)
            status = gr_tower_flat_reduce(res, F);
    }

    fmpz_mpoly_univar_clear(u, F->mctx);
    fmpz_mpoly_clear(c, F->mctx);
    fmpz_mpoly_clear(c0, F->mctx);
    fmpz_mpoly_clear(num, F->mctx);
    fmpz_mpoly_clear(xpow, F->mctx);
    return status;
}

/*
    Clears the algebraic variables from the denominator of x when the
    denominator is a monomial in them (the inverse of each generator is
    expressed through its modulus, recursively), returning GR_UNABLE
    otherwise (a general denominator needs an inversion in the nested
    field). The result is reduced.
*/
/*
    Clears the highest algebraic variable v (of level k) from the
    denominator of x, a polynomial in v of degree < d over the lower
    field: with M the matrix of multiplication by den modulo the modulus
    of v (a d x d matrix over the lower field), den * c = det(M) for
    c = adj(M) e_0, so that 1/den = (sum c_i v^i) / det(M), where the
    norm det(M) is free of v. The lower algebraic variables are cleared
    from the new denominator recursively.

    When den is a polynomial in w = v^g (g | d) and the multiples of den
    by the powers of w stay polynomials in w (as for a binomial modulus
    v^d - c: den = 2 sqrt(2) + i with sqrt(2) = 2^(1/16)^8, say), the
    same is done with w, of degree d/g: then det(M) is the norm from the
    field generated by w.
*/
static int
_flat_rationalize_top(fmpz_mpoly_q_t x, slong k, gr_tower_flat_t F)
{
    slong v = GR_TOWER_FLAT_VAR(F, k), d, i, j, g;
    gr_ctx_t Q;
    gr_mat_t M, adj;
    gr_ptr det;
    fmpz_mpoly_q_t p, c, t;
    fmpz_mpoly_t xv;
    fmpz_mpoly_univar_t u;
    int status = GR_SUCCESS;

    d = gr_tower_step_degree(F->T, k);

    /* g = gcd of d and the exponents of v in den */
    g = d;
    {
        fmpz_mpoly_univar_t w;
        fmpz_mpoly_univar_init(w, F->mctx);
        fmpz_mpoly_to_univar(w, fmpz_mpoly_q_denref(x), v, F->mctx);
        for (i = 0; i < fmpz_mpoly_univar_length(w, F->mctx) && g > 1; i++)
            g = n_gcd(g, fmpz_mpoly_univar_get_term_exp_si(w, i, F->mctx));
        fmpz_mpoly_univar_clear(w, F->mctx);
    }
    d /= g;

    if (d > 8)
        return GR_UNABLE;

    /* (the norm of a large denominator is enormous: not attempted) */
    if (fmpz_mpoly_q_denref(x)->length * d > GR_TOWER_OPTION(F->T, GR_TOWER_OPT_RATIONALIZE_LIMIT))
        return GR_UNABLE;

    /* (a context with the same variables: elements are interchangeable) */
    gr_ctx_init_fmpz_mpoly_q(Q, F->cap, ORD_LEX);
    gr_mat_init(M, d, d, Q);
    gr_mat_init(adj, d, d, Q);
    det = gr_heap_init(Q);
    fmpz_mpoly_q_init(p, F->mctx);
    fmpz_mpoly_q_init(c, F->mctx);
    fmpz_mpoly_q_init(t, F->mctx);
    fmpz_mpoly_init(xv, F->mctx);
    fmpz_mpoly_univar_init(u, F->mctx);

    /* column j of M: den * w^j reduced, split by powers of w = v^g */
    fmpz_mpoly_gen(xv, v, F->mctx);
    fmpz_mpoly_pow_ui(xv, xv, g, F->mctx);
    fmpz_mpoly_set(fmpz_mpoly_q_numref(p), fmpz_mpoly_q_denref(x), F->mctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(p), F->mctx);
    for (j = 0; j < d && status == GR_SUCCESS; j++)
    {
        slong n;

        if (j > 0)
            fmpz_mpoly_mul(fmpz_mpoly_q_numref(p), fmpz_mpoly_q_numref(p), xv, F->mctx);
        status |= gr_tower_flat_reduce(p, F);
        if (status != GR_SUCCESS)
            break;
        if (fmpz_mpoly_degree_si(fmpz_mpoly_q_denref(p), v, F->mctx) > 0)
        {
            status = GR_UNABLE;   /* (the reduction does not put v in denominators) */
            break;
        }

        fmpz_mpoly_to_univar(u, fmpz_mpoly_q_numref(p), v, F->mctx);
        n = fmpz_mpoly_univar_length(u, F->mctx);
        for (i = 0; i < n; i++)
        {
            slong e = fmpz_mpoly_univar_get_term_exp_si(u, i, F->mctx);
            fmpz_mpoly_q_struct * entry;
            if (e % g != 0 || e / g >= d)
            {
                status = GR_UNABLE;
                break;
            }
            e /= g;
            entry = gr_mat_entry_ptr(M, e, j, Q);
            fmpz_mpoly_univar_get_term_coeff(fmpz_mpoly_q_numref(entry), u, i, F->mctx);
            fmpz_mpoly_set(fmpz_mpoly_q_denref(entry), fmpz_mpoly_q_denref(p), F->mctx);
            fmpz_mpoly_q_canonicalise(entry, F->mctx);
        }
    }

    if (status == GR_SUCCESS)
        status = gr_mat_adjugate(adj, det, M, Q);

    /* (the determinant is computed in the polynomial ring: reduced
       modulo the ideal, it is the norm) */
    if (status == GR_SUCCESS)
        status = gr_tower_flat_reduce(det, F);

    if (status == GR_SUCCESS && fmpz_mpoly_q_is_zero(det, F->mctx))
        status = GR_DOMAIN;   /* den is a zero divisor: not a field */

    if (status == GR_SUCCESS)
    {
        /* x = num / den = num * (sum_i adj[i][0] v^i) / det */
        fmpz_mpoly_q_zero(c, F->mctx);
        for (i = d - 1; i >= 0; i--)
        {
            fmpz_mpoly_mul(fmpz_mpoly_q_numref(c), fmpz_mpoly_q_numref(c), xv, F->mctx);
            fmpz_mpoly_q_add(c, c, gr_mat_entry_ptr(adj, i, 0, Q), F->mctx);
        }
        fmpz_mpoly_one(fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_mpoly_q_mul(x, x, c, F->mctx);
        fmpz_mpoly_q_div(x, x, det, F->mctx);
        status = gr_tower_flat_reduce(x, F);
    }

    gr_mat_clear(M, Q);
    gr_mat_clear(adj, Q);
    gr_heap_clear(det, Q);
    gr_ctx_clear(Q);
    fmpz_mpoly_q_clear(p, F->mctx);
    fmpz_mpoly_q_clear(c, F->mctx);
    fmpz_mpoly_q_clear(t, F->mctx);
    fmpz_mpoly_clear(xv, F->mctx);
    fmpz_mpoly_univar_clear(u, F->mctx);
    return status;
}

int
gr_tower_flat_rationalize(fmpz_mpoly_q_t x_out, gr_tower_flat_t F)
{
    fmpz_mpoly_q_t x;
    fmpz_mpoly_struct * den;
    ulong * exp;
    slong k, e, i;
    fmpz_mpoly_t rest;
    fmpz_mpoly_q_t inv;
    fmpz_t cf;
    int status = GR_SUCCESS;

    gr_tower_flat_ensure(F);

    if (!gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(x_out), F))
        return GR_SUCCESS;

    if (fmpz_mpoly_q_denref(x_out)->length != 1)
    {
        /* a polynomial denominator: cleared variable by variable, from
           the highest algebraic one down */
        fmpz_mpoly_q_t y;
        int * used;
        slong iter;

        fmpz_mpoly_q_init(y, F->mctx);
        fmpz_mpoly_q_set(y, x_out, F->mctx);
        used = flint_malloc(sizeof(int) * F->cap);

        for (iter = 0; iter <= F->T->length && status == GR_SUCCESS; iter++)
        {
            slong top = 0;
            fmpz_mpoly_used_vars(used, fmpz_mpoly_q_denref(y), F->mctx);
            for (k = 1; k <= F->T->length; k++)
                if (used[GR_TOWER_FLAT_VAR(F, k)])
                    top = k;
            if (top == 0)
                break;
            status = _flat_rationalize_top(y, top, F);
        }

        if (status == GR_SUCCESS)
            fmpz_mpoly_q_swap(x_out, y, F->mctx);
        fmpz_mpoly_q_clear(y, F->mctx);
        flint_free(used);
        return status;
    }

    /* (x_out is left untouched unless the rationalization succeeds) */
    fmpz_mpoly_q_init(x, F->mctx);
    fmpz_mpoly_q_set(x, x_out, F->mctx);
    den = fmpz_mpoly_q_denref(x);

    exp = flint_malloc(sizeof(ulong) * F->cap);
    fmpz_mpoly_init(rest, F->mctx);
    fmpz_mpoly_q_init(inv, F->mctx);
    fmpz_init(cf);

    /* the denominator without its algebraic variables */
    fmpz_mpoly_get_term_exp_ui(exp, den, 0, F->mctx);
    fmpz_mpoly_get_term_coeff_fmpz(cf, den, 0, F->mctx);
    {
        ulong * exp2 = flint_malloc(sizeof(ulong) * F->cap);
        for (i = 0; i < F->cap; i++)
            exp2[i] = exp[i];
        for (k = 1; k <= F->T->length; k++)
            exp2[GR_TOWER_FLAT_VAR(F, k)] = 0;
        fmpz_mpoly_push_term_fmpz_ui(rest, cf, exp2, F->mctx);
        flint_free(exp2);
    }
    fmpz_mpoly_swap(den, rest, F->mctx);

    for (k = 1; k <= F->T->length && status == GR_SUCCESS; k++)
    {
        e = exp[GR_TOWER_FLAT_VAR(F, k)];
        if (e == 0)
            continue;
        status = _flat_gen_inverse(inv, k, F);
        for (i = 0; i < e && status == GR_SUCCESS; i++)
            fmpz_mpoly_q_mul(x, x, inv, F->mctx);
    }

    if (status == GR_SUCCESS)
        status = gr_tower_flat_reduce(x, F);
    if (status == GR_SUCCESS)
        fmpz_mpoly_q_swap(x_out, x, F->mctx);

    flint_free(exp);
    fmpz_mpoly_clear(rest, F->mctx);
    fmpz_mpoly_q_clear(inv, F->mctx);
    fmpz_mpoly_q_clear(x, F->mctx);
    fmpz_clear(cf);
    return status;
}

/*
    Numerical evaluation. Only the variables actually occurring are
    evaluated, so that the arguments of transcendental generators
    (which involve only earlier generators) can be evaluated recursively.
*/
int
gr_tower_flat_get_acb(acb_t res, const fmpz_mpoly_q_t x, slong prec, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong k, tot = F->cap;
    acb_ptr vals;
    int * used;
    int status = GR_SUCCESS;

    vals = _acb_vec_init(tot);
    used = flint_malloc(sizeof(int) * tot);
    fmpz_mpoly_q_used_vars(used, x, F->mctx);

    for (k = 1; k <= T->length && status == GR_SUCCESS; k++)
        if (used[GR_TOWER_FLAT_VAR(F, k)])
            status |= gr_tower_step_get_acb(vals + GR_TOWER_FLAT_VAR(F, k), T, k, prec);

    for (k = 1; k <= T->num_trans && status == GR_SUCCESS; k++)
        if (used[GR_TOWER_FLAT_TVAR(F, k)])
            status |= gr_tower_trans_get_acb(vals + GR_TOWER_FLAT_TVAR(F, k), T, k, prec);

    if (status == GR_SUCCESS)
        fmpz_mpoly_q_evaluate_acb(res, x, vals, prec, F->mctx);

    _acb_vec_clear(vals, tot);
    flint_free(used);
    return status;
}

slong
gr_tower_flat_level(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    return gr_tower_flat_def_order(x, F) + 1;
}

/* Highest algebraic generator occurring in x (0 if none). */
slong
gr_tower_flat_alg_level(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    return FLINT_MAX(gr_tower_flat_max_step(fmpz_mpoly_q_numref(x), F),
                     gr_tower_flat_max_step(fmpz_mpoly_q_denref(x), F));
}

/* the most significant variable (least index, lex order) occurring in P,
   or nvars if P is constant: it occurs in the leading term */
static slong
_lead_var(const fmpz_mpoly_t P, gr_tower_flat_t F)
{
    slong nvars = F->mctx->minfo->nvars, v, res = nvars;
    ulong e_small[16];
    ulong * e;

    if (P->length == 0)
        return nvars;

    /* packed exponents, lex order: the variable var is the field
       nvars - 1 - var, fields packed FLINT_BITS / bits per word (none
       straddling two words); the most significant variable occurring is
       the highest nonzero field of the leading monomial */
    if (P->bits <= FLINT_BITS && !F->mctx->minfo->deg && !F->mctx->minfo->rev)
    {
        slong N = mpoly_words_per_exp_sp(P->bits, F->mctx->minfo), w;
        ulong fpw = FLINT_BITS / P->bits;
        const ulong * m = P->exps;   /* the leading term */

        for (w = N - 1; w >= 0; w--)
        {
            if (m[w] != 0)
            {
                slong b = FLINT_BIT_COUNT(m[w]) - 1;
                slong idx = w * fpw + b / P->bits;
                return nvars - 1 - idx;
            }
        }
        return nvars;
    }

    e = (nvars <= 16) ? e_small : flint_malloc(sizeof(ulong) * nvars);
    fmpz_mpoly_get_term_exp_ui(e, P, 0, F->mctx);
    for (v = 0; v < nvars; v++)
    {
        if (e[v] != 0)
        {
            res = v;
            break;
        }
    }
    if (e != e_small)
        flint_free(e);
    return res;
}

slong
gr_tower_flat_def_order(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    /* (the generator of definition order d is the variable cap - 1 - d:
       the highest definition order occurring is that of the most
       significant variable, which occurs in a leading term) */
    slong v = FLINT_MIN(_lead_var(fmpz_mpoly_q_numref(x), F), _lead_var(fmpz_mpoly_q_denref(x), F));

    if (v >= F->mctx->minfo->nvars)
        return -1;

    return F->cap - 1 - v;
}

int
gr_tower_flat_compose(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, fmpz_mpoly_q_struct ** imgs, gr_tower_flat_t F)
{
    slong nvars = x_mctx->minfo->nvars, i;
    int * used;
    int all_poly = 1;
    int status = GR_SUCCESS;

    gr_tower_flat_ensure(F);

    used = flint_malloc(sizeof(int) * nvars);
    fmpz_mpoly_q_used_vars(used, x, x_mctx);

    for (i = 0; i < nvars; i++)
    {
        if (!used[i])
            continue;
        if (imgs[i] == NULL)
        {
            flint_free(used);
            return GR_DOMAIN;
        }
        if (!fmpz_mpoly_is_one(fmpz_mpoly_q_denref(imgs[i]), F->mctx))
            all_poly = 0;
    }

    /* renaming of variables: all images are generators */
    {
        int all_gens = 1;
        slong * c = flint_malloc(sizeof(slong) * nvars);
        ulong * exp = flint_malloc(sizeof(ulong) * F->cap);

        for (i = 0; i < nvars && all_gens; i++)
        {
            c[i] = 0;
            if (!used[i])
                continue;
            if (!fmpz_mpoly_is_gen(fmpz_mpoly_q_numref(imgs[i]), -1, F->mctx) || !fmpz_mpoly_is_one(fmpz_mpoly_q_denref(imgs[i]), F->mctx))
            {
                all_gens = 0;
            }
            else
            {
                slong j;
                fmpz_mpoly_get_term_exp_ui(exp, fmpz_mpoly_q_numref(imgs[i]), 0, F->mctx);
                for (j = 0; j < F->cap; j++)
                    if (exp[j] == 1)
                        c[i] = j;
            }
        }

        if (all_gens)
        {
            fmpz_mpoly_t t;
            fmpz_mpoly_init(t, F->mctx);
            fmpz_mpoly_compose_fmpz_mpoly_gen(t, fmpz_mpoly_q_numref(x), c, x_mctx, F->mctx);
            fmpz_mpoly_swap(fmpz_mpoly_q_numref(res), t, F->mctx);
            fmpz_mpoly_compose_fmpz_mpoly_gen(t, fmpz_mpoly_q_denref(x), c, x_mctx, F->mctx);
            fmpz_mpoly_swap(fmpz_mpoly_q_denref(res), t, F->mctx);
            fmpz_mpoly_clear(t, F->mctx);
            flint_free(c);
            flint_free(exp);
            flint_free(used);
            return gr_tower_flat_reduce(res, F);
        }

        flint_free(c);
        flint_free(exp);
    }

    if (all_poly)
    {
        fmpz_mpoly_struct ** C = flint_malloc(sizeof(fmpz_mpoly_struct *) * nvars);
        fmpz_mpoly_t zero;
        fmpz_mpoly_q_t t;

        fmpz_mpoly_init(zero, F->mctx);
        for (i = 0; i < nvars; i++)
            C[i] = (used[i] ? fmpz_mpoly_q_numref(imgs[i]) : zero);

        fmpz_mpoly_q_init(t, F->mctx);
        if (!fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_q_numref(t), fmpz_mpoly_q_numref(x), C, x_mctx, F->mctx) ||
            !fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_q_denref(t), fmpz_mpoly_q_denref(x), C, x_mctx, F->mctx))
        {
            status = GR_UNABLE;
        }
        else
        {
            fmpz_mpoly_q_canonicalise(t, F->mctx);
            fmpz_mpoly_q_swap(res, t, F->mctx);
            status = gr_tower_flat_reduce(res, F);
        }
        fmpz_mpoly_q_clear(t, F->mctx);
        fmpz_mpoly_clear(zero, F->mctx);
        flint_free(C);
    }
    else
    {
        /* term by term (rare path: images with denominators) */
        const fmpz_mpoly_struct * parts[2];
        fmpz_mpoly_q_t outs[2], acc, term;
        ulong * exp = flint_malloc(sizeof(ulong) * nvars);
        fmpz_t c;
        int p;

        parts[0] = fmpz_mpoly_q_numref(x);
        parts[1] = fmpz_mpoly_q_denref(x);
        fmpz_mpoly_q_init(acc, F->mctx);
        fmpz_mpoly_q_init(term, F->mctx);
        fmpz_mpoly_q_init(outs[0], F->mctx);
        fmpz_mpoly_q_init(outs[1], F->mctx);
        fmpz_init(c);

        for (p = 0; p < 2; p++)
        {
            slong j, len = fmpz_mpoly_length(parts[p], x_mctx);

            fmpz_mpoly_q_zero(acc, F->mctx);
            for (j = 0; j < len; j++)
            {
                fmpz_mpoly_get_term_exp_ui(exp, parts[p], j, x_mctx);
                fmpz_mpoly_get_term_coeff_fmpz(c, parts[p], j, x_mctx);
                fmpz_mpoly_q_set_fmpz(term, c, F->mctx);
                for (i = 0; i < nvars; i++)
                {
                    ulong e = exp[i];
                    while (e > 0)
                    {
                        fmpz_mpoly_q_mul(term, term, imgs[i], F->mctx);
                        e--;
                    }
                }
                fmpz_mpoly_q_add(acc, acc, term, F->mctx);
            }
            fmpz_mpoly_q_set(outs[p], acc, F->mctx);
        }

        fmpz_mpoly_q_div(res, outs[0], outs[1], F->mctx);
        status = gr_tower_flat_reduce(res, F);

        fmpz_clear(c);
        fmpz_mpoly_q_clear(acc, F->mctx);
        fmpz_mpoly_q_clear(term, F->mctx);
        fmpz_mpoly_q_clear(outs[0], F->mctx);
        fmpz_mpoly_q_clear(outs[1], F->mctx);
        flint_free(exp);
    }

    flint_free(used);
    return status;
}

/* String with the tower's generator names, for an element in the given context. */
char *
_gr_tower_flat_get_str(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t mctx, gr_tower_t T)
{
    slong cap, v;
    slong * var_gid;
    char ** vars;
    char * s, * s2, * res;

    if (!_gr_tower_flat_find_layout(&cap, &var_gid, mctx, &T->flat))
    {
        cap = mctx->minfo->nvars;
        var_gid = NULL;
    }

    vars = flint_malloc(sizeof(char *) * cap);
    for (v = 0; v < cap; v++)
    {
        slong d = (var_gid == NULL || var_gid[v] < 0) ? -1 : gr_tower_gid_order(T, var_gid[v]);
        vars[v] = (d >= 0) ? T->gens[d].name : (char *) "?";
    }

    s = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(x), (const char **) vars, mctx);
    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(x), mctx))
    {
        res = s;
    }
    else
    {
        /* parentheses only where needed: around a numerator of several
           terms, and around a denominator which is not a single factor
           (an integer, or a generator or a power of one) */
        const fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(x);
        const fmpz_mpoly_struct * den = fmpz_mpoly_q_denref(x);
        int pn = (num->length > 1);
        int pd = 1;

        if (fmpz_mpoly_is_fmpz(den, mctx))
            pd = 0;
        else if (den->length == 1 && fmpz_is_one(den->coeffs))
        {
            int * used = flint_calloc(mctx->minfo->nvars, sizeof(int));
            slong k, cnt = 0;
            fmpz_mpoly_used_vars(used, den, mctx);
            for (k = 0; k < mctx->minfo->nvars; k++)
                cnt += used[k];
            flint_free(used);
            pd = (cnt != 1);
        }

        s2 = fmpz_mpoly_get_str_pretty(den, (const char **) vars, mctx);
        res = flint_malloc(strlen(s) + strlen(s2) + 8);
        strcpy(res, pn ? "(" : "");
        strcat(res, s);
        strcat(res, pn ? ")/" : "/");
        strcat(res, pd ? "(" : "");
        strcat(res, s2);
        strcat(res, pd ? ")" : "");
        flint_free(s);
        flint_free(s2);
    }

    flint_free(vars);
    return res;
}

/*
    Whether the (reduced, nonzero) numerator of x is nonzero by
    Lindemann's theorem: it involves only the generators up to some
    transcendental t = exp(u), log(u) or pi, those before t being
    algebraic with proven moduli (a number field, in which the reduced
    representation is canonical); exp(u) and log(u) with u algebraic
    (u != 0, resp. u != 0, 1, as adjoined) are transcendental
    (Hermite-Lindemann), so a polynomial in t over the number field with
    a nonzero coefficient is nonzero. Decides exp(10^-10000) - 1 without
    a 33000-bit evaluation.
*/
static int
_flat_num_lindemann_nonzero(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong d, top = -1;
    int * used, ok = 1;

    used = flint_malloc(sizeof(int) * F->mctx->minfo->nvars);
    fmpz_mpoly_used_vars(used, fmpz_mpoly_q_numref(x), F->mctx);
    for (d = 0; d < T->num_gens; d++)
        if (used[GR_TOWER_FLAT_VAR_D(F, d)])
            top = d;
    flint_free(used);

    if (top < 0)
        return 0;
    {
        const gr_tower_gen_struct * g = T->gens + top;
        if (g->kind != GR_TOWER_EXP && g->kind != GR_TOWER_LOG && g->kind != GR_TOWER_PI)
            return 0;
    }
    for (d = 0; d < top && ok; d++)
        if (T->gens[d].kind != GR_TOWER_ALGEBRAIC || T->gens[d].status != GR_TOWER_STATUS_PROVEN)
            ok = 0;
    return ok;
}

/*
    The exact test of the numerator of x (reduced) via the nested
    representation, at the lowest possible level: a code of
    _gr_tower_field_is_zero_at (the tower may be refined). The monomial
    content of x is removed first when it is nonzero: x = sqrt(u) (exp(a)
    - exp(b)) is tested as exp(a) - exp(b), below sqrt(u), whose modulus
    can be large.
*/
int
_gr_tower_flat_nested_zero_code(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_t M;
    fmpz_mpoly_q_t y;
    gr_ctx_struct * Fk;
    gr_ptr t;
    slong k;
    int split = 0, code;

    fmpz_mpoly_init(M, F->mctx);
    fmpz_mpoly_q_init(y, F->mctx);
    fmpz_mpoly_term_content(M, fmpz_mpoly_q_numref(x), F->mctx);
    if (!fmpz_mpoly_is_fmpz(M, F->mctx))
    {
        acb_t w;
        acb_init(w);
        fmpz_mpoly_one(fmpz_mpoly_q_denref(y), F->mctx);
        fmpz_mpoly_set(fmpz_mpoly_q_numref(y), M, F->mctx);
        if (gr_tower_flat_get_acb(w, y, GR_TOWER_DEFAULT_PREC, F) == GR_SUCCESS && !acb_contains_zero(w) &&
            fmpz_mpoly_divides(fmpz_mpoly_q_numref(y), fmpz_mpoly_q_numref(x), M, F->mctx))
            split = 1;
        acb_clear(w);
    }
    if (!split)
        fmpz_mpoly_q_set(y, x, F->mctx);

    k = gr_tower_flat_alg_level(y, F);
    Fk = gr_tower_field_at(T, k);
    GR_TMP_INIT(t, Fk);
    if (gr_tower_flat_poly_get_nested_at(t, fmpz_mpoly_q_numref(y), k, F) == GR_SUCCESS)
    {
        fmpz_mpoly_q_clear(y, F->mctx);
        fmpz_mpoly_clear(M, F->mctx);
        code = _gr_tower_field_is_zero_at(t, k, T);
    }
    else
    {
        fmpz_mpoly_q_clear(y, F->mctx);
        fmpz_mpoly_clear(M, F->mctx);
        code = GR_TOWER_UNKNOWN;
    }
    GR_TMP_CLEAR(t, Fk);
    return code;
}

/*
    Complete zero test of the numerator of x (in place: x is reduced,
    and the tower may be refined).
*/
static truth_t
_flat_num_is_zero(fmpz_mpoly_q_t x, gr_tower_flat_t F, int fixed)
{
    gr_tower_struct * T = F->T;
    truth_t res;
    slong prec;
    acb_t z;

    if (gr_tower_flat_reduce(x, F) != GR_SUCCESS)
        return T_UNKNOWN;

    if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), F->mctx))
        return T_TRUE;

    /* with all moduli proved irreducible and all transcendental
       generators proved transcendental, the reduced representation is
       canonical */
    if (_gr_tower_all_proven(T))
        return T_FALSE;

    acb_init(z);
    for (prec = GR_TOWER_DEFAULT_PREC; prec <= GR_TOWER_OPTION(T, GR_TOWER_OPT_PREC_LIMIT); prec *= 2)
    {
        if (gr_tower_flat_get_acb(z, x, prec, F) == GR_SUCCESS && !acb_contains_zero(z))
        {
            acb_clear(z);
            return T_FALSE;
        }
    }
    acb_clear(z);

    if (_flat_num_lindemann_nonzero(x, F))
        return T_FALSE;

    /* square roots of integers which are Gauss sums in the roots of
       unity before them (sqrt(3) next to zeta_3 and i) become steps of
       degree one, rather than being discovered by dynamic evaluation
       (the flat element survives the change of the tower) */
    if (!fixed)
    {
        slong d;
        int changed = 0;
        for (d = 0; d < T->num_gens; d++)
            changed |= _gr_tower_gauss_sqrt_gen(T, d);
        if (changed)
        {
            if (gr_tower_flat_reduce(x, F) != GR_SUCCESS)
                return T_UNKNOWN;
            if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), F->mctx))
                return T_TRUE;
        }
    }

    /* exact test via the nested representation */
    {
        int code = _gr_tower_flat_nested_zero_code(x, F);

        if (code == GR_TOWER_ZERO)
            res = T_TRUE;
        else if (code == GR_TOWER_NONZERO)
            res = T_FALSE;
        else if (code == GR_TOWER_FIELD_NONZERO)
            res = _gr_tower_flat_involves_conjectural(x, F) ? T_UNKNOWN : T_FALSE;
        else
            res = T_UNKNOWN;
    }

    /* the tower may have been refined: x is reduced (a zero is set to
       zero directly, rather than reduced modulo moduli which may have
       become large: a generator found of degree one over the steps
       below it, sqrt(c) = a large expression) */
    if (res == T_TRUE)
        fmpz_mpoly_q_zero(x, F->mctx);
    else if (gr_tower_flat_reduce(x, F) != GR_SUCCESS)
        res = T_UNKNOWN;

    if (res == T_UNKNOWN && _gr_tower_has_conjectural(T))
    {
        if (fixed)
        {
            /* on a copy of the tower, so that neither the tower nor the
               polynomial context changes */
            gr_tower_t U;
            fmpz_mpoly_q_t xu;
            gr_tower_init(U, T->consts);
            gr_tower_set(U, T);
            U->keep_retired = 0;   /* (internal: no nested elements across rebuilds) */
            gr_tower_flat_ensure(&U->flat);
            fmpz_mpoly_q_init(xu, U->flat.mctx);
            _gr_tower_flat_transport(xu, x, F, &U->flat);
            res = _gr_tower_decide_zero_flat(xu, &U->flat, U->num_gens, 0);
            fmpz_mpoly_q_clear(xu, U->flat.mctx);
            gr_tower_clear(U);
        }
        else
        {
            /* search for relations between the transcendental generators,
               eliminating them in place (flat elements stay valid, but the
               polynomial context may change: x is not touched afterwards) */
            res = _gr_tower_decide_zero_flat(x, F, T->num_gens, 0);
        }
    }

    return res;
}

truth_t
gr_tower_flat_num_is_zero(fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    return _flat_num_is_zero(x, F, 0);
}

/* (for contexts whose elements cannot follow a change of the layout) */
truth_t
_gr_tower_flat_num_is_zero_fixed(fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    return _flat_num_is_zero(x, F, 1);
}
