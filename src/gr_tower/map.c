/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "acb.h"
#include "fmpq.h"
#include "fmpz_factor.h"
#include "fmpz_poly.h"
#include "arb_fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_mpoly_q.h"
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

int
gr_tower_promote(gr_ptr res, gr_srcptr x, slong j, slong k, gr_tower_t T)
{
    int status = GR_SUCCESS;

    if (j == k)
        return gr_set(res, x, gr_tower_field_at(T, k));

    if (j > k)
        return GR_DOMAIN;

    if (j + 1 == k)
        return gr_set_other(res, x, gr_tower_field_at(T, j), gr_tower_field_at(T, k));

    {
        gr_ctx_struct * mid = gr_tower_field_at(T, j + 1);
        gr_ptr t;
        GR_TMP_INIT(t, mid);
        status |= gr_set_other(t, x, gr_tower_field_at(T, j), mid);
        status |= gr_tower_promote(res, t, j + 1, k, T);
        GR_TMP_CLEAR(t, mid);
    }

    return status;
}

/* -------------------------------------------------------------------- */
/* maps                                                                  */
/* -------------------------------------------------------------------- */

void
gr_tower_map_init(gr_tower_map_t map, gr_tower_t source, gr_tower_t target)
{
    slong i;

    gr_tower_flat_ensure(&target->flat);

    map->source = source;
    map->target = target;
    map->length = 0;
    map->alloc = FLINT_MAX(source->num_gens, 1);
    map->mctx = target->flat.mctx;
    map->images = flint_malloc(map->alloc * sizeof(fmpz_mpoly_q_struct));
    for (i = 0; i < map->alloc; i++)
        fmpz_mpoly_q_init(map->images + i, map->mctx);
}

void
gr_tower_map_clear(gr_tower_map_t map)
{
    slong i;
    for (i = 0; i < map->alloc; i++)
        fmpz_mpoly_q_clear(map->images + i, map->mctx);
    flint_free(map->images);
}

void
gr_tower_map_sync(gr_tower_map_t map)
{
    gr_tower_flat_struct * F = &map->target->flat;
    slong i;

    gr_tower_flat_ensure(F);

    if (map->mctx == F->mctx)
        return;

    for (i = 0; i < map->alloc; i++)
    {
        fmpz_mpoly_q_t t;
        fmpz_mpoly_q_init(t, F->mctx);
        gr_tower_flat_convert(t, map->images + i, map->mctx, F);
        fmpz_mpoly_q_clear(map->images + i, map->mctx);
        fmpz_mpoly_q_init(map->images + i, F->mctx);
        fmpz_mpoly_q_swap(map->images + i, t, F->mctx);
        fmpz_mpoly_q_clear(t, F->mctx);
    }

    map->mctx = F->mctx;
}

/* Deep copy; res must be initialized (with any source and target). */
void
gr_tower_map_set(gr_tower_map_t res, const gr_tower_map_t map)
{
    slong i;

    if (res == map)
        return;

    gr_tower_map_clear(res);
    res->source = map->source;
    res->target = map->target;
    res->length = map->length;
    res->alloc = map->alloc;
    res->mctx = map->mctx;
    res->images = flint_malloc(map->alloc * sizeof(fmpz_mpoly_q_struct));
    for (i = 0; i < map->alloc; i++)
    {
        fmpz_mpoly_q_init(res->images + i, res->mctx);
        fmpz_mpoly_q_set(res->images + i, map->images + i, res->mctx);
    }
}

int
gr_tower_map_fit_length(gr_tower_map_t map, slong len)
{
    slong i;

    if (len <= map->alloc)
        return GR_SUCCESS;

    len = FLINT_MAX(len, 2 * map->alloc);
    map->images = flint_realloc(map->images, len * sizeof(fmpz_mpoly_q_struct));
    for (i = map->alloc; i < len; i++)
        fmpz_mpoly_q_init(map->images + i, map->mctx);
    map->alloc = len;
    return GR_SUCCESS;
}

int
gr_tower_map_set_image_flat(gr_tower_map_t map, slong d, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t mctx)
{
    if (d < 0)
        return GR_DOMAIN;

    gr_tower_map_sync(map);
    gr_tower_map_fit_length(map, d + 1);
    gr_tower_flat_convert(map->images + d, x, mctx, &map->target->flat);
    map->length = FLINT_MAX(map->length, d + 1);
    return GR_SUCCESS;
}

int
gr_tower_map_set_image(gr_tower_map_t map, slong d, gr_srcptr x, slong k)
{
    int status;

    if (d < 0)
        return GR_DOMAIN;

    gr_tower_map_sync(map);
    gr_tower_map_fit_length(map, d + 1);
    status = gr_tower_flat_set_nested_at(map->images + d, x, k, &map->target->flat);
    if (status == GR_SUCCESS)
        map->length = FLINT_MAX(map->length, d + 1);
    return status;
}

int
gr_tower_map_set_inclusion(gr_tower_map_t map)
{
    gr_tower_struct * S = map->source;
    gr_tower_struct * T = map->target;
    gr_tower_flat_struct * F = &T->flat;
    slong d;

    if (S->num_gens > T->num_gens)
        return GR_DOMAIN;

    gr_tower_map_sync(map);
    gr_tower_map_fit_length(map, S->num_gens);

    for (d = 0; d < S->num_gens; d++)
    {
        if ((GR_TOWER_GEN(S, d)->kind == GR_TOWER_ALGEBRAIC) != (GR_TOWER_GEN(T, d)->kind == GR_TOWER_ALGEBRAIC))
            return GR_DOMAIN;

        fmpz_mpoly_q_gen(map->images + d, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
    }

    map->length = S->num_gens;
    return GR_SUCCESS;
}

int
gr_tower_map_apply_flat(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, gr_tower_map_t map)
{
    gr_tower_struct * S = map->source;
    gr_tower_flat_struct * FS = &S->flat;
    fmpz_mpoly_q_struct ** imgs;
    fmpz_mpoly_q_t xs;
    slong d, nvars;
    int status;

    gr_tower_map_sync(map);
    gr_tower_flat_ensure(FS);

    /* bring x to the current context of the source */
    fmpz_mpoly_q_init(xs, FS->mctx);
    gr_tower_flat_convert(xs, x, x_mctx, FS);

    nvars = FS->cap;
    imgs = flint_calloc(nvars, sizeof(fmpz_mpoly_q_struct *));
    for (d = 0; d < map->length; d++)
        imgs[GR_TOWER_FLAT_VAR_D(FS, d)] = map->images + d;

    status = gr_tower_flat_compose(res, xs, FS->mctx, imgs, &map->target->flat);

    flint_free(imgs);
    fmpz_mpoly_q_clear(xs, FS->mctx);
    return status;
}

int
gr_tower_map_apply_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_map_t map)
{
    gr_tower_struct * S = map->source;
    gr_tower_struct * T = map->target;
    fmpz_mpoly_q_t xs, xt;
    fmpz_mpoly_ctx_struct * sctx, * tctx;
    int status;

    gr_tower_map_sync(map);
    gr_tower_flat_ensure(&S->flat);
    sctx = S->flat.mctx;
    tctx = map->mctx;

    fmpz_mpoly_q_init(xs, sctx);
    fmpz_mpoly_q_init(xt, tctx);
    status = gr_tower_flat_set_nested_at(xs, x, k, &S->flat);
    if (status == GR_SUCCESS)
        status = gr_tower_map_apply_flat(xt, xs, sctx, map);
    if (status == GR_SUCCESS)
        status = gr_tower_flat_get_nested_at(res, xt, T->length, &T->flat);
    fmpz_mpoly_q_clear(xs, sctx);
    fmpz_mpoly_q_clear(xt, tctx);

    return status;
}

int
gr_tower_map_apply(gr_ptr res, gr_srcptr x, gr_tower_map_t map)
{
    return gr_tower_map_apply_at(res, x, map->source->length, map);
}

int
gr_tower_map_apply_poly(gr_poly_t res, const gr_poly_t f, gr_tower_map_t map)
{
    gr_ctx_struct * top = gr_tower_field(map->target);
    gr_ctx_struct * stop = gr_tower_field(map->source);
    slong i, len = f->length;
    int status = GR_SUCCESS;

    gr_poly_fit_length(res, len, top);
    for (i = 0; i < len && status == GR_SUCCESS; i++)
        status |= gr_tower_map_apply(gr_poly_coeff_ptr(res, i, top), gr_poly_coeff_srcptr(f, i, stop), map);
    _gr_poly_set_length(res, len, top);
    _gr_poly_normalise(res, top);

    return status;
}

/* -------------------------------------------------------------------- */
/* absorption                                                            */
/* -------------------------------------------------------------------- */

/* Sets image d of the map to the generator with definition order du of
   the target. */
static void
_map_set_image_gen(gr_tower_map_t map, slong d, slong du)
{
    gr_tower_flat_struct * F = &map->target->flat;

    gr_tower_map_sync(map);
    gr_tower_map_fit_length(map, d + 1);
    fmpz_mpoly_q_gen(map->images + d, GR_TOWER_FLAT_VAR_D(F, du), F->mctx);
}

/* Absorbs the transcendental generator t_j of B (with definition order d). */
/*
    A generator of U (of any current kind) defined as the function kind
    with parameter param of the argument arg (an element of U in the
    context *actx), identified by an exact comparison of the arguments;
    returns its definition order, or -1. The zero tests may restructure
    U, in which case arg is converted to the new context of the map and
    *actx is updated.
*/
static slong
_absorb_find_definition(gr_tower_t U, gr_tower_map_t map, int kind, slong param,
    fmpz_mpoly_q_t arg, fmpz_mpoly_ctx_struct ** actx)
{
    gr_tower_flat_struct * F = &U->flat;
    slong du;

    for (du = 0; du < U->num_gens; du++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(U, du);
        slong gid = g->gid;
        fmpz_mpoly_q_t diff, a;
        truth_t eq;

        if (g->arg.mctx == NULL)
            continue;
        if (!((g->kind == kind) || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == kind)))
            continue;
        if (GR_TOWER_KIND_IS_SPECIAL(kind) && g->def_param != param)
            continue;

        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(diff, F->mctx);
        fmpz_mpoly_q_init(a, F->mctx);
        gr_tower_flat_convert(diff, &g->arg.data, g->arg.mctx, F);
        gr_tower_flat_convert(a, arg, *actx, F);
        fmpz_mpoly_q_sub(diff, diff, a, F->mctx);
        eq = gr_tower_flat_num_is_zero(diff, F);
        fmpz_mpoly_q_clear(diff, F->mctx);
        fmpz_mpoly_q_clear(a, F->mctx);

        if (eq == T_TRUE)
            return gr_tower_gid_order(U, gid);

        /* the zero test may have restructured U */
        gr_tower_map_sync(map);
        if (*actx != map->mctx)
        {
            fmpz_mpoly_q_t tmp;
            fmpz_mpoly_q_init(tmp, map->mctx);
            gr_tower_flat_convert(tmp, arg, *actx, &U->flat);
            fmpz_mpoly_q_clear(arg, *actx);
            *actx = map->mctx;
            fmpz_mpoly_q_init(arg, *actx);
            fmpz_mpoly_q_swap(arg, tmp, *actx);
            fmpz_mpoly_q_clear(tmp, *actx);
        }
    }

    return -1;
}

static int
_absorb_trans(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, slong j, slong d)
{
    gr_tower_gen_struct * t = GR_TOWER_TRANS(B, j - 1);
    slong ju;
    int status = GR_SUCCESS;

    /* already present in U (as any kind of generator)? */
    if (t->def_id != 0)
    {
        slong du = gr_tower_find_def_order(U, t->def_id);
        if (du >= 0)
        {
            _map_set_image_gen(map, d, du);
            return GR_SUCCESS;
        }
    }

    /* pi is pi; named constants likewise */
    if (t->kind == GR_TOWER_PI || t->kind == GR_TOWER_CONSTANT)
    {
        slong du;
        for (du = 0; du < U->num_gens; du++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_GEN(U, du);
            if ((g->kind == t->kind || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == t->kind)) &&
                (t->kind == GR_TOWER_PI || g->def_param == t->def_param))
            {
                _map_set_image_gen(map, d, du);
                return GR_SUCCESS;
            }
        }
    }

    {
        /* adjoin, with the image of the argument */
        fmpz_mpoly_q_t arg;
        fmpz_mpoly_ctx_struct * actx;
        gr_tower_gen_struct * nt;

        gr_tower_map_sync(map);
        actx = map->mctx;
        fmpz_mpoly_q_init(arg, actx);

        if (GR_TOWER_KIND_HAS_ARG(t->kind))
            status = gr_tower_map_apply_flat(arg, &t->arg.data, t->arg.mctx, map);

        /* the same function of the same argument already in U (possibly
           as a generator which became algebraic) */
        if (status == GR_SUCCESS && GR_TOWER_KIND_HAS_ARG(t->kind))
        {
            slong du = _absorb_find_definition(U, map, t->kind, t->def_param, arg, &actx);
            if (du >= 0)
            {
                fmpz_mpoly_q_clear(arg, actx);
                _map_set_image_gen(map, d, du);
                return GR_SUCCESS;
            }
        }

        if (status == GR_SUCCESS)
        {
            if (t->kind == GR_TOWER_PI)
                status = gr_tower_adjoin_pi(U, t->name);
            else if (t->kind == GR_TOWER_EXP)
                status = gr_tower_adjoin_exp_flat(U, arg, actx, t->name);
            else if (t->kind == GR_TOWER_LOG)
                status = gr_tower_adjoin_log_flat(U, arg, actx, t->name);
            else if (t->kind == GR_TOWER_TAN)
                status = gr_tower_adjoin_tan_flat(U, arg, actx, t->name);
            else if (t->kind == GR_TOWER_ATAN)
                status = gr_tower_adjoin_atan_flat(U, arg, actx, t->name);
            else if (GR_TOWER_KIND_IS_SPECIAL(t->kind))
                status = _gr_tower_adjoin_special_flat_nocheck(U, t->kind, t->def_param, arg, actx, t->name);
            else
                status = GR_UNABLE;
        }

        if (status == GR_DOMAIN && GR_TOWER_KIND_IS_SPECIAL(t->kind) &&
            !(t->kind == GR_TOWER_ERF || (t->kind == GR_TOWER_LAMBERTW && t->def_param == 0) || t->kind == GR_TOWER_POLYLOG))
        {
            /* (a pole: cannot occur for a valid generator) */
            fmpz_mpoly_q_clear(arg, actx);
            return GR_UNABLE;
        }

        if (status == GR_DOMAIN && t->kind != GR_TOWER_PI)
        {
            /* the image of the argument is trivial in U: exp(0) = 1 or
               log(1) = 0, tan(0) = atan(0) = 0, erf(0) = W_0(0) =
               Li_s(0) = 0 (log(0) cannot occur for a valid generator) */
            gr_tower_map_sync(map);
            gr_tower_map_fit_length(map, d + 1);
            if (t->kind == GR_TOWER_EXP)
                fmpz_mpoly_q_one(map->images + d, map->mctx);
            else
                fmpz_mpoly_q_zero(map->images + d, map->mctx);
            fmpz_mpoly_q_clear(arg, actx);
            return GR_SUCCESS;
        }

        fmpz_mpoly_q_clear(arg, actx);

        if (status != GR_SUCCESS)
            return status;

        ju = U->num_trans;
        nt = GR_TOWER_TRANS(U, ju - 1);
        nt->def_id = t->def_id;
        nt->def_param = t->def_param;
        nt->status = t->status;
        if (t->enclosure_prec > nt->enclosure_prec)
        {
            acb_set(&nt->enclosure, &t->enclosure);
            nt->enclosure_prec = t->enclosure_prec;
        }
    }

    _map_set_image_gen(map, d, GR_TOWER_TRANS(U, ju - 1)->def_order);
    return GR_SUCCESS;
}

/*
    The root of the monic quadratic q = X^2 + q1 X + q0 over the top field
    of U selected by the enclosure z, if it lies in that field:
    -q1/2 +- sqrt(q1^2/4 - q0).
*/
/* (the target z and the two roots, for the selection) */
typedef struct { const acb_struct * z; gr_srcptr r1; gr_srcptr r2; gr_tower_struct * U; } _quadratic_root_arg;

static int
_quadratic_root_get(acb_t z, int which, slong prec, void * arg)
{
    _quadratic_root_arg * a = (_quadratic_root_arg *) arg;
    if (which < 0)
    {
        acb_set(z, a->z);
        return GR_SUCCESS;
    }
    return gr_tower_get_acb(z, (which == 0) ? a->r1 : a->r2, prec, a->U);
}

static int
_quadratic_root(gr_ptr res, const gr_poly_t q, const acb_t z, gr_tower_t U)
{
    gr_ctx_struct * top = gr_tower_field(U);
    gr_ptr disc, s, c, r1, r2;
    int status;

    GR_TMP_INIT5(disc, s, c, r1, r2, top);

    status = gr_div_ui(c, gr_poly_coeff_srcptr(q, 1, top), 2, top);
    status |= gr_neg(c, c, top);                    /* c = -q1/2 */
    status |= gr_sqr(disc, c, top);
    status |= gr_sub(disc, disc, gr_poly_coeff_srcptr(q, 0, top), top);
    if (status == GR_SUCCESS)
        status = gr_tower_sqrt(s, disc, U);
    if (status == GR_SUCCESS)
    {
        status = gr_add(r1, c, s, top);
        status |= gr_sub(r2, c, s, top);
    }

    if (status == GR_SUCCESS)
    {
        _quadratic_root_arg a = { z, r1, r2, U };
        int which = _gr_tower_select_candidate(_quadratic_root_get, &a, U);
        if (which >= 0)
            status = gr_set(res, (which == 0) ? r1 : r2, top);
        else
            status = GR_UNABLE;
    }

    GR_TMP_CLEAR5(disc, s, c, r1, r2, top);
    return status;
}

/*
    Whether the generator is a root of unity exp(2 pi i / n) (returns 1)
    or the principal n-th root of a positive integer p (returns 2, and
    sets p), as recorded by its definition.
*/
int
_gr_tower_gen_const_root(const gr_tower_gen_struct * g, fmpz_t p, ulong * n)
{
    if (g->kind != GR_TOWER_ALGEBRAIC)
        return 0;

    if (g->def_kind == GR_TOWER_ROOT_OF_UNITY)
    {
        *n = g->def_param;
        return 1;
    }

    if (g->def_kind == GR_TOWER_ROOT && g->arg.mctx != NULL &&
        fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(&g->arg.data), g->arg.mctx) &&
        fmpz_mpoly_is_one(fmpz_mpoly_q_denref(&g->arg.data), g->arg.mctx))
    {
        fmpz_mpoly_get_fmpz(p, fmpz_mpoly_q_numref(&g->arg.data), g->arg.mctx);
        if (fmpz_sgn(p) <= 0)
            return 0;
        *n = g->def_param;
        return 2;
    }

    return 0;
}

/*
    The generator of U of the given kind (root of unity, or root of the
    integer p) whose order is the largest power of the prime l; returns
    its definition order or -1, and sets *f to the exponent of l.
*/
static slong
_find_prime_power_root(ulong * f, gr_tower_t U, int kind, const fmpz_t p, ulong l)
{
    fmpz_t pu;
    ulong m, best_f = 0;
    slong du, best = -1;

    fmpz_init(pu);
    for (du = 0; du < U->num_gens; du++)
    {
        if (_gr_tower_gen_const_root(GR_TOWER_GEN(U, du), pu, &m) == kind && (kind == 1 || fmpz_equal(p, pu)))
        {
            ulong e = 0;
            while (m % l == 0)
            {
                m /= l;
                e++;
            }
            if (m == 1 && e > best_f)
            {
                best_f = e;
                best = du;
            }
        }
    }
    fmpz_clear(pu);

    *f = best_f;
    return best;
}

/* res = res * (generator du of F)^e, reduced */
static void
_mul_gen_pow(fmpz_mpoly_q_t res, gr_tower_flat_t F, slong du, ulong e)
{
    fmpz_mpoly_t t;
    fmpz_mpoly_init(t, F->mctx);
    fmpz_mpoly_gen(t, GR_TOWER_FLAT_VAR_D(F, du), F->mctx);
    if (e != 1)
        fmpz_mpoly_pow_ui(t, t, e, F->mctx);
    fmpz_mpoly_mul(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), t, F->mctx);
    GR_MUST_SUCCEED(gr_tower_flat_reduce(res, F));
    fmpz_mpoly_clear(t, F->mctx);
}

/*
    The principal square root of the nonzero integer c through the roots
    of unity among the generators of F with definition order < limit
    (quadratic Gauss sums: sqrt(p*) = sum_a (a/p) zeta_p^a for odd primes
    p, p* = (-1)^((p-1)/2) p; sqrt(2) = zeta_8 + zeta_8^-1; sqrt(-1) =
    zeta_4). Writes it to res (in F's current context) and returns 1, or
    returns 0 if the roots of unity needed are missing. Without this,
    sqrt(3) next to zeta_3 and i, or sqrt(2) next to zeta_8, is a
    generator of degree 2 over a field containing it. (Primes up to
    GR_TOWER_OPT_GAUSS_SUM_LIMIT.)
*/
static slong
_find_root_of_unity_multiple(ulong * n, gr_tower_t T, ulong m, slong limit)
{
    fmpz_t pu;
    slong d, best = -1;
    ulong k, best_n = 0;

    fmpz_init(pu);
    for (d = 0; d < limit && d < T->num_gens; d++)
    {
        if (_gr_tower_gen_const_root(GR_TOWER_GEN(T, d), pu, &k) == 1 && k % m == 0 &&
            (best < 0 || k < best_n))
        {
            best = d;
            best_n = k;
        }
    }
    fmpz_clear(pu);
    *n = best_n;
    return best;
}


/* (the principal square root of c, and +/- the candidate res, for the
   sign selection) */
typedef struct { const fmpz_mpoly_q_struct * res; const fmpz * c; gr_tower_flat_struct * F; } _cyclo_sqrt_sign_arg;

static int
_cyclo_sqrt_sign_get(acb_t z, int which, slong prec, void * arg)
{
    _cyclo_sqrt_sign_arg * a = (_cyclo_sqrt_sign_arg *) arg;
    int st = GR_SUCCESS;
    if (which < 0)
    {
        acb_set_fmpz(z, a->c);
        acb_sqrt(z, z, prec);
    }
    else
    {
        st = gr_tower_flat_get_acb(z, a->res, prec, a->F);
        if (st == GR_SUCCESS && which == 1)
            acb_neg(z, z);
    }
    return st;
}

int
_gr_tower_cyclotomic_sqrt(fmpz_mpoly_q_t res, const fmpz_t c, gr_tower_flat_t F, slong limit)
{
    gr_tower_struct * T = F->T;
    fmpz_factor_t fac;
    fmpz_mpoly_t g, t;
    slong i, d;
    ulong n;
    int neg, ok = 1;

    if (fmpz_is_zero(c))
        return 0;

    gr_tower_flat_ensure(F);

    /* (only primes up to the Gauss sum limit can occur to odd powers:
       trial division by them, the cofactor having to be a square, rather
       than a complete factorization of a possibly huge radicand) */
    fmpz_factor_init(fac);
    {
        fmpz_t r, s;
        ulong lim = FLINT_MIN(GR_TOWER_OPTION(F->T, GR_TOWER_OPT_GAUSS_SUM_LIMIT), 1000000);
        fmpz_init(r);
        fmpz_init(s);
        _gr_tower_fmpz_factor_trial(fac, r, c, n_prime_pi(lim));
        fac->sign = fmpz_sgn(c);
        if (!fmpz_is_one(r))
        {
            if (!fmpz_is_square(r))
            {
                fmpz_clear(r);
                fmpz_clear(s);
                fmpz_factor_clear(fac);
                return 0;
            }
            fmpz_sqrt(s, r);
            _fmpz_factor_append(fac, s, 2);
        }
        fmpz_clear(r);
        fmpz_clear(s);
    }
    neg = (fac->sign < 0);

    fmpz_mpoly_init(g, F->mctx);
    fmpz_mpoly_init(t, F->mctx);
    fmpz_mpoly_q_one(res, F->mctx);

    for (i = 0; i < fac->num && ok; i++)
    {
        if (fac->exp[i] >= 2)
        {
            fmpz_t pp;
            fmpz_init(pp);
            fmpz_pow_ui(pp, fac->p + i, fac->exp[i] / 2);
            fmpz_mpoly_q_mul_fmpz(res, res, pp, F->mctx);
            fmpz_clear(pp);
        }
        if (fac->exp[i] % 2 == 0)
            continue;

        if (fmpz_cmp_ui(fac->p + i, GR_TOWER_OPTION(F->T, GR_TOWER_OPT_GAUSS_SUM_LIMIT)) > 0)
        {
            ok = 0;
            break;
        }

        if (fmpz_equal_ui(fac->p + i, 2))
        {
            /* zeta_8 + zeta_8^7 */
            d = _find_root_of_unity_multiple(&n, T, 8, limit);
            if (d < 0)
            {
                ok = 0;
                break;
            }
            fmpz_mpoly_gen(g, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            fmpz_mpoly_pow_ui(g, g, n / 8, F->mctx);
            fmpz_mpoly_pow_ui(t, g, 7, F->mctx);
            fmpz_mpoly_add(g, g, t, F->mctx);
        }
        else
        {
            ulong p = fmpz_get_ui(fac->p + i), a;
            fmpz_mpoly_t w;

            d = _find_root_of_unity_multiple(&n, T, p, limit);
            if (d < 0)
            {
                ok = 0;
                break;
            }
            fmpz_mpoly_init(w, F->mctx);
            fmpz_mpoly_gen(w, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            fmpz_mpoly_pow_ui(w, w, n / p, F->mctx);
            fmpz_mpoly_zero(g, F->mctx);
            fmpz_mpoly_one(t, F->mctx);
            for (a = 1; a < p; a++)
            {
                fmpz_mpoly_mul(t, t, w, F->mctx);
                if (n_jacobi(a, p) > 0)
                    fmpz_mpoly_add(g, g, t, F->mctx);
                else
                    fmpz_mpoly_sub(g, g, t, F->mctx);
            }
            fmpz_mpoly_clear(w, F->mctx);
            if (p % 4 == 3)
                neg = !neg;
        }

        fmpz_mpoly_mul(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), g, F->mctx);
        if (gr_tower_flat_reduce(res, F) != GR_SUCCESS)
            ok = 0;
    }

    if (ok && neg)
    {
        d = _find_root_of_unity_multiple(&n, T, 4, limit);
        if (d < 0)
            ok = 0;
        else
        {
            fmpz_mpoly_gen(g, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            fmpz_mpoly_pow_ui(g, g, n / 4, F->mctx);
            fmpz_mpoly_mul(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), g, F->mctx);
            if (gr_tower_flat_reduce(res, F) != GR_SUCCESS)
                ok = 0;
        }
    }

    /* the principal root: res^2 = c exactly, the sign is chosen
       numerically (the two candidates are well separated) */
    if (ok)
    {
        _cyclo_sqrt_sign_arg a = { res, c, F };
        int which = _gr_tower_select_candidate(_cyclo_sqrt_sign_get, &a, F->T);
        if (which == 1)
            fmpz_mpoly_q_neg(res, res, F->mctx);
        ok = (which >= 0);
    }

    fmpz_mpoly_clear(g, F->mctx);
    fmpz_mpoly_clear(t, F->mctx);
    fmpz_factor_clear(fac);
    return ok;
}

/*
    If the generator with definition order i of T is a square root of an
    integer with an unsplit dynamic modulus x^2 - p and the roots of
    unity before it express sqrt(p) (a Gauss sum), makes it algebraic
    of degree one (sqrt(3) next to zeta_3 and i). Returns 1 if so.
*/
int
_gr_tower_gauss_sqrt_gen(gr_tower_t T, slong i)
{
    gr_tower_gen_struct * h = T->gens + i;
    fmpz_t p;
    ulong n;
    int found = 0;

    if (h->kind != GR_TOWER_ALGEBRAIC || h->status == GR_TOWER_STATUS_PROVEN)
        return 0;

    fmpz_init(p);
    if (_gr_tower_gen_const_root(h, p, &n) == 2 && n == 2 &&
        gr_tower_step_minpoly(T, h->index)->length == 3)
    {
        fmpz_mpoly_q_struct mm[2];
        gr_tower_flat_ensure(&T->flat);
        fmpz_mpoly_q_init(mm + 0, T->flat.mctx);
        fmpz_mpoly_q_init(mm + 1, T->flat.mctx);
        found = _gr_tower_cyclotomic_sqrt(mm + 0, p, &T->flat, i);
        if (found)
        {
            fmpz_mpoly_q_neg(mm + 0, mm + 0, T->flat.mctx);
            fmpz_mpoly_q_one(mm + 1, T->flat.mctx);
            found = _gr_tower_make_algebraic(T, i, mm, 2, T->flat.mctx, GR_TOWER_STATUS_PROVEN);
        }
        fmpz_mpoly_q_clear(mm + 0, T->flat.mctx);
        fmpz_mpoly_q_clear(mm + 1, T->flat.mctx);
    }
    fmpz_clear(p);
    return found;
}

/*
    Moves the generators tan(pi/M) (def_kind GR_TOWER_TAN_PI) whose
    moduli have base field coefficients to the front of U, then proves
    (or refines, by Trager's method) the steps after them; returns 1 if
    any was moved.
*/
static int
_tan_pi_to_front(gr_tower_t U)
{
    slong d, k;
    int moved = 0;

    for (d = 1; d < U->num_gens; d++)
    {
        gr_tower_gen_struct * g = GR_TOWER_GEN(U, d);
        const gr_poly_struct * m;
        gr_ctx_struct * below;
        slong i;
        int base = 1;

        if (g->kind != GR_TOWER_ALGEBRAIC || g->def_kind != GR_TOWER_TAN_PI)
            continue;

        m = gr_tower_step_minpoly(U, g->index);
        below = gr_tower_field_at(U, g->index - 1);
        for (i = 0; i < m->length && base; i++)
        {
            fmpz_mpoly_q_t f;
            gr_tower_flat_ensure(&U->flat);
            fmpz_mpoly_q_init(f, U->flat.mctx);
            base = (gr_tower_flat_set_nested_at(f, gr_poly_coeff_srcptr(m, i, below), g->index - 1, &U->flat) == GR_SUCCESS) &&
                fmpz_mpoly_q_is_fmpq(f, U->flat.mctx);
            fmpz_mpoly_q_clear(f, U->flat.mctx);
        }
        if (!base)
            continue;

        _gr_tower_move_to_front(U, d);
        if (GR_TOWER_HAS_CAP(U, GR_TOWER_CAP_RATIONAL))
            _gr_tower_gen_set_status(GR_TOWER_GEN(U, 0), GR_TOWER_STATUS_PROVEN);
        moved = 1;
    }

    if (moved)
    {
        for (k = 1; k <= U->length; k++)
        {
            gr_tower_gen_struct * g = GR_TOWER_STEP(U, k - 1);
            if (g->status != GR_TOWER_STATUS_PROVEN)
            {
                if (GR_TOWER_HAS_CAP(U, GR_TOWER_CAP_RATIONAL))
                    gr_tower_prove_step_modular(U, k, GR_TOWER_OPTION(U, GR_TOWER_OPT_MODULAR_TRIES));
                if (g->status != GR_TOWER_STATUS_PROVEN)
                    gr_tower_prove_step_trager(U, k, GR_TOWER_OPTION(U, GR_TOWER_OPT_TRAGER_DEGREE_LIMIT));
            }
        }
    }

    return moved;
}

/*
    zeta_n through a root of unity generator of U of an order N divisible
    by n (zeta_n = zeta_N^(N/n)): returns 1 with the image set. Otherwise,
    if adjoin is set (composite roots of unity), adjoins zeta_L for L the
    least common multiple of n and the orders of the roots of unity of U,
    in front, makes those roots its powers (degree one), and returns 1
    with the image set; returns 0 if neither applies (-1 on failure).
*/
static int
_composite_root_image(fmpz_mpoly_q_t img, fmpz_mpoly_ctx_struct ** img_ctx, gr_tower_t U, ulong n, int adjoin, ulong new_def_id, const char * name, int gauss)
{
    gr_tower_flat_struct * F = &U->flat;
    ulong N, L;
    slong du, i, nold = 0, * old_gid, gid_new;

    du = _find_root_of_unity_multiple(&N, U, n, U->num_gens);

    if (du < 0)
    {
        gr_tower_gen_struct * ng;

        if (!adjoin)
            return 0;

        old_gid = flint_malloc(sizeof(slong) * FLINT_MAX(U->num_gens, 1));
        L = n;
        for (i = 0; i < U->num_gens; i++)
        {
            const gr_tower_gen_struct * h = GR_TOWER_GEN(U, i);
            if (h->kind == GR_TOWER_ALGEBRAIC && h->def_kind == GR_TOWER_ROOT_OF_UNITY && h->def_param > 0)
            {
                if (L != 0)
                {
                    /* (capped: 0 when beyond the order limit) */
                    ulong a = L / n_gcd(L, h->def_param);
                    L = (a > (ulong) GR_TOWER_OPTION(U, GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT) / (ulong) h->def_param) ? 0 : a * h->def_param;
                }
                old_gid[nold++] = h->gid;
            }
        }

        /* (beyond the degree limit, the prime power decomposition) */
        if (L == 0 || L > (ulong) GR_TOWER_OPTION(U, GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT) || n_euler_phi(L) > (ulong) GR_TOWER_OPTION(U, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT))
        {
            flint_free(old_gid);
            return 0;
        }

        if (gr_tower_adjoin_root_of_unity(U, L, (L == n) ? name : NULL) != GR_SUCCESS)
        {
            flint_free(old_gid);
            return -1;
        }
        ng = GR_TOWER_GEN(U, U->num_gens - 1);
        gid_new = ng->gid;
        if (L == n)
            ng->def_id = new_def_id;
        _gr_tower_move_to_front(U, U->num_gens - 1);
        if (GR_TOWER_HAS_CAP(U, GR_TOWER_CAP_RATIONAL))
            _gr_tower_gen_set_status(GR_TOWER_GEN(U, 0), GR_TOWER_STATUS_PROVEN);
        if (!gauss)
            _tan_pi_to_front(U);

        /* the old roots of unity are powers of the new one */
        for (i = 0; i < nold; i++)
        {
            fmpz_mpoly_q_struct mm[2];
            slong dold = gr_tower_gid_order(U, old_gid[i]);
            ulong No = GR_TOWER_GEN(U, dold)->def_param;

            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_init(mm + 0, F->mctx);
            fmpz_mpoly_q_init(mm + 1, F->mctx);
            fmpz_mpoly_q_one(mm + 0, F->mctx);
            _mul_gen_pow(mm + 0, F, gr_tower_gid_order(U, gid_new), L / No);
            fmpz_mpoly_q_neg(mm + 0, mm + 0, F->mctx);
            fmpz_mpoly_q_one(mm + 1, F->mctx);
            _gr_tower_make_algebraic(U, dold, mm, 2, F->mctx, GR_TOWER_STATUS_PROVEN);
            fmpz_mpoly_q_clear(mm + 0, F->mctx);
            fmpz_mpoly_q_clear(mm + 1, F->mctx);
        }
        flint_free(old_gid);

        /* generators whose proofs were dropped by the move may lie in the
           new field */
        _gr_tower_express_downgraded(U, 1, U->num_gens - 1);

        du = gr_tower_gid_order(U, gid_new);
        N = L;
    }

    gr_tower_flat_ensure(F);
    if (*img_ctx != F->mctx)
    {
        fmpz_mpoly_q_clear(img, *img_ctx);
        *img_ctx = F->mctx;
        fmpz_mpoly_q_init(img, *img_ctx);
    }
    fmpz_mpoly_q_one(img, F->mctx);
    _mul_gen_pow(img, F, du, N / n);
    return 1;
}

/*
    A root of unity of order n, or the principal n-th root of a positive
    integer p, in terms of roots of the same kind of prime power orders
    q_i = l_i^{e_i} (n = prod q_i): with c_i = (n/q_i)^{-1} mod q_i,
    sum c_i / q_i = 1/n + k for an integer k >= 0, so that
    root_n = p^{-k} prod root_{q_i}^{c_i} (p = 1 for roots of unity).

    Writes the image of the generator g in F to img (in F's current
    context), where the prime power roots are found among the
    generators of F (any generator of the same kind whose order is a
    power of l_i divisible by q_i). Returns 0 if some prime power root is
    missing (and, if adjoin is set, adjoins the missing ones -- or
    replaces a generator of too small an order l^f by one of order q_i,
    re-expressed as its power -- instead, so that the image can be
    written; the new generators are placed in front, where their minimal
    polynomials over the base field are known to be irreducible).
*/
static int _composite_root_image(fmpz_mpoly_q_t img, fmpz_mpoly_ctx_struct ** img_ctx, gr_tower_t U, ulong n, int adjoin, ulong new_def_id, const char * name, int gauss);

static int
_structured_image(fmpz_mpoly_q_t img, fmpz_mpoly_ctx_struct ** img_ctx, gr_tower_t U, const gr_tower_gen_struct * g, int adjoin, ulong new_def_id, int gauss, int composite)
{
    gr_tower_flat_struct * F = &U->flat;
    fmpz_t p;
    ulong n;
    n_factor_t fac;
    slong i;
    int kind, ok = 1;
    ulong ksum = 0;

    fmpz_init(p);
    kind = _gr_tower_gen_const_root(g, p, &n);
    if (kind == 0)
    {
        fmpz_clear(p);
        return 0;
    }

    /* a root of unity of an order dividing that of a generator of U is
       a power of it (also when the generator has a composite order);
       with composite roots of unity, a missing one is adjoined as one
       generator */
    if (kind == 1)
    {
        int r = _composite_root_image(img, img_ctx, U, n, adjoin && composite, new_def_id, g->name, gauss);
        if (r != 0)
        {
            fmpz_clear(p);
            return (r == 1);
        }
    }

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    gr_tower_flat_ensure(F);
    if (*img_ctx != F->mctx)
    {
        fmpz_mpoly_q_clear(img, *img_ctx);
        *img_ctx = F->mctx;
        fmpz_mpoly_q_init(img, *img_ctx);
    }
    fmpz_mpoly_q_one(img, F->mctx);

    for (i = 0; i < fac.num && ok; i++)
    {
        ulong l = fac.p[i], e = fac.exp[i], q = n_pow(l, e), f, c;
        slong du;

        du = _find_prime_power_root(&f, U, kind, p, l);

        /* sqrt(p) through roots of unity of U (a Gauss sum) */
        if (gauss && du < 0 && kind == 2 && q == 2)
        {
            fmpz_mpoly_q_t s;
            int found;
            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_init(s, F->mctx);
            found = _gr_tower_cyclotomic_sqrt(s, p, F, U->num_gens);
            if (found)
            {
                /* (c = (n/2)^{-1} mod 2 = 1) */
                if (*img_ctx != F->mctx)
                {
                    fmpz_mpoly_q_t t;
                    fmpz_mpoly_q_init(t, F->mctx);
                    gr_tower_flat_convert(t, img, *img_ctx, F);
                    fmpz_mpoly_q_clear(img, *img_ctx);
                    *img_ctx = F->mctx;
                    fmpz_mpoly_q_init(img, *img_ctx);
                    fmpz_mpoly_q_swap(img, t, *img_ctx);
                    fmpz_mpoly_q_clear(t, *img_ctx);
                }
                fmpz_mpoly_q_mul(img, img, s, F->mctx);
                GR_MUST_SUCCEED(gr_tower_flat_reduce(img, F));
                ksum += n / q;
            }
            fmpz_mpoly_q_clear(s, F->mctx);
            if (found)
                continue;
        }

        if (du < 0 || f < e)
        {
            slong gid_old = (du >= 0) ? GR_TOWER_GEN(U, du)->gid : -1, gid_new;
            gr_tower_gen_struct * ng;
            int status;

            /* a root of an integer with no root of that integer present
               is left to the general absorption (it may be expressible
               through roots of the factors of the integer, e.g.
               sqrt(6) = sqrt(2) sqrt(3)); roots of unity are always
               adjoined here (the representation stays canonical) */
            if (!adjoin || (kind == 2 && du < 0))
            {
                ok = 0;
                break;
            }

            /* (a generator of prime power order is g itself: same name) */
            if (kind == 1)
                status = gr_tower_adjoin_root_of_unity(U, q, (fac.num == 1) ? g->name : NULL);
            else
                status = gr_tower_adjoin_root_fmpz(U, p, q, (fac.num == 1) ? g->name : NULL);

            if (status != GR_SUCCESS)
            {
                ok = 0;
                break;
            }

            ng = GR_TOWER_GEN(U, U->num_gens - 1);
            gid_new = ng->gid;
            if (fac.num == 1)
                ng->def_id = new_def_id;   /* the new generator is g itself */
            _gr_tower_move_to_front(U, U->num_gens - 1);
            ng = GR_TOWER_GEN(U, 0);
            if (GR_TOWER_HAS_CAP(U, GR_TOWER_CAP_RATIONAL))
                _gr_tower_gen_set_status(ng, GR_TOWER_STATUS_PROVEN);
            /* (real views) the tangent generators tan(pi/M), defined over
               the base field, are kept in front of the roots of unity
               rather than re-expressed through them below: a real
               subfield stays represented by its real generator */
            if (!gauss)
                _tan_pi_to_front(U);
            /* generators whose proofs were dropped by the move may lie in
               the new field (sqrt(2) in Q(zeta_8)) */
            _gr_tower_express_downgraded(U, 1, U->num_gens - 1);

            if (gid_old >= 0)
            {
                /* the old generator is the (q / l^f)-th power of the new one */
                fmpz_mpoly_q_struct mm[2];

                gr_tower_flat_ensure(F);
                fmpz_mpoly_q_init(mm + 0, F->mctx);
                fmpz_mpoly_q_init(mm + 1, F->mctx);
                fmpz_mpoly_q_one(mm + 0, F->mctx);
                _mul_gen_pow(mm + 0, F, 0, q / n_pow(l, f));
                fmpz_mpoly_q_neg(mm + 0, mm + 0, F->mctx);
                fmpz_mpoly_q_one(mm + 1, F->mctx);

                _gr_tower_make_algebraic(U, gr_tower_gid_order(U, gid_old), mm, 2, F->mctx, GR_TOWER_STATUS_PROVEN);

                fmpz_mpoly_q_clear(mm + 0, F->mctx);
                fmpz_mpoly_q_clear(mm + 1, F->mctx);
            }

            /* the layout of U changed: img must follow */
            gr_tower_flat_ensure(F);
            if (*img_ctx != F->mctx)
            {
                fmpz_mpoly_q_t t;
                fmpz_mpoly_q_init(t, F->mctx);
                gr_tower_flat_convert(t, img, *img_ctx, F);
                fmpz_mpoly_q_clear(img, *img_ctx);
                *img_ctx = F->mctx;
                fmpz_mpoly_q_init(img, *img_ctx);
                fmpz_mpoly_q_swap(img, t, *img_ctx);
                fmpz_mpoly_q_clear(t, *img_ctx);
            }

            du = gr_tower_gid_order(U, gid_new);
            f = e;
        }

        /* c = (n/q)^{-1} mod q; the generator has order l^f >= q */
        c = (fac.num == 1) ? 1 : n_invmod((n / q) % q, q);
        ksum += c * (n / q);
        _mul_gen_pow(img, F, du, c * n_pow(l, f - e));
    }

    if (ok && kind == 2 && fac.num > 1)
    {
        /* p^{-k} with k = (sum c_i n/q_i - 1)/n */
        ulong k = (ksum - 1) / n;
        if (k > 0)
        {
            fmpz_t pk;
            fmpz_init(pk);
            fmpz_pow_ui(pk, p, k);
            fmpz_mpoly_q_div_fmpz(img, img, pk, F->mctx);
            fmpz_clear(pk);
        }
    }

    fmpz_clear(p);
    return ok;
}

/* Whether the roots of prime power orders needed to express g are
   present in U (without building the image). */
int
_gr_tower_structured_present(gr_tower_t U, const gr_tower_gen_struct * g, int gauss)
{
    fmpz_t p;
    ulong n;
    n_factor_t fac;
    slong i;
    int kind, ok = 1;

    fmpz_init(p);
    kind = _gr_tower_gen_const_root(g, p, &n);
    if (kind == 0)
    {
        fmpz_clear(p);
        return 0;
    }

    if (kind == 1)
    {
        ulong N;
        if (_find_root_of_unity_multiple(&N, U, n, U->num_gens) >= 0)
        {
            fmpz_clear(p);
            return 1;
        }
    }

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    for (i = 0; i < fac.num && ok; i++)
    {
        ulong f;
        if (_find_prime_power_root(&f, U, kind, p, fac.p[i]) < 0 || f < fac.exp[i])
        {
            ok = 0;
            if (gauss && kind == 2 && fac.p[i] == 2 && fac.exp[i] == 1)
            {
                /* a Gauss sum (see _structured_image) */
                fmpz_mpoly_q_t s;
                gr_tower_flat_ensure(&U->flat);
                fmpz_mpoly_q_init(s, U->flat.mctx);
                ok = _gr_tower_cyclotomic_sqrt(s, p, &U->flat, U->num_gens);
                fmpz_mpoly_q_clear(s, U->flat.mctx);
            }
        }
    }

    fmpz_clear(p);
    return ok;
}

int
_gr_tower_structured_image(fmpz_mpoly_q_t img, gr_tower_t U, const gr_tower_gen_struct * g, int gauss)
{
    fmpz_mpoly_ctx_struct * ctx = U->flat.mctx;
    return _structured_image(img, &ctx, U, g, 0, 0, gauss, 0);
}

/*
    Absorption by definition (see _structured_image). Returns 1 if the
    image was set.
*/
static int
_absorb_structured(gr_tower_t U, gr_tower_map_t map, const gr_tower_gen_struct * step, slong d, int gauss, int composite)
{
    fmpz_mpoly_q_t img;
    int ok;

    if (step->def_kind != GR_TOWER_ROOT_OF_UNITY && step->def_kind != GR_TOWER_ROOT)
        return 0;

    fmpz_mpoly_ctx_struct * img_ctx;

    gr_tower_flat_ensure(&U->flat);
    img_ctx = U->flat.mctx;
    fmpz_mpoly_q_init(img, img_ctx);
    ok = _structured_image(img, &img_ctx, U, step, 1, step->def_id, gauss, composite);
    if (ok)
    {
        gr_tower_map_sync(map);
        gr_tower_map_fit_length(map, d + 1);
        gr_tower_flat_convert(map->images + d, img, img_ctx, &U->flat);
    }
    fmpz_mpoly_q_clear(img, img_ctx);
    return ok;
}

/*
    A root of the same order of the same argument (the image of the
    argument in U being equal to the argument of a root generator of U)
    is that generator. Returns 1 if the image was set.
*/
static int
_absorb_same_root(gr_tower_t U, gr_tower_map_t map, const gr_tower_gen_struct * step, slong d)
{
    gr_tower_flat_struct * F = &U->flat;
    fmpz_mpoly_q_t arg;
    slong du;
    int found = 0;

    if (step->def_kind != GR_TOWER_ROOT || step->arg.mctx == NULL)
        return 0;

    gr_tower_map_sync(map);
    fmpz_mpoly_q_init(arg, map->mctx);
    if (gr_tower_map_apply_flat(arg, &step->arg.data, step->arg.mctx, map) != GR_SUCCESS)
    {
        fmpz_mpoly_q_clear(arg, map->mctx);
        return 0;
    }

    for (du = 0; du < U->num_gens && !found; du++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(U, du);
        slong gid = g->gid;
        fmpz_mpoly_q_t diff;

        if (g->kind != GR_TOWER_ALGEBRAIC || g->def_kind != GR_TOWER_ROOT || g->def_param != step->def_param || g->arg.mctx == NULL)
            continue;

        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(diff, F->mctx);
        gr_tower_flat_convert(diff, &g->arg.data, g->arg.mctx, F);
        {
            fmpz_mpoly_q_t a;
            fmpz_mpoly_q_init(a, F->mctx);
            gr_tower_flat_convert(a, arg, map->mctx, F);
            fmpz_mpoly_q_sub(diff, diff, a, F->mctx);
            fmpz_mpoly_q_clear(a, F->mctx);
        }
        if (gr_tower_flat_num_is_zero(diff, F) == T_TRUE)
        {
            _map_set_image_gen(map, d, gr_tower_gid_order(U, gid));
            found = 1;
        }
        fmpz_mpoly_q_clear(diff, F->mctx);
    }

    gr_tower_map_sync(map);
    fmpz_mpoly_q_clear(arg, map->mctx);
    return found;
}

/*
    Divides the monic polynomial q over the top field of U by x - beta for
    every algebraic generator beta of U which is a root of q (filtered
    numerically, then decided exactly).
*/
static int
_divide_out_conjugates(gr_poly_t q, const acb_t z, slong * match, const fmpz_poly_t origin, gr_tower_t U)
{
    gr_ctx_struct * top = gr_tower_field(U);
    slong k;
    int status = GR_SUCCESS;

    *match = -1;

    /* the conjugates already present usually form a chain, each one's
       modulus being the origin polynomial divided by x - (the previous
       ones): then q is the modulus of the last conjugate divided by
       x - (that conjugate), a single synthetic division */
    if (origin != NULL)
    {
        slong count = 0, last = -1;

        for (k = 1; k <= U->length; k++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_STEP(U, k - 1);
            if (g->origin != NULL && fmpz_poly_equal(origin, g->origin))
            {
                count++;
                last = k;
            }
        }

        if (last > 0 && q->length == origin->length &&
            gr_tower_step_minpoly(U, last)->length == q->length - (count - 1))
        {
            const gr_poly_struct * m = gr_tower_step_minpoly(U, last);
            gr_ctx_struct * below = gr_tower_field_at(U, last - 1);
            gr_ctx_struct * Fk = gr_tower_field_at(U, last);
            gr_ptr a, beta, u;
            slong i, n = m->length - 1;
            acb_t w;
            int which = 0;

            /* is the root being absorbed the last conjugate itself?
               (decided from the cached enclosures, which isolate roots
               of the origin polynomial; refining them can be expensive,
               so an ambiguous case goes to the general path below) */
            acb_init(w);
            acb_set(w, &GR_TOWER_STEP(U, last - 1)->enclosure);
            if (acb_contains(z, w))
                which = 1;
            else if (!acb_overlaps(z, w))
                which = -1;
            acb_clear(w);

            if (which == 1)
            {
                *match = last;
                return GR_SUCCESS;
            }

          if (which == -1)
          {
            GR_TMP_INIT(a, Fk);
            GR_TMP_INIT2(beta, u, top);
            status = gr_gen(a, Fk);
            status |= gr_tower_promote(beta, a, last, U->length, U);

            /* q = m / (x - beta) over the top field */
            gr_poly_fit_length(q, n + 1, top);
            for (i = 0; i <= n && status == GR_SUCCESS; i++)
                status |= gr_set_other(gr_poly_coeff_ptr(q, i, top), gr_poly_coeff_srcptr(m, i, below), below, top);
            _gr_poly_set_length(q, n + 1, top);
            if (status == GR_SUCCESS)
                status = gr_poly_div_root(q, u, q, beta, top);

            GR_TMP_CLEAR(a, Fk);
            GR_TMP_CLEAR2(beta, u, top);
            return status;
          }
        }
    }

    for (k = 1; k <= U->length && q->length >= 3 && status == GR_SUCCESS; k++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_STEP(U, k - 1);
        int candidate = 0, known = 0;

        /* a generator with the same origin polynomial is a root of q,
           as long as q is the origin polynomial itself (afterwards, a
           root already divided out is not a root of q any more) */
        if (origin != NULL && g->origin != NULL && fmpz_poly_equal(origin, g->origin))
        {
            candidate = 1;
            known = (q->length == origin->length);
        }

        /* otherwise a numerical filter */
        if (!candidate && origin == NULL)
        {
            acb_poly_t f;
            acb_t zb, v;
            slong i, prec = GR_TOWER_DEFAULT_PREC;

            acb_poly_init(f);
            acb_init(zb);
            acb_init(v);
            for (i = 0; i < q->length && status == GR_SUCCESS; i++)
            {
                status |= gr_tower_get_acb(v, gr_poly_coeff_srcptr(q, i, top), prec, U);
                acb_poly_set_coeff_acb(f, i, v);
            }
            if (status == GR_SUCCESS && gr_tower_step_get_acb(zb, U, k, prec) == GR_SUCCESS)
            {
                acb_poly_evaluate(v, f, zb, prec);
                candidate = acb_contains_zero(v);
            }
            acb_poly_clear(f);
            acb_clear(zb);
            acb_clear(v);
        }

        if (status != GR_SUCCESS || !candidate)
            continue;

        /* exact test (unless known) and synthetic division */
        {
            gr_ctx_struct * Fk = gr_tower_field_at(U, k);
            gr_ptr a, beta, val;
            int is_root = known;
            GR_TMP_INIT(a, Fk);
            GR_TMP_INIT2(beta, val, top);
            status = gr_gen(a, Fk);
            status |= gr_tower_promote(beta, a, k, U->length, U);
            if (status == GR_SUCCESS && !known)
            {
                status |= gr_poly_evaluate(val, q, beta, top);
                is_root = (status == GR_SUCCESS && gr_tower_is_zero(val, U) == T_TRUE);
            }
            if (status == GR_SUCCESS && is_root)
            {
                /* beta is a root of q: the root being absorbed (its
                   isolating enclosure z contains beta), or a conjugate
                   to divide out */
                gr_ptr u;
                acb_t w;
                slong prec;
                int which = 0;   /* 1: beta is the root, -1: conjugate */

                acb_init(w);
                for (prec = GR_TOWER_DEFAULT_PREC; prec <= 4096 && which == 0; prec *= 2)
                {
                    if (gr_tower_step_get_acb(w, U, k, prec) != GR_SUCCESS)
                        break;
                    if (acb_contains(z, w))
                        which = 1;
                    else if (!acb_overlaps(z, w))
                        which = -1;
                }
                acb_clear(w);

                if (which == 1)
                {
                    *match = k;
                    GR_TMP_CLEAR(a, Fk);
                    GR_TMP_CLEAR2(beta, val, top);
                    break;
                }
                if (which == 0)
                {
                    status = GR_UNABLE;
                    GR_TMP_CLEAR(a, Fk);
                    GR_TMP_CLEAR2(beta, val, top);
                    break;
                }

                GR_TMP_INIT(u, top);
                status |= gr_poly_div_root(q, u, q, beta, top);
                GR_TMP_CLEAR(u, top);
            }
            GR_TMP_CLEAR(a, Fk);
            GR_TMP_CLEAR2(beta, val, top);
        }
    }

    return status;
}

/* Absorbs the algebraic generator a_k of B (with definition order d). */
/*
    For the algebraic generators with definition orders in [from, to]
    whose moduli are not known irreducible (typically after a generator
    was moved in front of them), tries to express them in the field
    generated by the preceding generators (lattice reduction, for
    fields of moderate degree); a generator found there gets a linear
    modulus. Elements in the nested representation become invalid.
*/
void
_gr_tower_express_downgraded(gr_tower_t T, slong from, slong to)
{
    slong i;

    for (i = from; i <= to && i < T->num_gens; i++)
    {
        gr_tower_gen_struct * h = T->gens + i;
        gr_tower_t P;
        gr_ctx_struct * ptop;
        gr_poly_t q;
        gr_ptr img;
        slong k, kp, j, deg, D;
        int found = 0;

        if (h->kind != GR_TOWER_ALGEBRAIC || h->status != GR_TOWER_STATUS_DYNAMIC)
            continue;

        if (gr_tower_prefix_length(T, i) != gr_tower_prefix_length(T, i + 1) - 1)
            continue;

        k = gr_tower_prefix_length(T, i);   /* algebraic level of the prefix */
        deg = gr_poly_quotient_ctx_degree(GR_TOWER_STEP_CTX(T, k));
        if (deg < 2)
            continue;

        /* the prefix as a tower of its own */
        gr_tower_init(P, T->consts);
        gr_tower_set_prefix(P, T, i);
        P->keep_retired = 0;   /* (internal: no nested elements across rebuilds) */
        D = gr_tower_degree(P);
        if (D > GR_TOWER_OPTION(T, GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT) || P->num_trans != T->num_trans)
        {
            gr_tower_clear(P);
            continue;
        }
        kp = P->length;
        ptop = gr_tower_field(P);
        gr_poly_init(q, ptop);
        GR_TMP_INIT(img, ptop);

        /* the modulus of h, over the prefix field */
        {
            const gr_poly_struct * m = gr_poly_quotient_ctx_modulus(GR_TOWER_STEP_CTX(T, k));
            gr_ctx_struct * below = gr_tower_field_at(T, k);
            fmpz_mpoly_q_t cf, cp;
            int ok = 1;

            gr_tower_flat_ensure(&T->flat);
            gr_tower_flat_ensure(&P->flat);
            fmpz_mpoly_q_init(cf, T->flat.mctx);
            fmpz_mpoly_q_init(cp, P->flat.mctx);
            for (j = 0; j < m->length && ok; j++)
            {
                gr_ptr c;
                if (gr_tower_flat_set_nested_at(cf, gr_poly_coeff_srcptr(m, j, below), k, &T->flat) != GR_SUCCESS)
                {
                    ok = 0;
                    break;
                }
                _gr_tower_flat_transport(cp, cf, &T->flat, &P->flat);
                GR_TMP_INIT(c, ptop);
                if (gr_tower_flat_get_nested_at(c, cp, kp, &P->flat) != GR_SUCCESS ||
                    gr_poly_set_coeff_scalar(q, j, c, ptop) != GR_SUCCESS)
                    ok = 0;
                GR_TMP_CLEAR(c, ptop);
            }
            fmpz_mpoly_q_clear(cf, T->flat.mctx);
            fmpz_mpoly_q_clear(cp, P->flat.mctx);

            if (ok && q->length == m->length && !gr_tower_poly_no_roots_modular(q, P, GR_TOWER_OPTION(T, GR_TOWER_OPT_NO_ROOTS_TRIES)))
                found = (gr_tower_express_limit(img, q, &h->enclosure,
                            GR_TOWER_MERGE_EXPRESS_PREC(T, D), P) == GR_SUCCESS);
        }

        if (found)
        {
            fmpz_mpoly_q_struct mm[2];
            fmpz_mpoly_q_t cp;

            gr_tower_flat_ensure(&P->flat);
            fmpz_mpoly_q_init(cp, P->flat.mctx);
            if (gr_tower_flat_set_nested_at(cp, img, kp, &P->flat) == GR_SUCCESS)
            {
                gr_tower_flat_ensure(&T->flat);
                fmpz_mpoly_q_init(mm + 0, T->flat.mctx);
                fmpz_mpoly_q_init(mm + 1, T->flat.mctx);
                _gr_tower_flat_transport(mm + 0, cp, &P->flat, &T->flat);
                fmpz_mpoly_q_neg(mm + 0, mm + 0, T->flat.mctx);
                fmpz_mpoly_q_one(mm + 1, T->flat.mctx);
                _gr_tower_make_algebraic(T, i, mm, 2, T->flat.mctx, GR_TOWER_STATUS_PROVEN);
                fmpz_mpoly_q_clear(mm + 0, T->flat.mctx);
                fmpz_mpoly_q_clear(mm + 1, T->flat.mctx);
            }
            fmpz_mpoly_q_clear(cp, P->flat.mctx);
        }

        GR_TMP_CLEAR(img, ptop);
        gr_poly_clear(q, ptop);
        gr_tower_clear(P);
    }
}

/*
    Whether g is an n-th root of a prime p (any n >= 2) and the algebraic
    generators of U with definition order < limit are roots of unity of
    orders not divisible by p and roots of integers a of orders m with
    p not dividing a m (with proven moduli): the field K they generate
    (with any transcendental generators, over which algebraic elements
    stay in K) is then unramified at p (the discriminants of the moduli
    divide powers of m and a m), so that X^n - p is Eisenstein at every
    prime of K above p, hence irreducible over K. In particular the root
    does not lie in K (dft example: sqrt(11) over Q(zeta_16, sqrt 2, ...,
    sqrt 7), where the lattice search would otherwise be tried). With
    require_proven = 0 the statuses of the steps of U are not checked:
    the root is then still not in the field generated by the numbers
    (only the irreducibility over the quotient ring needs them).
*/
int
_gr_tower_unramified_radical(const gr_tower_gen_struct * g, const gr_tower_t U, slong limit, int require_proven)
{
    fmpz_t p, pu;
    ulong n, m;
    slong du;
    int ok = 0;

    if (!GR_TOWER_HAS_CAP(U, GR_TOWER_CAP_RATIONAL))
        return 0;

    fmpz_init(p);
    fmpz_init(pu);

    if (_gr_tower_gen_const_root(g, p, &n) == 2 && n >= 2 && fmpz_sgn(p) > 0 &&
        fmpz_is_prime(p) == 1)
    {
        ok = 1;
        for (du = 0; du < limit && ok; du++)
        {
            const gr_tower_gen_struct * h = GR_TOWER_GEN(U, du);
            int kind;

            if (h->kind != GR_TOWER_ALGEBRAIC)
                continue;

            if (require_proven && h->status != GR_TOWER_STATUS_PROVEN)
            {
                ok = 0;
                break;
            }

            /* (p divides the order m iff p <= m and m = 0 mod p) */
            kind = _gr_tower_gen_const_root(h, pu, &m);
            if (kind == 1)
                ok = (fmpz_cmp_ui(p, m) > 0 || m % fmpz_get_ui(p) != 0);
            else if (kind == 2 && fmpz_sgn(pu) > 0)
                ok = !fmpz_divisible(pu, p) && (fmpz_cmp_ui(p, m) > 0 || m % fmpz_get_ui(p) != 0);
            else
                ok = 0;
        }
    }

    fmpz_clear(p);
    fmpz_clear(pu);
    return ok;
}

/*
    Structural irreducibility of the minimal polynomial of an algebraic
    generator g over the field generated by the generators of U with
    definition order < limit:

    - g and all algebraic generators of the prefix (i excepted) are
      principal roots of positive integers, with moduli x^n - a
      irreducible over Q (proven status), the radicands pairwise coprime
      (or equal, with coprime orders): the real radicals are then
      multiplicatively independent
      modulo Q^*, so that the degree of the field they generate is the
      product of the orders (Mordell's theorem on real radicals, which
      generalizes Besicovitch's theorem on roots of distinct primes);

    - a root of unity of order N over a field generated by roots of
      unity of orders coprime to N: the cyclotomic polynomial stays
      irreducible (linear disjointness).
*/
int
_gr_tower_structural_proven(const gr_tower_gen_struct * g, const gr_tower_t U, slong limit)
{
    int ok = 0;

    if (!GR_TOWER_HAS_CAP(U, GR_TOWER_CAP_RATIONAL))
        return 0;

    if (g->def_kind == GR_TOWER_ROOT)
    {
        fmpz_t p, pu, t;
        ulong n, m;
        fmpz_init(p);
        fmpz_init(pu);
        fmpz_init(t);
        if (g->status == GR_TOWER_STATUS_PROVEN && _gr_tower_gen_const_root(g, p, &n) == 2 && fmpz_sgn(p) > 0)
        {
            slong du;
            ok = 1;
            for (du = 0; du < limit && ok; du++)
            {
                const gr_tower_gen_struct * h = GR_TOWER_GEN(U, du);
                if (h->kind != GR_TOWER_ALGEBRAIC)
                    continue;
                /* i does not interfere: the radicals generate a real
                   field R, i is not in R, so [R(i) : Q] = 2 [R : Q] and
                   the degrees over Q(i) are those over Q */
                if (h->def_kind == GR_TOWER_ROOT_OF_UNITY && h->def_param == 4)
                    continue;
                if (h->status != GR_TOWER_STATUS_PROVEN || _gr_tower_gen_const_root(h, pu, &m) != 2 || fmpz_sgn(pu) <= 0)
                    ok = 0;
                else if (fmpz_equal(p, pu))
                    ok = (n_gcd(n, m) == 1);
                else
                {
                    fmpz_gcd(t, p, pu);
                    ok = fmpz_is_one(t);
                }
            }
        }
        fmpz_clear(p);
        fmpz_clear(pu);
        fmpz_clear(t);
    }
    else if (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param > 0)
    {
        /* (also for a composite order N: Q(zeta_N) and Q(zeta_M) with
           gcd(N, M) = 1 are linearly disjoint) */
        slong du;
        ok = 1;
        for (du = 0; du < limit && ok; du++)
        {
            const gr_tower_gen_struct * h = GR_TOWER_GEN(U, du);
            if (h->kind == GR_TOWER_ALGEBRAIC &&
                !(h->def_kind == GR_TOWER_ROOT_OF_UNITY && h->def_param > 0 &&
                  n_gcd(h->def_param, g->def_param) == 1))
                ok = 0;
        }
    }

    /* a root of a prime p over a field unramified at p */
    if (!ok && g->def_kind == GR_TOWER_ROOT && g->status == GR_TOWER_STATUS_PROVEN &&
        _gr_tower_unramified_radical(g, U, limit, 1))
        ok = 1;

    return ok;
}

/*
    Structural irreducibility of step k of T, from the values of the
    generators before it (whatever their statuses), when its modulus is
    the full minimal polynomial over Q of a structured generator:

    - x^n - p, p prime: the algebraic generators before it (degree-one
      steps excepted) are roots of unity of orders prime to p and roots
      of positive integers prime to p of orders prime to p; they
      generate a field unramified at p, over which x^n - p is
      Eisenstein;

    - Phi_N (a root of unity, also of composite order): the generators
      before it are roots of unity of orders prime to N and roots of
      positive integers prime to N of orders prime to N; they generate
      a field unramified at the primes dividing N, linearly disjoint
      from Q(zeta_N) (the intersection is unramified everywhere).

    (i next to zeta_111 and root_37(3), which the modular proofs miss
    since primes = 1 mod 4 split x^2 + 1.) Transcendental generators do
    not matter.
*/
int
_gr_tower_structural_step_proven(gr_tower_t T, slong k)
{
    const gr_tower_gen_struct * g = GR_TOWER_STEP(T, k - 1);
    fmpz_t p, pu;
    ulong n, m, l;
    slong du;
    int kind, ok = 1;

    if (!GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL))
        return 0;

    fmpz_init(p);
    fmpz_init(pu);
    kind = _gr_tower_gen_const_root(g, p, &n);

    if (kind == 1)
    {
        /* (Phi_N for a composite N: the generators before it generate a
           field unramified at the primes dividing N, whose intersection
           with Q(zeta_N) is unramified everywhere, hence Q; l is then the
           radical of N, tested by gcds below) */
        ok = (n >= 1 && gr_tower_step_minpoly(T, k)->length == (slong) n_euler_phi(n) + 1);
        l = n;
    }
    else if (kind == 2)
    {
        ok = (fmpz_abs_fits_ui(p) && n_is_prime(fmpz_get_ui(p)) &&
              gr_tower_step_minpoly(T, k)->length == (slong) n + 1);
        l = ok ? fmpz_get_ui(p) : 0;
    }
    else
        ok = 0;

    for (du = 0; du < g->def_order && ok; du++)
    {
        const gr_tower_gen_struct * h = GR_TOWER_GEN(T, du);
        int hk;

        if (h->kind != GR_TOWER_ALGEBRAIC || gr_tower_step_minpoly(T, h->index)->length <= 2)
            continue;

        hk = _gr_tower_gen_const_root(h, pu, &m);
        if (kind == 1)
        {
            /* (l = N: no prime of N divides m or the radicand) */
            if (hk == 1)
                ok = (n_gcd(m, l) == 1);
            else if (hk == 2)
            {
                fmpz_t gg;
                fmpz_init(gg);
                fmpz_gcd_ui(gg, pu, l);
                ok = (n_gcd(m, l) == 1) && fmpz_is_one(gg);
                fmpz_clear(gg);
            }
            else
                ok = 0;
        }
        else if (hk == 1)
            ok = (m % l != 0);
        else if (hk == 2)
            ok = (m % l != 0) && !fmpz_divisible_ui(pu, l);
        else
            ok = 0;
    }

    fmpz_clear(p);
    fmpz_clear(pu);
    return ok;
}

/* The least definition order of an algebraic generator with a nonreal
   enclosure used by the element img of the top field of U, or -1. */
static slong
_nonreal_position(gr_srcptr img, gr_tower_t U)
{
    gr_tower_flat_struct * F = &U->flat;
    fmpz_mpoly_q_t f;
    int * used;
    slong d, pos = -1, nvars;

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(f, F->mctx);
    if (gr_tower_flat_set_nested_at(f, img, U->length, F) != GR_SUCCESS)
    {
        fmpz_mpoly_q_clear(f, F->mctx);
        return -1;
    }
    nvars = F->mctx->minfo->nvars;
    used = flint_calloc(nvars, sizeof(int));
    fmpz_mpoly_q_used_vars(used, f, F->mctx);
    for (d = 0; d < U->num_gens && pos < 0; d++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(U, d);
        slong v = GR_TOWER_FLAT_VAR_D(F, d);
        if (v < nvars && used[v] && g->kind == GR_TOWER_ALGEBRAIC && !arb_contains_zero(acb_imagref(&g->enclosure)))
            pos = d;
    }
    flint_free(used);
    fmpz_mpoly_q_clear(f, F->mctx);
    return pos;
}

/*
    The real generator a_k of B, found expressible in U only through a
    nonreal generator at definition order pos (sqrt(5) in terms of a
    fifth root of unity): instead of the expression, a_k is adjoined to U
    and moved before pos, when its minimal polynomial q (mapped to U)
    only involves generators before pos; the steps after it become
    dynamic (the fifth root of unity then has a quadratic modulus over
    Q(sqrt 5), found by Trager's method when the norms are small). Real
    subfields are thus represented by real generators. Returns the gid of
    the new generator (to be moved by _absorb_real_first_finish), or -1.
*/
static slong
_absorb_real_first(gr_tower_t U, gr_tower_map_t map, gr_tower_gen_struct * step, const gr_poly_t q, const acb_t z, slong pos)
{
    gr_ctx_struct * top = gr_tower_field(U);
    gr_tower_flat_struct * F = &U->flat;
    slong i, gid;
    int status;

    /* the coefficients of q lie below pos */
    for (i = 0; i < q->length; i++)
    {
        fmpz_mpoly_q_t f;
        slong lev;
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(f, F->mctx);
        status = gr_tower_flat_set_nested_at(f, gr_poly_coeff_srcptr(q, i, top), U->length, F);
        lev = gr_tower_flat_level(f, F);
        fmpz_mpoly_q_clear(f, F->mctx);
        if (status != GR_SUCCESS || lev > pos)
            return -1;
    }

    if (gr_tower_adjoin_algebraic(U, q, z, GR_TOWER_STATUS_DYNAMIC, step->name) != GR_SUCCESS)
        return -1;

    {
        gr_tower_gen_struct * ng = GR_TOWER_STEP(U, U->length - 1);
        ng->def_id = step->def_id;
        if (step->origin != NULL)
            _gr_tower_gen_set_origin(ng, step->origin);
        if (step->def_kind != GR_TOWER_ALGEBRAIC)
        {
            ng->def_kind = step->def_kind;
            ng->def_param = step->def_param;
            if (step->arg.mctx != NULL)
            {
                gr_tower_map_sync(map);
                ng->arg.mctx = map->mctx;
                fmpz_mpoly_q_init(&ng->arg.data, ng->arg.mctx);
                if (gr_tower_map_apply_flat(&ng->arg.data, &step->arg.data, step->arg.mctx, map) != GR_SUCCESS)
                {
                    fmpz_mpoly_q_clear(&ng->arg.data, ng->arg.mctx);
                    ng->arg.mctx = NULL;
                    ng->def_kind = GR_TOWER_ALGEBRAIC;
                }
            }
        }
        gid = ng->gid;
    }

    return gid;
}

/* (after the caller has released its elements of the nested fields of U,
   which the move rebuilds) */
static void
_absorb_real_first_finish(gr_tower_t U, gr_tower_map_t map, slong gid, slong d, slong pos)
{
    slong k;

    _gr_tower_move_gen(U, gr_tower_gid_order(U, gid), pos);

    /* proofs (and refinements) of the moved generator and of the steps
       after it, while the norms are small */
    for (k = 1; k <= U->length; k++)
    {
        gr_tower_gen_struct * g = GR_TOWER_STEP(U, k - 1);
        if (g->def_order >= pos && g->status != GR_TOWER_STATUS_PROVEN)
        {
            if (GR_TOWER_HAS_CAP(U, GR_TOWER_CAP_RATIONAL))
                gr_tower_prove_step_modular(U, k, GR_TOWER_OPTION(U, GR_TOWER_OPT_MODULAR_TRIES));
            if (g->status != GR_TOWER_STATUS_PROVEN)
                gr_tower_prove_step_trager(U, k, GR_TOWER_OPTION(U, GR_TOWER_OPT_TRAGER_DEGREE_LIMIT));
        }
    }

    _map_set_image_gen(map, d, gr_tower_gid_order(U, gid));
}

/*
    When the coefficients of q (over the top field of U) lie in the base
    field (Q or Q(t)), replaces q by its irreducible factor over the base
    vanishing at z (identified numerically), e.g. a modulus which splits
    into linear factors over Q(t). Leaves q unchanged otherwise. (Exact
    factorization over U by Trager's method would also apply over towers
    with transcendental generators, but the multivariate norms are
    costly compared to what it saves.)
*/
static int
_base_factor_root(gr_poly_t q, const acb_t z, gr_tower_t U)
{
    gr_tower_flat_struct * F = &U->flat;
    gr_ctx_struct * top = gr_tower_field(U);
    gr_ctx_struct * B = U->base;
    gr_ctx_t pctx;
    gr_poly_t N;
    gr_vec_t fac;
    fmpz_vec_t mult;
    fmpz_mpoly_q_t t;
    slong i, n = q->length - 1, pick = -1;
    int status = GR_SUCCESS, base = 1;

    if (n < 2 || (B->which_ring != GR_CTX_FMPQ && B->which_ring != GR_CTX_FMPZ_MPOLY_Q))
        return GR_SUCCESS;

    /* (cheap test on the nested representation first: an element of
       F_k is a reduced polynomial over F_{k-1}) */
    for (i = 0; i <= n; i++)
    {
        gr_srcptr x = gr_poly_coeff_srcptr(q, i, top);
        slong k;
        for (k = U->length; k > 0; k--)
        {
            const gr_poly_struct * p = x;
            if (p->length > 1)
                return GR_SUCCESS;
            if (p->length == 0)
                break;
            x = p->coeffs;
        }
    }

    gr_tower_flat_ensure(F);
    gr_poly_init2(N, n + 1, B);
    fmpz_mpoly_q_init(t, F->mctx);
    for (i = 0; i <= n && base && status == GR_SUCCESS; i++)
    {
        status |= gr_tower_flat_set_nested_at(t, gr_poly_coeff_srcptr(q, i, top), U->length, F);
        if (status == GR_SUCCESS &&
            (gr_tower_flat_has_alg_var(fmpz_mpoly_q_numref(t), F) || gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(t), F)))
            base = 0;
        if (status == GR_SUCCESS && base)
            status |= gr_tower_flat_get_nested_at(gr_poly_coeff_ptr(N, i, B), t, 0, F);
    }
    fmpz_mpoly_q_clear(t, F->mctx);
    _gr_poly_set_length(N, n + 1, B);
    _gr_poly_normalise(N, B);

    if (status != GR_SUCCESS || !base || N->length != n + 1)
    {
        gr_poly_clear(N, B);
        return status;
    }

    gr_ctx_init_gr_poly(pctx, B);
    gr_vec_init(fac, 0, pctx);
    fmpz_vec_init(mult, 0);

    if (_gr_tower_base_poly_factor(fac, mult, N, U) == GR_SUCCESS &&
        (fac->length > 1 || (fac->length == 1 && !fmpz_is_one(mult->entries + 0))))
    {
        /* the factor (over the base field) vanishing at z */
        pick = _gr_tower_select_poly_factor(fac, 0, _gr_tower_get_z_fixed, (void *) z, U);

        if (pick >= 0 && status == GR_SUCCESS)
        {
            gr_poly_struct * h = gr_vec_entry_ptr(fac, pick, pctx);
            gr_poly_fit_length(q, h->length, top);
            for (i = 0; i < h->length && status == GR_SUCCESS; i++)
                status |= gr_tower_promote(gr_poly_coeff_ptr(q, i, top), gr_poly_coeff_srcptr(h, i, B), 0, U->length, U);
            _gr_poly_set_length(q, h->length, top);
        }
    }

    fmpz_vec_clear(mult);
    gr_vec_clear(fac, pctx);
    gr_ctx_clear(pctx);
    gr_poly_clear(N, B);
    return status;
}

/*
    An algebraic generator of U with the same origin polynomial as step k
    of B (an integer polynomial both are roots of; when step has none,
    that of a generator of U of which step is shown to be a root) which
    is the same root:
    the roots of the origin polynomial are isolated by disjoint disks, and
    two generators whose enclosures meet the same disk only are the same
    root. This identifies a generator defined in different ways
    (root_4(sqrt(2)) and root_8(2), with sqrt(2) written in terms of a
    root of unity) without exact computations in U. Returns its index k,
    or -1.
*/

static slong
_origin_root_index(acb_srcptr roots, slong deg, const acb_t z)
{
    slong i, idx = -1;

    for (i = 0; i < deg; i++)
    {
        if (acb_overlaps(roots + i, z))
        {
            if (idx >= 0)
                return -1;
            idx = i;
        }
    }

    return idx;
}

/* step k of B (with modulus m) is a root of the integer polynomial P:
   the remainder of P modulo m vanishes (structurally, which is a proof) */
static int
_step_is_root_of(gr_tower_t B, slong k, const fmpz_poly_t P)
{
    const gr_poly_struct * m = gr_tower_step_minpoly(B, k);
    gr_ctx_struct * below = gr_tower_field_at(B, k - 1);
    gr_poly_t p, r;
    slong i;
    int status = GR_SUCCESS, res;

    if (m->length < 2 || fmpz_poly_length(P) < m->length)
        return 0;

    gr_poly_init(p, below);
    gr_poly_init(r, below);
    gr_poly_fit_length(p, fmpz_poly_length(P), below);
    for (i = 0; i < fmpz_poly_length(P) && status == GR_SUCCESS; i++)
        status |= gr_set_fmpz(gr_poly_coeff_ptr(p, i, below), P->coeffs + i, below);
    _gr_poly_set_length(p, fmpz_poly_length(P), below);
    if (status == GR_SUCCESS)
        status = gr_poly_rem(r, p, m, below);
    res = (status == GR_SUCCESS && gr_poly_is_zero(r, below) == T_TRUE);
    gr_poly_clear(p, below);
    gr_poly_clear(r, below);
    return res;
}

static slong
_absorb_same_origin_root(gr_tower_t U, gr_tower_t B, slong kb)
{
    const gr_tower_gen_struct * step = GR_TOWER_STEP(B, kb - 1);
    const fmpz_poly_struct * origin = step->origin;
    slong k, deg, prec, res = -1, si;
    acb_ptr roots;
    acb_t w;

    if (origin == NULL)
    {
        /* the origin polynomial of a generator of U which may be the
           same root, if step is a root of it */
        for (k = 1; k <= U->length && origin == NULL; k++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_STEP(U, k - 1);
            if (g->origin != NULL && fmpz_poly_degree(g->origin) <= GR_TOWER_OPTION(U, GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT) &&
                acb_overlaps(&g->enclosure, &step->enclosure) &&
                _step_is_root_of(B, kb, g->origin))
                origin = g->origin;
        }
        if (origin == NULL)
            return -1;
    }
    else
    {
        int any = 0;
        for (k = 1; k <= U->length && !any; k++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_STEP(U, k - 1);
            any = (g->origin != NULL && fmpz_poly_equal(g->origin, origin) &&
                   acb_overlaps(&g->enclosure, &step->enclosure));
        }
        if (!any)
            return -1;
    }

    if (fmpz_poly_degree(origin) < 2 || fmpz_poly_degree(origin) > GR_TOWER_OPTION(U, GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT) ||
        !fmpz_poly_is_squarefree(origin))
        return -1;

    deg = fmpz_poly_degree(origin);
    roots = _acb_vec_init(deg);
    acb_init(w);

    prec = FLINT_MAX(GR_TOWER_DEFAULT_PREC, acb_rel_accuracy_bits(&step->enclosure));
    arb_fmpz_poly_complex_roots(roots, origin, 0, prec);

    si = _origin_root_index(roots, deg, &step->enclosure);

    for (k = 1; k <= U->length && res < 0 && si >= 0; k++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_STEP(U, k - 1);

        if (g->origin == NULL || !fmpz_poly_equal(g->origin, origin) ||
            !acb_overlaps(&g->enclosure, &step->enclosure))
            continue;

        if (gr_tower_step_get_acb(w, U, k, prec) == GR_SUCCESS &&
            _origin_root_index(roots, deg, w) == si)
            res = k;
    }

    _acb_vec_clear(roots, deg);
    acb_clear(w);
    return res;
}

static int
_absorb_alg(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, slong k, slong d, int flags)
{
    gr_ctx_struct * top;
    gr_tower_gen_struct * step = GR_TOWER_STEP(B, k - 1);
    gr_poly_t q;
    acb_t z;
    int found = 0, search = 1;
    int status = GR_SUCCESS;

    /* the same definition already present in U (as any kind of generator) */
    if (step->def_id != 0)
    {
        slong du = gr_tower_find_def_order(U, step->def_id);
        if (du >= 0)
        {
            _map_set_image_gen(map, d, du);
            return GR_SUCCESS;
        }
    }

    /* an exp, log or special function value which became algebraic:
       the same function of the same argument already in U */
    if (step->arg.mctx != NULL && step->def_kind != GR_TOWER_ROOT &&
        step->def_kind != GR_TOWER_ROOT_OF_UNITY && GR_TOWER_KIND_HAS_ARG(step->def_kind))
    {
        fmpz_mpoly_q_t arg;
        fmpz_mpoly_ctx_struct * actx;
        slong du;

        gr_tower_map_sync(map);
        actx = map->mctx;
        fmpz_mpoly_q_init(arg, actx);
        if (gr_tower_map_apply_flat(arg, &step->arg.data, step->arg.mctx, map) == GR_SUCCESS)
            du = _absorb_find_definition(U, map, step->def_kind, step->def_param, arg, &actx);
        else
            du = -1;
        fmpz_mpoly_q_clear(arg, actx);

        if (du >= 0)
        {
            _map_set_image_gen(map, d, du);
            return GR_SUCCESS;
        }
    }

    if (_absorb_structured(U, map, step, d, !(flags & GR_TOWER_MERGE_REAL_FIRST), GR_TOWER_OPTION(U, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT) > 0))
        return GR_SUCCESS;

    {
        slong ku = _absorb_same_origin_root(U, B, k);
        if (ku > 0)
        {
            _map_set_image_gen(map, d, GR_TOWER_STEP(U, ku - 1)->def_order);
            return GR_SUCCESS;
        }
    }

    if (_absorb_same_root(U, map, step, d))
        return GR_SUCCESS;

    top = gr_tower_field(U);
    gr_poly_init(q, top);
    acb_init(z);

    /* the minimal polynomial of a_k over F_{k-1} of B, mapped to U
       (not needed when the origin polynomial replaces it below) */
    if (!(step->origin != NULL && gr_tower_step_minpoly(B, k)->length >= 3))
    {
        const gr_poly_struct * m = gr_tower_step_minpoly(B, k);
        gr_ctx_struct * below = gr_tower_field_at(B, k - 1);
        slong i;

        gr_poly_fit_length(q, m->length, top);
        for (i = 0; i < m->length && status == GR_SUCCESS; i++)
            status |= gr_tower_map_apply_at(gr_poly_coeff_ptr(q, i, top), gr_poly_coeff_srcptr(m, i, below), k - 1, map);
        _gr_poly_set_length(q, m->length, top);
        _gr_poly_normalise(q, top);
    }
    acb_set(z, &step->enclosure);

    /* conjugates already in U: every algebraic generator of U which is
       a root of q is divided out (q(beta) = 0 is filtered numerically,
       then decided exactly), so that conjugate roots of one polynomial
       build up a splitting tower with the Cauchy moduli
       p(x) / prod (x - beta) rather than towers of reducible moduli.
       With a known origin polynomial, the division starts from it (the
       modulus in B may already be a Cauchy factor with respect to other
       conjugates), and the conjugates in U are recognized by it. */
    if (status == GR_SUCCESS && step->origin != NULL && gr_tower_step_minpoly(B, k)->length >= 3)
    {
        /* made monic over QQ (an inversion in the top field would be
           expensive) */
        fmpq_poly_t oq;
        fmpq_t c;
        slong i;
        fmpq_poly_init(oq);
        fmpq_init(c);
        fmpq_poly_set_fmpz_poly(oq, step->origin);
        fmpq_poly_make_monic(oq, oq);
        gr_poly_fit_length(q, oq->length, top);
        for (i = 0; i < oq->length && status == GR_SUCCESS; i++)
        {
            fmpq_poly_get_coeff_fmpq(c, oq, i);
            status |= gr_set_fmpq(gr_poly_coeff_ptr(q, i, top), c, top);
        }
        _gr_poly_set_length(q, oq->length, top);
        _gr_poly_normalise(q, top);
        fmpq_clear(c);
        fmpq_poly_clear(oq);
    }

    if (status == GR_SUCCESS && q->length >= 3)
    {
        slong match;
        status = _divide_out_conjugates(q, z, &match, step->origin, U);
        if (status == GR_SUCCESS && match > 0)
        {
            /* the root is a generator of U */
            _map_set_image_gen(map, d, GR_TOWER_STEP(U, match - 1)->def_order);
            gr_poly_clear(q, top);
            acb_clear(z);
            return GR_SUCCESS;
        }
        if (status == GR_SUCCESS && step->origin != NULL && q->length < step->origin->length)
        {
            /* conjugates were divided out: the step is pushed with the
               Cauchy factor q as its (only conjecturally irreducible)
               modulus. The enclosure of the step isolates a root of its
               modulus in B, which need not isolate it among the roots of
               q (a conjugate divided out in B but not in U may be close),
               so the root is certified against q. */
            if (q->length >= 3 && gr_tower_adjoin_algebraic(U, q, z, GR_TOWER_STATUS_DYNAMIC, step->name) == GR_SUCCESS)
            {
                gr_tower_gen_struct * ng = GR_TOWER_STEP(U, U->length - 1);
                ng->def_id = step->def_id;
                _gr_tower_gen_set_origin(ng, step->origin);
                _map_set_image_gen(map, d, ng->def_order);
                gr_poly_clear(q, top);
                acb_clear(z);
                return GR_SUCCESS;
            }
        }
    }

    /* a modulus over the base field of U which factors there */
    if (status == GR_SUCCESS && q->length >= 3)
        status = _base_factor_root(q, z, U);

    if (status == GR_SUCCESS)
    {
        gr_ptr img;
        GR_TMP_INIT(img, top);

        if (q->length == 2)
        {
            /* linear: the generator is -q[0] */
            status |= gr_neg(img, gr_poly_coeff_srcptr(q, 0, top), top);
            found = 1;
        }
        else if ((flags & GR_TOWER_MERGE_EXPRESS) && q->length == 3)
        {
            /* quadratic: exact search through the quadratic steps of U;
               a negative answer is a proof, so the lattice search is
               skipped then */
            int st = _quadratic_root(img, q, z, U);
            found = (st == GR_SUCCESS);
            if (st == GR_DOMAIN)
                search = 0;
        }

        /* a structurally irreducible modulus (roots of primes over fields
           unramified at them, independent real radicals, cyclotomic
           steps) has no root in U */
        if (!found && search && q->length >= 3 &&
            ((step->status == GR_TOWER_STATUS_PROVEN &&
              ((step->def_kind == GR_TOWER_ROOT && q->length - 1 == step->def_param) ||
               (step->def_kind == GR_TOWER_ROOT_OF_UNITY && step->def_param > 0 &&
                  q->length - 1 == (slong) n_euler_phi(step->def_param))) &&
              _gr_tower_structural_proven(step, U, U->num_gens)) ||
             (step->def_kind == GR_TOWER_ROOT && _gr_tower_unramified_radical(step, U, U->num_gens, 0))))
        {
            /* (the numbers generating U lie in a field unramified at p,
               whatever the statuses of the steps, and x^n - p has no root
               there; for a proven step the modulus is irreducible) */
            search = 0;
        }

        /* (a place of U at which q has no root excludes roots in U) */
        if (!found && search && (flags & GR_TOWER_MERGE_EXPRESS) && gr_tower_degree(U) <= GR_TOWER_OPTION(U, GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT) &&
            !gr_tower_poly_no_roots_modular(q, U, GR_TOWER_OPTION(U, GR_TOWER_OPT_NO_ROOTS_TRIES)))
        {
            found = (gr_tower_express_limit(img, q, z,
                        GR_TOWER_MERGE_EXPRESS_PREC(U, gr_tower_degree(U)), U) == GR_SUCCESS);
        }

        /* a real generator expressed through nonreal ones: real first */
        if (found && (flags & GR_TOWER_MERGE_REAL_FIRST) && q->length >= 3 &&
            arb_is_zero(acb_imagref(&step->enclosure)))
        {
            slong pos = _nonreal_position(img, U), gid = -1;
            if (pos >= 0)
                gid = _absorb_real_first(U, map, step, q, z, pos);
            if (gid >= 0)
            {
                GR_TMP_CLEAR(img, top);
                gr_poly_clear(q, top);
                acb_clear(z);
                _absorb_real_first_finish(U, map, gid, d, pos);
                return GR_SUCCESS;
            }
        }

        if (found)
            status |= gr_tower_map_set_image(map, d, img, U->length);

        GR_TMP_CLEAR(img, top);

        if (!found)
        {
            /* the proof of irreducibility in B carries over when the
               field it was over maps isomorphically into U: with U
               having no algebraic step, and no conjecturally
               transcendental generator on either side (a conjectural
               generator of B may be identified non-injectively, e.g.
               exp(2) with exp(1)^2, under which an irreducible modulus
               can become reducible) */
            int st = (step->status == GR_TOWER_STATUS_PROVEN && U->length == 0 &&
                      !_gr_tower_has_conjectural_below(U, U->num_gens) &&
                      !_gr_tower_has_conjectural_below(B, step->def_order))
                        ? GR_TOWER_STATUS_PROVEN : GR_TOWER_STATUS_DYNAMIC;

            /* structural irreducibility rules (Besicovitch, cyclotomic) */
            if (st == GR_TOWER_STATUS_DYNAMIC && step->status == GR_TOWER_STATUS_PROVEN &&
                _gr_tower_structural_proven(step, U, U->num_gens))
                st = GR_TOWER_STATUS_PROVEN;

            /* a radical X^n - a: Capelli's criterion, certified at a
               place of U (root_53(2) over Q(i, zeta_53), say, which is
               ramified at 2) */
            if (st == GR_TOWER_STATUS_DYNAMIC && step->def_kind == GR_TOWER_ROOT && step->arg.mctx != NULL &&
                step->def_param >= 2 && q->length == step->def_param + 1)
            {
                fmpz_mpoly_q_t a;
                gr_tower_map_sync(map);
                fmpz_mpoly_q_init(a, map->mctx);
                if (gr_tower_map_apply_flat(a, &step->arg.data, step->arg.mctx, map) == GR_SUCCESS &&
                    gr_tower_binomial_irreducible_modular(U, a, map->mctx, step->def_param, GR_TOWER_OPTION(U, GR_TOWER_OPT_MODULAR_TRIES)))
                    st = GR_TOWER_STATUS_PROVEN;
                fmpz_mpoly_q_clear(a, map->mctx);
            }

            status |= gr_tower_adjoin_algebraic(U, q, z, st, step->name);

            if (status == GR_SUCCESS)
            {
                gr_tower_gen_struct * ng = GR_TOWER_STEP(U, U->length - 1);

                ng->def_id = step->def_id;
                if (step->origin != NULL)
                    _gr_tower_gen_set_origin(ng, step->origin);

                /* keep the definition of a generator which became algebraic */
                if (step->def_kind != GR_TOWER_ALGEBRAIC)
                {
                    ng->def_kind = step->def_kind;
                    ng->def_param = step->def_param;
                    if (step->arg.mctx != NULL)
                    {
                        gr_tower_map_sync(map);
                        ng->arg.mctx = map->mctx;
                        fmpz_mpoly_q_init(&ng->arg.data, ng->arg.mctx);
                        if (gr_tower_map_apply_flat(&ng->arg.data, &step->arg.data, step->arg.mctx, map) != GR_SUCCESS)
                        {
                            fmpz_mpoly_q_clear(&ng->arg.data, ng->arg.mctx);
                            ng->arg.mctx = NULL;
                            ng->def_kind = GR_TOWER_ALGEBRAIC;
                        }
                    }
                }

                _map_set_image_gen(map, d, ng->def_order);
            }
        }
    }

    gr_poly_clear(q, top);
    acb_clear(z);
    return status;
}

/*
    Processes the generators of B (in definition order) numbered from
    map->length to p - 1, so that a map can be extended after its source
    has grown.
*/
int
gr_tower_absorb_prefix(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, slong p, int flags)
{
    slong d;
    int status = GR_SUCCESS;

    if (map->target != U || map->source != B)
        return GR_DOMAIN;

    if (p < 0 || p > B->num_gens)
        p = B->num_gens;

    for (d = map->length; d < p && status == GR_SUCCESS; d++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(B, d);

        if (g->kind == GR_TOWER_ALGEBRAIC)
            status = _absorb_alg(U, map, B, g->index, d, flags);
        else
            status = _absorb_trans(U, map, B, g->index, d);

        if (status == GR_SUCCESS)
            map->length = d + 1;
    }

    return status;
}

int
gr_tower_absorb(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, int flags)
{
    return gr_tower_absorb_prefix(U, map, B, -1, flags);
}

int
gr_tower_merge(gr_tower_t U, gr_tower_map_t mapA, gr_tower_map_t mapB, gr_tower_t A, gr_tower_t B, int flags)
{
    int status = GR_SUCCESS;

    gr_tower_set(U, A);
    gr_tower_map_init(mapA, A, U);
    status |= gr_tower_map_set_inclusion(mapA);
    gr_tower_map_init(mapB, B, U);
    status |= gr_tower_absorb(U, mapB, B, flags);
    gr_tower_map_sync(mapA);

    return status;
}

int
gr_tower_eliminate(gr_tower_t U, gr_tower_map_t map, gr_tower_t T)
{
    int status;

    gr_tower_clear(U);
    gr_tower_init(U, T->consts);
    U->keep_retired = 0;   /* (internal: no nested elements across rebuilds) */
    U->options = T->options;
    gr_tower_map_init(map, T, U);
    status = gr_tower_absorb(U, map, T, 0);
    return status;
}
