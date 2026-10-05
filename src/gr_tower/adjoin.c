/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "qqbar.h"
#include "fmpz_mpoly_q.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

int
gr_tower_adjoin_algebraic(gr_tower_t T, const gr_poly_t m, const acb_t z, int status, const char * name)
{
    gr_ctx_struct * top = gr_tower_field(T);
    acb_t zz;
    slong prec;
    int st;

    if (m->length < 2)
        return GR_DOMAIN;

    if (gr_is_one(gr_poly_coeff_srcptr(m, m->length - 1, top), top) != T_TRUE)
        return GR_DOMAIN;

    acb_init(zz);

    prec = FLINT_MAX(GR_TOWER_DEFAULT_PREC, acb_rel_accuracy_bits(z));
    if (prec > 100000)
        prec = 100000;

    st = _gr_tower_certify_root(zz, m, T->length, z, prec, T);

    if (st == GR_SUCCESS)
        _gr_tower_push_step(T, m, zz, prec, status, name);

    acb_clear(zz);
    return st;
}

int
gr_tower_adjoin_qqbar(gr_tower_t T, const qqbar_t x, const char * name)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_poly_t m;
    acb_t z;
    int status;

    gr_poly_init(m, top);
    acb_init(z);

    status = gr_poly_set_fmpz_poly(m, QQBAR_POLY(x), top);
    status |= gr_poly_make_monic(m, m, top);
    /* an enclosure isolating x among the roots: qqbar's own enclosure is
       isolating; a refinement at a precision above its accuracy stays
       isolating also for clustered roots (the two roots of
       x^20 - 2 (101 x - 1)^2 near 1/101 differ by 1.3e-22) */
    {
        slong acc = acb_rel_accuracy_bits(QQBAR_ENCLOSURE(x));
        qqbar_get_acb(z, x, FLINT_MAX(GR_TOWER_DEFAULT_PREC, 2 * FLINT_MAX(acc, 0) + 64));
        if (!acb_contains(QQBAR_ENCLOSURE(x), z))
            acb_set(z, QQBAR_ENCLOSURE(x));
    }

    if (status == GR_SUCCESS)
    {
        /* the minimal polynomial over QQ is irreducible over QQ, but
           may factor over a nontrivial tower */
        int step_status = (T->length == 0 && GR_TOWER_BASE_IS_CONSTS(T) && GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL))
                        ? GR_TOWER_STATUS_PROVEN : GR_TOWER_STATUS_DYNAMIC;

        status = gr_tower_adjoin_algebraic(T, m, z, step_status, name);
        if (status == GR_UNABLE && !acb_equal(z, QQBAR_ENCLOSURE(x)))
            status = gr_tower_adjoin_algebraic(T, m, QQBAR_ENCLOSURE(x), step_status, name);

        if (status == GR_SUCCESS)
            _gr_tower_gen_set_origin(GR_TOWER_STEP(T, T->length - 1), QQBAR_POLY(x));
    }

    gr_poly_clear(m, top);
    acb_clear(z);
    return status;
}

/* (zx: an enclosure of x, or NULL; the caller may know that x is real
   when its enclosure straddles the real axis, which decides the
   principal root of a negative x) */
int
_gr_tower_adjoin_root_ui_enclosure(gr_tower_t T, gr_srcptr x, ulong n, const acb_t zx, const char * name)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_poly_t m;
    acb_t z;
    int status;

    if (n == 0)
        return GR_DOMAIN;

    gr_poly_init(m, top);
    acb_init(z);

    /* m = X^n - x */
    status = gr_poly_set_coeff_si(m, n, 1, top);
    status |= gr_neg(gr_poly_coeff_ptr(m, 0, top), x, top);
    _gr_poly_normalise(m, top);

    if (status == GR_SUCCESS)
        status = _gr_tower_principal_root_enclosure(z, x, n, zx, T);

    if (status == GR_SUCCESS)
    {
        fmpz_mpoly_q_t arg;
        fmpz_mpoly_ctx_struct * actx;

        /* the definition (for identification of repeated roots) */
        gr_tower_flat_ensure(&T->flat);
        actx = T->flat.mctx;
        fmpz_mpoly_q_init(arg, actx);
        status = gr_tower_flat_set_nested_at(arg, x, T->length, &T->flat);

        if (status == GR_SUCCESS)
        {
            /* Capelli: X^n - x is irreducible when x is not a p-th power
               for any prime p | n (and not in -4 K^4 when 4 | n), which
               a place of the tower can certify */
            int st = gr_tower_binomial_irreducible_modular(T, arg, actx, n, GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES))
                        ? GR_TOWER_STATUS_PROVEN : GR_TOWER_STATUS_DYNAMIC;
            status = gr_tower_adjoin_algebraic(T, m, z, st, name);
        }

        if (status == GR_SUCCESS)
        {
            gr_tower_gen_struct * g = GR_TOWER_STEP(T, T->length - 1);
            g->def_kind = GR_TOWER_ROOT;
            g->def_param = n;
            g->arg.mctx = actx;
            fmpz_mpoly_q_init(&g->arg.data, actx);
            fmpz_mpoly_q_set(&g->arg.data, arg, actx);
        }

        fmpz_mpoly_q_clear(arg, actx);
    }

    gr_poly_clear(m, top);
    acb_clear(z);
    return status;
}

int
gr_tower_adjoin_root_ui(gr_tower_t T, gr_srcptr x, ulong n, const char * name)
{
    return _gr_tower_adjoin_root_ui_enclosure(T, x, n, NULL, name);
}

int
gr_tower_adjoin_root_of_unity(gr_tower_t T, ulong n, const char * name)
{
    qqbar_t z;
    int status;

    if (n == 0)
        return GR_DOMAIN;

    qqbar_init(z);
    qqbar_root_of_unity(z, 1, n);
    status = gr_tower_adjoin_qqbar(T, z, name);
    if (status == GR_SUCCESS)
    {
        gr_tower_gen_struct * g = GR_TOWER_STEP(T, T->length - 1);
        g->def_kind = GR_TOWER_ROOT_OF_UNITY;
        g->def_param = n;
    }
    qqbar_clear(z);
    return status;
}

int
gr_tower_adjoin_root_fmpz(gr_tower_t T, const fmpz_t p, ulong n, const char * name)
{
    qqbar_t z;
    int status;

    if (n == 0 || fmpz_is_zero(p))
        return GR_DOMAIN;

    int proven = 0;

    qqbar_init(z);
    qqbar_set_fmpz(z, p);
    qqbar_root_ui(z, z, n);

    /* over a nontrivial tower, X^n - p is irreducible when a place of
       the tower certifies it (Capelli); otherwise the step is dynamic */
    if (T->length > 0 && qqbar_degree(z) == (slong) n && n > 1)
    {
        fmpz_mpoly_q_t a;
        gr_tower_flat_ensure(&T->flat);
        fmpz_mpoly_q_init(a, T->flat.mctx);
        fmpz_mpoly_q_set_fmpz(a, p, T->flat.mctx);
        proven = gr_tower_binomial_irreducible_modular(T, a, T->flat.mctx, n, GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES));
        fmpz_mpoly_q_clear(a, T->flat.mctx);
    }

    status = gr_tower_adjoin_qqbar(T, z, name);
    if (status == GR_SUCCESS)
    {
        gr_tower_gen_struct * g = GR_TOWER_STEP(T, T->length - 1);
        if (proven)
            _gr_tower_gen_set_status(g, GR_TOWER_STATUS_PROVEN);
        gr_tower_flat_ensure(&T->flat);
        g->def_kind = GR_TOWER_ROOT;
        g->def_param = n;
        g->arg.mctx = T->flat.mctx;
        fmpz_mpoly_q_init(&g->arg.data, g->arg.mctx);
        fmpz_mpoly_q_set_fmpz(&g->arg.data, p, g->arg.mctx);
    }
    qqbar_clear(z);
    return status;
}
