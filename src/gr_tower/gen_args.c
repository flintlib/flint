/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Arguments of generators. A transcendental generator (and a generator
    which became algebraic, keeping its definition) has the argument arg
    of its definition, and functions of several arguments (pFq) the
    additional arguments xargs. They are flat elements of the tower, each
    in the context it was stored in (converted on use). The argument with
    index 0 is arg, the argument with index i >= 1 is xargs[i - 1].
*/

#include <string.h>
#include "acb.h"
#include "fmpz_mpoly_q.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

slong
_gr_tower_gen_num_args(const gr_tower_gen_struct * g)
{
    return (g->arg.mctx != NULL) + g->num_xargs;
}

gr_tower_flat_elem_struct *
_gr_tower_gen_arg_ptr(const gr_tower_gen_struct * g, slong i)
{
    return (i == 0) ? (gr_tower_flat_elem_struct *) &g->arg : g->xargs + i - 1;
}

void
_gr_tower_gen_xargs_clear(gr_tower_gen_struct * g)
{
    slong i;

    for (i = 0; i < g->num_xargs; i++)
        if (g->xargs[i].mctx != NULL)
            fmpz_mpoly_q_clear(&g->xargs[i].data, g->xargs[i].mctx);

    flint_free(g->xargs);
    g->xargs = NULL;
    g->num_xargs = 0;
}

/*
    The arguments as a string, in the order in which the function prints
    them: for pFq "a_1, ..., a_p, b_1, ..., b_q, z" (as in hypgeom_2f1(a,
    b, c, z)) for p <= 2, q = 1, otherwise "[a_1, ...], [b_1, ...], z";
    the argument alone for the other kinds.
*/
char *
_gr_tower_gen_args_str(const gr_tower_gen_struct * g, gr_tower_t T)
{
    char ** parts;
    char * s;
    slong i, n, len, p = 0, q = 0;
    int brackets = 0;

    if (g->arg.mctx == NULL)
        return NULL;

    if (g->num_xargs == 0)
        return _gr_tower_flat_get_str(&g->arg.data, g->arg.mctx, T);

    n = 1 + g->num_xargs;
    parts = flint_malloc(sizeof(char *) * n);
    for (i = 0; i < n; i++)
    {
        const gr_tower_flat_elem_struct * a = _gr_tower_gen_arg_ptr(g, i);
        parts[i] = _gr_tower_flat_get_str(&a->data, a->mctx, T);
    }

    if (g->def_kind == GR_TOWER_JACOBI_THETA || g->def_kind == GR_TOWER_HURWITZ_ZETA)
    {
        s = flint_malloc(strlen(parts[0]) + strlen(parts[1]) + 4);
        strcpy(s, parts[0]);
        strcat(s, ", ");
        strcat(s, parts[1]);
        for (i = 0; i < n; i++)
            flint_free(parts[i]);
        flint_free(parts);
        return s;
    }

    if (g->def_kind == GR_TOWER_HYPGEOM)
    {
        p = GR_TOWER_HYPGEOM_P(g->def_param);
        q = GR_TOWER_HYPGEOM_Q(g->def_param);
        brackets = !(p <= 2 && q == 1);
    }

    len = 16;
    for (i = 0; i < n; i++)
        len += strlen(parts[i]) + 2;

    s = flint_malloc(len);
    s[0] = '\0';

    /* the parameters (arguments 1, ..., n - 1), then z (argument 0) */
    if (brackets)
        strcat(s, "[");
    for (i = 1; i < n; i++)
    {
        if (brackets && i == p + 1)
            strcat(s, "], [");
        else if (i > 1)
            strcat(s, ", ");
        strcat(s, parts[i]);
    }
    if (brackets)
        strcat(s, (q == 0 && p == n - 1) ? "], [], " : "], ");
    else
        strcat(s, ", ");
    strcat(s, parts[0]);

    for (i = 0; i < n; i++)
        flint_free(parts[i]);
    flint_free(parts);
    return s;
}

/* Enclosures of all the arguments (in the order arg, xargs). */
int
_gr_tower_gen_args_get_acb(acb_ptr res, const gr_tower_gen_struct * g, slong prec, gr_tower_flat_struct * F)
{
    slong i, n = _gr_tower_gen_num_args(g);
    int status = GR_SUCCESS;

    gr_tower_flat_ensure(F);

    for (i = 0; i < n && status == GR_SUCCESS; i++)
    {
        const gr_tower_flat_elem_struct * a = _gr_tower_gen_arg_ptr(g, i);
        fmpz_mpoly_q_t u;
        fmpz_mpoly_q_init(u, F->mctx);
        gr_tower_flat_convert(u, &a->data, a->mctx, F);
        status = gr_tower_flat_get_acb(res + i, u, prec, F);
        fmpz_mpoly_q_clear(u, F->mctx);
    }

    return status;
}

/*
    Copies the definition (def_kind, def_param and the arguments) of the
    generator g of the source of the map to the generator ng of its
    target, applying the map to the arguments. If an argument cannot be
    mapped, ng is left as a plain algebraic generator; returns 0 then.
*/
int
_gr_tower_gen_copy_def_map(gr_tower_gen_struct * ng, const gr_tower_gen_struct * g, gr_tower_map_t map)
{
    slong i, n, done;

    ng->def_kind = g->def_kind;
    ng->def_param = g->def_param;

    if (g->arg.mctx == NULL)
        return 1;

    n = _gr_tower_gen_num_args(g);

    if (g->num_xargs > 0)
    {
        ng->num_xargs = g->num_xargs;
        ng->xargs = flint_malloc(sizeof(gr_tower_flat_elem_struct) * g->num_xargs);
        for (i = 0; i < g->num_xargs; i++)
            ng->xargs[i].mctx = NULL;
    }

    for (done = 0; done < n; done++)
    {
        const gr_tower_flat_elem_struct * a = _gr_tower_gen_arg_ptr(g, done);
        gr_tower_flat_elem_struct * b = _gr_tower_gen_arg_ptr(ng, done);
        fmpz_mpoly_q_t t;
        fmpz_mpoly_ctx_struct * mctx;
        int status;

        gr_tower_map_sync(map);
        mctx = map->mctx;
        fmpz_mpoly_q_init(t, mctx);
        status = gr_tower_map_apply_flat(t, &a->data, a->mctx, map);
        if (status != GR_SUCCESS || map->mctx != mctx)
        {
            fmpz_mpoly_q_clear(t, mctx);
            break;
        }
        b->mctx = mctx;
        fmpz_mpoly_q_init(&b->data, mctx);
        fmpz_mpoly_q_swap(&b->data, t, mctx);
        fmpz_mpoly_q_clear(t, mctx);
    }

    if (done < n)
    {
        if (ng->arg.mctx != NULL)
            fmpz_mpoly_q_clear(&ng->arg.data, ng->arg.mctx);
        ng->arg.mctx = NULL;
        _gr_tower_gen_xargs_clear(ng);
        ng->def_kind = GR_TOWER_ALGEBRAIC;
        return 0;
    }

    return 1;
}

/* the flags of acb_hypgeom_2f1 for a 2F1 generator (0 otherwise; see
   _gr_tower_hypgeom_2f1_flags), from its exact arguments */
int
_gr_tower_gen_hypgeom_flags(const gr_tower_gen_struct * g, gr_tower_flat_struct * F)
{
    fmpz_mpoly_q_struct u[3];
    slong i;
    int flags;

    if (g->def_kind != GR_TOWER_HYPGEOM || GR_TOWER_HYPGEOM_P(g->def_param) != 2 ||
        GR_TOWER_HYPGEOM_Q(g->def_param) != 1 || _gr_tower_gen_num_args(g) != 4)
        return 0;

    gr_tower_flat_ensure(F);
    for (i = 0; i < 3; i++)
    {
        const gr_tower_flat_elem_struct * a = _gr_tower_gen_arg_ptr(g, i + 1);
        fmpz_mpoly_q_init(u + i, F->mctx);
        gr_tower_flat_convert(u + i, &a->data, a->mctx, F);
    }
    flags = _gr_tower_hypgeom_2f1_flags(u + 0, u + 1, u + 2, F->mctx);
    for (i = 0; i < 3; i++)
        fmpz_mpoly_q_clear(u + i, F->mctx);
    return flags;
}

POP_OPTIONS
