/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* Generic-ring counterparts of acb_dirichlet_dft_index and
   acb_dirichlet_dft: the group of Dirichlet characters mod q is the
   product of the cyclic component groups given by the Conrey
   decomposition, so its DFT is a product DFT over the component
   sizes, with a gather/scatter between number indexing and
   lexicographic Conrey indexing for the non-index variant. */

static void
_gr_dft_dirichlet_cyc(ulong * cyc, const dirichlet_group_t G)
{
    slong k;

    for (k = 0; k < G->num; k++)
        cyc[k] = G->P[k].phi.n;
}

/* dft, lexicographic conrey indexing, array size G->phi_q */
int
gr_dft_dirichlet_index(gr_ptr w, gr_srcptr v, const dirichlet_group_t G,
        gr_ctx_t ctx)
{
    int status;
    ulong * cyc;

    if (G->phi_q == 1)
        return gr_set(w, v, ctx);

    cyc = flint_malloc(G->num * sizeof(ulong));
    _gr_dft_dirichlet_cyc(cyc, G);
    status = gr_dft_prod(w, v, cyc, G->num, ctx);
    flint_free(cyc);

    return status;
}

/* dft, number indexing, array size G->q */
int
gr_dft_dirichlet(gr_ptr w, gr_srcptr v, const dirichlet_group_t G,
        gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong i, len = G->phi_q;
    gr_ptr t;
    dirichlet_char_t x;

    /* the trivial group: the Conrey number of the unique character
       is 1, outside an array of size q = 1, so the generic
       gather/scatter below must not run */
    if (G->q == 1)
        return gr_set(w, v, ctx);

    t = gr_heap_init_vec(len, ctx);
    dirichlet_char_init(x, G);

    dirichlet_char_one(x, G);
    for (i = 0; i < len; i++)
    {
        status |= gr_set(GR_ENTRY(t, i, sz), GR_ENTRY(v, x->n, sz), ctx);
        dirichlet_char_next(x, G);
    }

    status |= gr_dft_dirichlet_index(t, t, G, ctx);

    dirichlet_char_one(x, G);
    for (i = 0; i < len; i++)
    {
        status |= gr_set(GR_ENTRY(w, x->n, sz), GR_ENTRY(t, i, sz), ctx);
        dirichlet_char_next(x, G);
    }

    dirichlet_char_clear(x);
    gr_heap_clear_vec(t, len, ctx);

    return status;
}

/* acb variants routed through the fixed-point product transform */
void
gr_dft_acb_dirichlet_index(acb_ptr w, acb_srcptr v,
        const dirichlet_group_t G, slong prec)
{
    ulong * cyc;

    if (G->phi_q == 1)
    {
        acb_set(w, v);
        return;
    }

    cyc = flint_malloc(G->num * sizeof(ulong));
    _gr_dft_dirichlet_cyc(cyc, G);
    gr_dft_acb_prod(w, v, cyc, G->num, prec);
    flint_free(cyc);
}

void
gr_dft_acb_dirichlet(acb_ptr w, acb_srcptr v, const dirichlet_group_t G,
        slong prec)
{
    ulong i, len = G->phi_q;
    acb_ptr t;
    dirichlet_char_t x;

    if (G->q == 1)
    {
        acb_set(w, v);
        return;
    }

    t = _acb_vec_init(len);
    dirichlet_char_init(x, G);

    dirichlet_char_one(x, G);
    for (i = 0; i < len; i++)
    {
        acb_set(t + i, v + x->n);
        dirichlet_char_next(x, G);
    }

    gr_dft_acb_dirichlet_index(t, t, G, prec);

    dirichlet_char_one(x, G);
    for (i = 0; i < len; i++)
    {
        acb_set(w + x->n, t + i);
        dirichlet_char_next(x, G);
    }

    dirichlet_char_clear(x);
    _acb_vec_clear(t, len);
}
