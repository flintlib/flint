/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_tower.h"

/* Heap-allocated towers and non-inline accessors, for language bindings. */

gr_tower_struct *
gr_tower_heap_init(gr_ctx_t base)
{
    gr_tower_struct * T = flint_malloc(sizeof(gr_tower_struct));
    gr_tower_init(T, base);
    return T;
}

void
gr_tower_heap_clear(gr_tower_struct * T)
{
    gr_tower_clear(T);
    flint_free(T);
}

int
gr_tower_get_str(char ** s, const gr_tower_t T)
{
    gr_stream_t out;
    int status;
    gr_stream_init_str(out);
    status = gr_tower_write(out, T);
    *s = out->s;
    return status;
}

gr_ctx_struct * gr_tower_field_ptr(const gr_tower_t T) { return gr_tower_field(T); }
slong gr_tower_length_si(const gr_tower_t T) { return gr_tower_length(T); }
slong gr_tower_num_gens_si(const gr_tower_t T) { return gr_tower_num_gens(T); }
slong gr_tower_num_trans_si(const gr_tower_t T) { return T->num_trans; }
ulong gr_tower_version(const gr_tower_t T) { return T->version; }
const char * gr_tower_gen_name(const gr_tower_t T, slong d) { return GR_TOWER_GEN(T, d)->name; }
int gr_tower_gen_status(const gr_tower_t T, slong d) { return GR_TOWER_GEN(T, d)->status; }
int gr_tower_gen_kind(const gr_tower_t T, slong d) { return GR_TOWER_GEN(T, d)->kind; }

/* The generator with definition order d as an element of the top field. */
int
gr_tower_gen_get(gr_ptr res, const gr_tower_t T, slong d)
{
    const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
    gr_tower_struct * Tm = (gr_tower_struct *) T;

    if (g->kind == GR_TOWER_ALGEBRAIC)
    {
        gr_ctx_struct * Fk = gr_tower_field_at(T, g->index);
        gr_ptr t;
        int status;
        GR_TMP_INIT(t, Fk);
        status = gr_gen(t, Fk);
        status |= gr_tower_promote(res, t, g->index, T->length, Tm);
        GR_TMP_CLEAR(t, Fk);
        return status;
    }
    else
    {
        gr_vec_t gens;
        int status;
        gr_vec_init(gens, 0, T->base);
        status = gr_gens(gens, T->base);
        if (status == GR_SUCCESS && g->index - 1 < gens->length)
            status = gr_tower_promote(res, gr_vec_entry_ptr(gens, g->index - 1, T->base), 0, T->length, Tm);
        else
            status = GR_UNABLE;
        gr_vec_clear(gens, T->base);
        return status;
    }
}
