/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"


static void
_gr_uninitialized_ctx_clear(gr_ctx_t ctx)
{
    /* For debugging: flint_printf("CLEARING UNINTIIALIZED\n"); */
}

static int _gr_uninitialized_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    return gr_stream_write(out, "Uninitialized context (init failed)");
}

static int
_gr_uninitialized_init(gr_ptr res, gr_ctx_t ctx)
{
    flint_throw(FLINT_ERROR, "Cannot initialize element of uninitialized context");
    return GR_UNABLE;
}

int _uninitialized_methods_initialized = 0;

gr_static_method_table _uninitialized_methods;

gr_method_tab_input _uninitialized_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,     (gr_funcptr) _gr_uninitialized_ctx_clear},
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_uninitialized_ctx_write},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_uninitialized_init},
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_uninitialized(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_UNINITIALIZED;
    ctx->sizeof_elem = 1;
    ctx->size_limit = WORD_MAX;

    ctx->methods = _uninitialized_methods;

    if (!_uninitialized_methods_initialized)
    {
        gr_method_tab_init(_uninitialized_methods, _uninitialized_methods_input);
        _uninitialized_methods_initialized = 1;
    }
}

