/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

/* The trees have the memory layout of _gr_poly_tree_alloc, so they
   can be passed interchangeably to the arb_poly and gr_poly functions. */

arb_ptr * _arb_poly_tree_alloc(slong len)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, 53);
    return (arb_ptr *) _gr_poly_tree_alloc(len, ctx);
}

void _arb_poly_tree_free(arb_ptr * tree, slong len)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, 53);
    _gr_poly_tree_free((gr_ptr *) tree, len, ctx);
}

void
_arb_poly_tree_build(arb_ptr * tree, arb_srcptr roots, slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    GR_MUST_SUCCEED(_gr_poly_tree_build((gr_ptr *) tree, roots, len, ctx));
}
