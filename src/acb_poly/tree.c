/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "gr_poly.h"

/* The trees have the memory layout of _gr_poly_tree_alloc, so they
   can be passed interchangeably to the acb_poly and gr_poly functions. */

acb_ptr * _acb_poly_tree_alloc(slong len)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, 53);
    return (acb_ptr *) _gr_poly_tree_alloc(len, ctx);
}

void _acb_poly_tree_free(acb_ptr * tree, slong len)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, 53);
    _gr_poly_tree_free((gr_ptr *) tree, len, ctx);
}

void
_acb_poly_tree_build(acb_ptr * tree, acb_srcptr roots, slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    GR_MUST_SUCCEED(_gr_poly_tree_build((gr_ptr *) tree, roots, len, ctx));
}
