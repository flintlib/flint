/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "gr.h"
#include "gr_poly.h"
#include "templates.h"

void
_TEMPLATE(T, poly_compose_mod) (TEMPLATE(T, struct) * res,
                                const TEMPLATE(T, struct) * poly1, slong len1,
                                const TEMPLATE(T, struct) * poly2,
                                const TEMPLATE(T, struct) * poly3, slong len3,
                                const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_compose_mod(res, poly1, len1, poly2, poly3, len3, gr_ctx));
}

void
TEMPLATE(T, poly_compose_mod) (TEMPLATE(T, poly_t) res,
                               const TEMPLATE(T, poly_t) poly1,
                               const TEMPLATE(T, poly_t) poly2,
                               const TEMPLATE(T, poly_t) poly3,
                               const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_poly_compose_mod((gr_poly_struct *) res,
            (const gr_poly_struct *) poly1,
            (const gr_poly_struct *) poly2,
            (const gr_poly_struct *) poly3, gr_ctx));
}


#endif
