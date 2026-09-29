/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"
#include "gr.h"
#include "gr_poly.h"

int fmpz_mod_poly_is_irreducible_rabin(const fmpz_mod_poly_t f,
                                                      const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    truth_t res;

    /* as special cases, zero and constants are considered irreducible */
    if (f->length <= 2)
        return 1;

    /* the context borrows ctx and must not be cleared */
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    res = gr_poly_is_irreducible_rabin((const gr_poly_struct *) f, gr_ctx);

    if (res == T_UNKNOWN)
        flint_throw(FLINT_ERROR, "fmpz_mod_poly_is_irreducible_rabin: "
                                 "unable to decide\n");

    return res == T_TRUE;
}
