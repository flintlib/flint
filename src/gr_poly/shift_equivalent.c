/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_poly.h"

truth_t
gr_poly_shift_equivalent(fmpz_t shift, const gr_poly_t p, const gr_poly_t q,
			 gr_ctx_t ctx)
{
    gr_ptr elt;
    fmpz_t _shift;

    GR_TMP_INIT(elt, ctx);
    fmpz_init(_shift);

    int status = GR_SUCCESS;

    status |= gr_poly_leading_taylor_shift(elt, p, q, ctx);
    status |= gr_get_fmpz(_shift, elt, ctx);

    /* when status == GR_SUCCESS, deg(p) is well-defined */
    if (status == GR_SUCCESS && p->length > 2)
    {
	gr_poly_t pshift;
	gr_poly_init(pshift, ctx);

	status |= gr_poly_taylor_shift(pshift, p, elt, ctx);
	status |= gr_check(gr_poly_equal(pshift, q, ctx));

	gr_poly_clear(pshift, ctx);
    }

    if (status == GR_SUCCESS && shift != NULL)
	fmpz_swap(shift, _shift);

    fmpz_clear(_shift);
    GR_TMP_CLEAR(elt, ctx);

    return gr_in_domain(status);
}
