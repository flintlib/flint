/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_poly_factor.h"

/* The distinct degree factorization uses several threads by itself when
   they are available (see gr_poly_factor_distinct_deg). */
void
nmod_poly_factor_distinct_deg_threaded(nmod_poly_factor_t res,
    const nmod_poly_t poly, slong * const * degs)
{
    nmod_poly_factor_distinct_deg(res, poly, degs);
}
