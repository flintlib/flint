/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


int fq_zech_mpoly_univar_content_mpoly(
    fq_zech_mpoly_t g,
    const fq_zech_mpoly_univar_t A,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i;

    fq_zech_mpoly_zero(g, ctx);
    for (i = 0; i < A->length; i++)
    {
		if (!fq_zech_mpoly_gcd(g, g, A->coeffs + i, ctx))
			return 0;
    }

    return 1;
}

