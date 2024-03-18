/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fq_zech.h"
#include "fq_zech_embed.h"

#define T fq_zech
#define B nmod

void TEMPLATE(T, modulus_pow_series_inv)(TEMPLATE(B, poly_t) res,
                                         const TEMPLATE(T, ctx_t) ctx,
                                         slong trunc)
{
    TEMPLATE(B, poly_reverse)(res,
                              TEMPLATE(T, ctx_modulus)(ctx),
                              TEMPLATE(T, ctx_degree)(ctx) + 1);
    TEMPLATE(B, poly_inv_series)(res, res, trunc);
}

#undef B
#undef T
