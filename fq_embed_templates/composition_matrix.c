/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#ifdef B

#include "templates.h"

/* The naming of these macros is an abuse */
#define __nmod_poly_get_coeff(p, i) ((p)->coeffs[(i)])
#define __fmpz_mod_poly_get_coeff(p, i) ((p)->coeffs + (i))

void TEMPLATE(T, embed_composition_matrix_sub)(TEMPLATE(B, mat_t) matrix,
                                         const TEMPLATE(T, t) gen,
                                         const TEMPLATE(T, ctx_t) ctx,
                                         slong trunc)
{
    slong i, j, len = TEMPLATE(T, ctx_degree)(ctx);
    TEMPLATE(T, t) tmp;

    TEMPLATE(T, init)(tmp, ctx);
    TEMPLATE(T, one)(tmp, ctx);

    TEMPLATE(B, mat_zero)(matrix);
    for (j = 0; j < trunc; j++)
    {
        for (i = 0; i < tmp->length; i++)
            TEMPLATE(B, mat_set_entry)(matrix, i, j, 
                                       __TEMPLATE(B, poly_get_coeff)(tmp, i));
        if (j < len - 1)
            TEMPLATE(T, mul)(tmp, tmp, gen, ctx);
    }

    TEMPLATE(T, clear)(tmp, ctx);
}

#endif
#endif
