/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_discriminant(fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A,
                                     slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    fq_nmod_mpoly_univar_t Ax;

    fq_nmod_mpoly_univar_init(Ax, ctx);

    fq_nmod_mpoly_to_univar(Ax, A, var, ctx);

    success = fq_nmod_mpoly_univar_discriminant(R, Ax, ctx);

    fq_nmod_mpoly_univar_clear(Ax, ctx);

    return success;
}

