/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_g2_covariant(acb_poly_t r, const fmpz_mpoly_t pol,
    acb_poly_struct* basic, const fmpz_mpoly_ctx_t ctx, slong prec)
{
    acb_poly_eval_fmpz_mpoly(r, pol, basic, ctx, prec);
}
