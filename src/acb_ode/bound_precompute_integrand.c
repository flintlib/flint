/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "acb_ode.h"
#include "arb_poly.h"

/* [M19], Algorithm 6.1, step 4. */

void
acb_ode_bound_precompute_integrand(arb_poly_t itg_pol, arb_poly_t itg_num,
                                   const acb_ode_group_t group,
                                   const acb_ode_bound_t bound,
                                   const acb_ode_group_bound_t gbound,
                                   slong n0, slong nlogs, slong prec)
{
    mag_ptr maj_all_nums = _mag_vec_init(bound->all_nums_len);

    /* In view of [M19], Eq. (5.9), one can take here for nlogs the actual
       length wrt log(x) of the solutions of interest truncated after their term
       in x^n, instead of τ(n0) as Algo 6.1, step 3 would prescribe. */
    FLINT_ASSERT(gbound->nums_len == bound->all_nums_len - 1);
    acb_ode_bound_rat_vec(maj_all_nums, gbound->nums, gbound->nums_len,
                          gbound->ind, group, gbound->ind_lbound, gbound->stairs,
                          n0, nlogs, prec);

    arb_poly_fit_length(itg_pol, bound->pol_part_len - 1);
    for (slong i = 0; i < bound->pol_part_len - 1; i++)
        arf_set_mag(arb_midref(itg_pol->coeffs + i), maj_all_nums + i);
    _arb_poly_set_length(itg_pol, bound->pol_part_len - 1);
    _arb_poly_normalise(itg_pol);

    slong num_len = bound->all_nums_len - bound->pol_part_len;
    arb_poly_fit_length(itg_num, num_len);
    _arb_poly_set_length(itg_num, num_len);
    for (slong i = 0; i < num_len; i++)
        arf_set_mag(arb_midref(itg_num->coeffs + i),
                    maj_all_nums + bound->pol_part_len - 1 + i);
    _arb_poly_normalise(itg_num);

    _mag_vec_clear(maj_all_nums, bound->all_nums_len);
}
