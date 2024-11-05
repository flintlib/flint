/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_DIRICHLET_IMPL_H
#define ACB_DIRICHLET_IMPL_H

#include "acb_dirichlet.h"

#ifdef __cplusplus
extern "C" {
#endif

void get_smk_points(fmpz * res, slong A, slong B);

slong platt_get_smk_index(slong B, const fmpz_t j, slong prec);

void
_platt_smk(acb_ptr table, acb_ptr startvec, acb_ptr stopvec,
        const fmpz * smk_points, const arb_t t0, slong A, slong B,
        const fmpz_t jstart, const fmpz_t jstop, slong mstart, slong mstop,
        slong K, slong prec);

void
_acb_dirichlet_platt_multieval(arb_ptr out, acb_srcptr S_table,
        const arb_t t0, slong A, slong B, const arb_t h, const fmpz_t J,
        slong K, slong sigma, slong prec);

void acb_dirichlet_stieltjes_integral(acb_t res, const fmpz_t n, const acb_t a, slong prec);
void acb_dirichlet_stieltjes_em(acb_t res, const fmpz_t n, const acb_t a, slong prec);

/* FIXME: This function is unused, untested and undocumented. */
void
acb_dirichlet_vec_mellin_arb(acb_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, const arb_t t, slong n, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* ACB_DIRICHLET_IMPL_H */
