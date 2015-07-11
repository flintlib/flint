/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

#include "aprcl.h"

void
unity_zp_mont_init(unity_zp_mont f, ulong p, ulong exp, const fmpz_t n)
{
    f->p = p;
    f->exp = exp;
    fmpz_init_set(f->n, n);
    fmpz_init_set_ui(f->r, 0);
    f->r = fmpz_bits(n);
    fmpz_init(f->ninv);
    fmpz_init(f->rinv);
    fmpz_mod_poly_init(f->poly, n);
}

void
unity_zp_mont_clear(unity_zp_mont f)
{
    fmpz_clear(f->n);
    fmpz_clear(f->ninv);
    fmpz_clear(f->rinv);
    fmpz_mod_poly_clear(f->poly);
}

void
unity_zp_mont_reduction(unity_zp_mont f)
{
    slong i;
    fmpz_t m, t;
    fmpz_init(m);
    fmpz_init(t);

    for (i = 0; i < f->poly->length; i++)
    {
        fmpz_fdiv_r_2exp(m, f->poly->coeffs + i, f->r);
        fmpz_mul(t, m, f->ninv);
        fmpz_fdiv_r_2exp(m, t, f->r);
        fmpz_addmul(f->poly->coeffs + i, f->n, m);
        fmpz_fdiv_q_2exp(f->poly->coeffs + i, f->poly->coeffs + i, f->r);
        if (fmpz_cmp(f->n, f->poly->coeffs + i) <= 0)
            fmpz_sub(f->poly->coeffs + i, f->poly->coeffs + i, f->n);
    }

    fmpz_clear(m);
    fmpz_clear(t);
}

