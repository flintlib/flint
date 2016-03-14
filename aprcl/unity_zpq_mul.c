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
unity_zpq_mul(unity_zpq f, const unity_zpq g, const unity_zpq h)
{
    ulong i, j, k, p, q;
    fmpz_mod_poly_t temp;

    q = f->q;
    p = f->p;
    fmpz_mod_poly_init(temp, f->n);

    for (i = 0; i < p; i++)
    {
        fmpz_mod_poly_zero(f->polys[i]);
    }

    for (i = 0; i < p; i++)
    {
        for (j = 0; j < p; j++)
        {
            ulong qpow;

            qpow = n_addmod(i, j, p);
            fmpz_mod_poly_mul(temp, g->polys[i], h->polys[j]);

            if (temp->length == 0)
                continue;

            for (k = temp->length - 1; k >= q; k--)
            {
                fmpz_add(temp->coeffs + k - q,
                        temp->coeffs + k - q, temp->coeffs + k);
                fmpz_set_ui(temp->coeffs + k, 0);
                fmpz_mod(temp->coeffs + k - q, temp->coeffs + k - q, f->n);
            }
            _fmpz_mod_poly_normalise(temp);

            fmpz_mod_poly_add(f->polys[qpow], f->polys[qpow], temp);
        }
    }

    fmpz_mod_poly_clear(temp);
}

