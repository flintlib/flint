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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

void
elem_set_si(elem_t elem, long v, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_set_si(&elem->z, v);
            break;

        case TYPE_LIMB:
            elem->n = v;
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_FMPZ:
                        fmpz_set_si(&elem->z, v);
                        fmpz_mod(&elem->z, &elem->z, &RING_MODULUS(ring)->z);
                        break;

                    case TYPE_LIMB:
                        elem->n = n_mod2_preinv(FLINT_ABS(v),
                            ring->nmod.n, ring->nmod.ninv);
                        if (v < 0)
                            elem->n = n_negmod(elem->n, ring->nmod.n);
                        break;

                    default:
                        NOT_IMPLEMENTED("set_si (mod)", ring);
                }
            }
            break;

        default:
            NOT_IMPLEMENTED("set_si", ring);
    }
}

void
gen_set_si(gen_t x, long v)
{
    elem_set_si(&x->elem, v, x->ring);
}
