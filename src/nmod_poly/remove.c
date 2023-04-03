/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"

ulong nmod_poly_remove(nmod_poly_t f, const nmod_poly_t p)
{
    nmod_poly_t q, r;
    ulong i = 0;

    nmod_poly_init_preinv(q, p->mod.n, p->mod.ninv);
    nmod_poly_init_preinv(r, p->mod.n, p->mod.ninv);

    while (1)
    {
        if (f->length < p->length)
            break;
        nmod_poly_divrem(q, r, f, p);
        if (r->length == 0)
            nmod_poly_swap(q, f);
        else
            break;
        i++;
    }

    nmod_poly_clear(q);
    nmod_poly_clear(r);

    return i;
}
