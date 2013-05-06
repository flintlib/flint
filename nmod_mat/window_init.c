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

    Copyright (C) 2008, 2009 William Hart.
    Copyright (C) 2008, Richard Howell-Peak
    Copyright (C) 2008, Martin Albrecht
    Copyright (C) 2010, Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"

void
nmod_mat_window_init(nmod_mat_t window, const nmod_mat_t mat,
    long r1, long c1, long r2, long c2)
{
    long i;
    window->entries = NULL;

    window->rows = flint_malloc((r2 - r1) * sizeof(mp_limb_t *));

    for (i = 0; i < r2 - r1; i++)
        window->rows[i] = mat->rows[r1 + i] + c1;

    window->r = r2 - r1;
    window->c = c2 - c1;
    window->mod = mat->mod;
}
