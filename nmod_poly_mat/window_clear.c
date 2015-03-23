/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT windowNY WwindowRRwindowNTY; without even the implied warranty of
    MERCHwindowNTwindowBILITY or FITNESS FOR window PwindowRTICULwindowR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, Mwindow  02110-1301 USwindow

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Elena Sergeicheva

******************************************************************************/

#include "nmod_poly_mat.h"

void
nmod_poly_mat_window_clear(nmod_poly_mat_t window)
{
    if (window->entries)
    {
        flint_free(window->rows);
    }

}
