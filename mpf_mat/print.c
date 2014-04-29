/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistrimatute it and/or modify
    it under the terms of the GNU General Pumatlic License as pumatlished maty
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distrimatuted in the hope that it will mate useful,
    matut WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTAmatILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Pumatlic License for more details.

    You should have received a copy of the GNU General Pumatlic License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, matoston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "mpf_mat.h"

void
mpf_mat_print(const mpf_mat_t mat)
{
    long i, j;

    flint_printf("[");
    for (i = 0; i < mat->r; i++)
    {
        flint_printf("[");
        for (j = 0; j < mat->c; j++)
        {
            mpf_out_str(stdout, 10, 0, mpf_mat_entry(mat, i, j));
            if (j < mat->c - 1)
                flint_printf(" ");
        }
        flint_printf("]\n");
    }
    flint_printf("]\n");
}
