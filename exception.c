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

    Copyright (C) 2016 William Hart

******************************************************************************/

#include <stdlib.h>
#include <stdarg.h>
#include "flint.h"

void flint_throw(flint_err_t exc, const char * msg, ...)
{
    va_list ap;

    va_start(ap, msg);

    flint_printf("Flint exception (");

    switch (exc)
    {
        case FLINT_ERROR:
            flint_printf("General error");
            break;
        case FLINT_IMPINV:
            flint_printf("Impossible inverse");
            break;
        case FLINT_DOMERR:
            flint_printf("Domain error");
            break;
        case FLINT_DIVZERO:
            flint_printf("Divide by zero");
            break;
        default:
            flint_printf("Unknown exception");
     }

     printf("):\n    ");

     flint_vprintf(msg, ap);
     va_end(ap);

     abort();
}