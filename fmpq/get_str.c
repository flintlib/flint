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

   Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpq.h"

char * _fmpq_get_str(char * str, int b, const fmpz_t num, const fmpz_t den)
{
    char *s;

    if (str == NULL)
    {
        str = flint_malloc(fmpz_sizeinbase(num, b) + fmpz_sizeinbase(den, b) + 3);

        if (str == NULL)
        {
            printf("Exception (_fmpq_get_str). Not enough memory.\n");
            abort();
        }
    }

    fmpz_get_str(str, b, num);

    if (!fmpz_is_one(den))
    {
        s = str;
        while (*s != '\0')
            s++;

        *s = '/';
        s++;
        fmpz_get_str(s, b, den);
    }

    return str;
}

char * fmpq_get_str(char * str, int b, const fmpq_t f)
{
    return _fmpq_get_str(str, b, fmpq_numref(f), fmpq_denref(f));
}

