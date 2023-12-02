/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>

#include "fmpz_poly_q.h"

/**
 * \ingroup  StringConversions
 *
 * Returns the string representation of the rational function \c op.
 */
char * fmpz_poly_q_get_str(const fmpz_poly_q_t op)
{
    int i, j;
    char * str;
    char * numstr;
    char * denstr;

    if (fmpz_poly_is_one(op->den))
    {
        numstr = fmpz_poly_get_str(op->num);
        i = strlen(numstr) - 1;
        if (numstr[i] == ' ')
        {
            numstr[i] = '\0';
        }
        return numstr;
    }

    numstr = fmpz_poly_get_str(op->num);
    denstr = fmpz_poly_get_str(op->den);

    i = strlen(numstr) - 1;
    if (numstr[i] == ' ')
        numstr[i] = '\0';
    i = strlen(denstr) - 1;
    if (denstr[i] == ' ')
        denstr[i] = '\0';

    str = flint_malloc(strlen(numstr) + strlen(denstr) + 2);
    if (str == NULL)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_q_get_str). Memory allocation failed.\n");
    }

    for (i = 0; i < strlen(numstr); i++)
        str[i] = numstr[i];
    str[i++] = '/';
    for (j = 0; j < strlen(denstr); j++)
        str[i++] = denstr[j];
    str[i] = '\0';

    flint_free(numstr);
    flint_free(denstr);

    return str;
}
