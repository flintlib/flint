/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz_vec.h"

#define _isdigit(x) ((x) >= '0' && (x) <= '9')
#define NUMBER_LIMIT 512

int
_fmpz_vec_set_str(fmpz * res, slong * len, const char * str)
{
    ulong ix, kt;
    char * t1 = (char *) str;
    char * t2;
    char tmp[NUMBER_LIMIT];

    /* If the first character is not a digit, it is has failed. */
    if (!_isdigit(*t1))
        return -1;

    /* Read length */
    *len = strtol(t1, &t1, 10);

    /* In case the vector has zero length */
    if (*len <= 0)
        return -1;

    /* Has to have two spaces after the length. (The other one is checked later.) */
    if (*t1 != ' ')
        return -1;

    t1 += 2; /* t1 should point to the first digit in the number */
    t2 = t1; /* t2 will point to the first character after the number */
    ix = *len;
    
    do
    {
        /* If the character behind us was not a space, it has failed. */
        if (t1[-1] != ' ')
            return -1;

        /* Increase width until a non-digit is found. */
        kt = NUMBER_LIMIT - 1;
        do
        {
            t2++;
            kt--;
            if (kt == 0)
                return -1;
        } while (_isdigit(*t2));

        kt = t2 - t1;

        /* Copy the digits into tmp */
        memcpy(tmp, t1, sizeof(char) * kt);

        /* Make it null-terminated string and read the string into res. */ 
        tmp[kt] = '\0';
        if (fmpz_set_str(res, tmp, 10))
            return -1;
        
        t1 = ++t2;
        ix--, res++;
    } while (ix > 0);

    if (t2[-1] != '\0')
        return -1;

    return 0;
}

#undef _isdigit
