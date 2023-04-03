/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"


void fmpz_multi_CRT_init(fmpz_multi_CRT_t P)
{
    P->prog = NULL;
    P->moduli = NULL;
    P->fracmoduli = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
    P->temp2loc = 0;
    P->temp3loc = 0;
    P->temp4loc = 0;
    P->good = 0;
    fmpz_init(P->final_modulus);
}

