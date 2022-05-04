/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mini.h"

void
_fmpz_vec_demote(fmpz * vec, slong len)
{
    vec += len;
    len = -len;
    for (; len != 0; len++)
       _fmpz_demote(vec + len); 
}
