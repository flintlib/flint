/*
    Copyright (C) 2023 Mathieu Gouttenoire
    
    This file is part of FLINT.
    
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

int
fmpzi_is_probabprime(const fmpzi_t n)
{
    fmpz_t t;
    const fmpz *a, *b;
    
    int res = 0;
    
    fmpz_init(t);
    
    a = fmpzi_realref(n);
    b = fmpzi_imagref(n);
    
    if (fmpz_is_zero(b)) {
        if (fmpz_tdiv_ui(a, 4) == 3)
            fmpz_abs(t, a);
    } else if (fmpz_is_zero(a)) {
        if (fmpz_tdiv_ui(b, 4) == 3)
            fmpz_abs(t, b);
    } else {
        fmpzi_norm(t, n);
    }
    
    res = fmpz_is_probabprime(t);
    
    fmpz_clear(t);
    
    return res;
}
