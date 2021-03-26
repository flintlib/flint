/*
    Copyright (C) 2008, 2009, William Hart 
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"


void fmpz_comb_temp_clear(fmpz_comb_temp_t CT)
{
    _fmpz_vec_clear(CT->A, CT->Alen);
    _fmpz_vec_clear(CT->T, CT->Tlen);
}


void fmpz_comb_clear(fmpz_comb_t C)
{
    flint_free(C->step);
    flint_free(C->packed_multipliers);
    flint_free(C->crt_lu);
    flint_free(C->mod_lu);
    flint_free(C->crt_offsets);
    flint_free(C->mod_offsets);
    fmpz_multi_CRT_clear(C->crt_P);
    fmpz_multi_mod_clear(C->mod_P);
}
