/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NF_H
#define NF_H

#include "nf_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************

    Initialisation

******************************************************************************/

void nf_init(nf_t nf, const fmpq_poly_t pol);

void nf_init_randtest(nf_t nf, flint_rand_t state, slong len,  mp_bitcnt_t bits_in);

void nf_clear(nf_t nf);

void nf_print(const nf_t nf);

#ifdef __cplusplus
}
#endif

#endif
