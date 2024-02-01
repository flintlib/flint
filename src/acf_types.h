/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACF_TYPES_H
#define ACF_TYPES_H

#include "arf_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    arf_struct real;
    arf_struct imag;
}
acf_struct;

typedef acf_struct acf_t[1];
typedef acf_struct * acf_ptr;
typedef const acf_struct * acf_srcptr;

#ifdef __cplusplus
}
#endif

#endif /* ACF_TYPES_H */
