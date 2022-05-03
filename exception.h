/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef EXCEPTION_H
#define EXCEPTION_H

#include "flint-config.h"

typedef enum
{
    FLINT_ERROR,    /* general error */
    FLINT_ALLOC,    /* allocation error */
    FLINT_MEMMGR,   /* memory manager error */
    FLINT_IMPINV,   /* impossible inverse */
    FLINT_DOMERR,   /* domain error */
    FLINT_DIVZERO,  /* divide by zero */
    FLINT_EXPOF,    /* exponent overflow */
    FLINT_INEXACT   /* inexact error */
} flint_err_t;

FLINT_DLL void flint_throw(flint_err_t exc, const char * msg, ...) FLINT_NORETURN;

#endif
