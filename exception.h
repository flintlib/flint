/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef EXCEPTION_H
#define EXCEPTION_H

typedef enum
{
   FLINT_ERROR,  /* general error */
   FLINT_IMPINV, /* impossible inverse */
   FLINT_DOMERR, /* domain error */
   FLINT_DIVZERO, /* divide by zero */
   FLINT_EXPOF,   /* exponent overflow */
   FLINT_INEXACT /* inexact error */
} flint_err_t;

/* for all except NONE the "extra" is flint-mallocated object */
typedef enum
{
    FLINT_EXC_EXTRA_NONE,
    FLINT_EXC_EXTRA_STR,
    FLINT_EXC_EXTRA_SI,
    FLINT_EXC_EXTRA_UI,
    FLINT_EXC_EXTRA_FMPZ,
    FLINT_EXC_EXTRA_FMPQ
} flint_exc_extra_t;

FLINT_DLL void flint_throw(flint_err_t exc, const char * msg, ...);

FLINT_DLL FLINT_NORETURN void flint_exception(flint_err_t exc, const char * msg,
                                   flint_exc_extra_t extra_type, void * extra);

#endif
