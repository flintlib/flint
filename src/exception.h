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
   FLINT_ERROR,     /* general error */
   FLINT_OVERFLOW,  /* overflow */
   FLINT_IMPINV,    /* impossible inverse */
   FLINT_DOMERR,    /* domain error */
   FLINT_DIVZERO,   /* divide by zero */
   FLINT_EXPOF,     /* exponent overflow */
   FLINT_INEXACT,   /* inexact error */
   FLINT_TEST_FAIL  /* test fail */
} flint_err_t;

FLINT_NORETURN void flint_throw(flint_err_t exc, const char * msg, ...);

#endif
