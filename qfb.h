/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart

******************************************************************************/

#ifndef QFB_H
#define QFB_H

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct qfb
{
    fmpz_t a;
    fmpz_t b;
    fmpz_t c;
} qfb;

typedef qfb qfb_t[1];

long qfb_reduced_forms(qfb ** forms, long d);

long qfb_reduced_forms_large(qfb ** forms, long d);

static __inline__
void qfb_array_clear(qfb ** forms, long num)
{
   long k;

   for (k = 0; k < num; k++)
   {
      fmpz_clear((*forms)[k].a);
      fmpz_clear((*forms)[k].b);
      fmpz_clear((*forms)[k].c);
   }
   flint_free(*forms);
}

#ifdef __cplusplus
}
#endif

#endif
