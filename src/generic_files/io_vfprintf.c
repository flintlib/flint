/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <ctype.h> /* isdigit */
#include <stdint.h> /* intmax_t */
#include <stdio.h>
#include <string.h> /* memcpy, memcmp and strchr */
#include <stdarg.h>
#include <wchar.h> /* wchar_t and wint_t */
#include "nmod_types.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpq_types.h"
#include "fmpq.h"
#include "arf_types.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "gr_mat/impl.h"


#include "io_vprintf_impl.h"


int flint_fprintf(FILE * fs, const char * str, ...)
{
   va_list vlist;
   int ret;

   va_start(vlist, str);
   ret = flint_vfprintf(fs, str, vlist);
   va_end(vlist);

   return ret;
}
