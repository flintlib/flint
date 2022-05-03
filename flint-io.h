/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_IO_H
#define FLINT_IO_H

#include <stdio.h>
#include "flint-config.h"

#ifdef __cplusplus
extern "C" {
#endif

FLINT_DLL int parse_fmt(int * floating, const char * fmt);

FLINT_DLL int flint_printf(const char * str, ...);
FLINT_DLL int flint_vprintf(const char * str, va_list ap);
FLINT_DLL int flint_fprintf(FILE * f, const char * str, ...);
FLINT_DLL int flint_sprintf(char * s, const char * str, ...);

FLINT_DLL int flint_scanf(const char * str, ...);
FLINT_DLL int flint_fscanf(FILE * f, const char * str, ...);
FLINT_DLL int flint_sscanf(const char * s, const char * str, ...);

#ifdef __cplusplus
}
#endif

#endif
