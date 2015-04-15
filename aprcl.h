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

    Copyright (C) 2015, Vladimir Glazchev

******************************************************************************/

#ifndef APRCL_H
#define APRCL_H

#include <gmp.h>
#include "flint.h"
#include "fmpz_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct 
{
    fmpz_poly_t poly;
    ulong power;
} unnnnn;

typedef unnnnn unity_root[1];

void unity_init(unity_root element, ulong n);
void unity_clear(unity_root element);
void unity_print(unity_root element);
void unity_nth_root(unity_root element, ulong n);
void unity_roots_add(unity_root res, const unity_root element1, const unity_root element2);
void unity_roots_mul(unity_root res, const unity_root element1, const unity_root element2);
void unity_roots_mul_sub(unity_root res, const unity_root element1, const unity_root element2);

mp_ptr f_table(const ulong q);

#endif


