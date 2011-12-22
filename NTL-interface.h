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

    NTL-interface.h: Header file for NTL-interface.cpp

    Copyright (C) 2007 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#ifndef FLINT_NTL_INT_H
#define FLINT_NTL_INT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

NTL_CLIENT

/*
   Returns the number of limbs taken up by an NTL ZZ.
 */
mp_size_t ZZ_limbs(const ZZ& op);

/* 
   Converts an NTL ZZ to an fmpz_t.

   Assumes the fmpz_t has already been allocated to have sufficient space.
 */
void fmpz_set_ZZ(fmpz_t rop, const ZZ& op);

/* 
   Converts an fmpz_t to an NTL ZZ. Allocation is automatically handled.
 */
void fmpz_get_ZZ(ZZ& rop, const fmpz_t op);

/*
   Converts an fmpz_poly_t to an NTL ZZX.
 */
void fmpz_poly_get_ZZX(ZZX& rop, const fmpz_poly_t op);

/*
   Converts an NTL ZZX to an fmpz_poly_t.
*/
void fmpz_poly_set_ZZX(fmpz_poly_t rop, const ZZX& op);

#ifdef __cplusplus
}
#endif

#endif

