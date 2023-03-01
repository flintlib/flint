/*
    Copyright (C) 2007 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_NTL_INT_H
#define FLINT_NTL_INT_H

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pEX.h>
#include <NTL/vec_ZZ.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "fq.h"
#include "fq_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
   Converts an NTL ZZ to an fmpz_t.

   Assumes the fmpz_t has already been allocated to have sufficient space.
*/
FLINT_DLL void fmpz_set_ZZ(fmpz_t rop, const NTL_NNS ZZ& op);

/*
   Converts an fmpz_t to an NTL ZZ. Allocation is automatically handled.
 */
FLINT_DLL void fmpz_get_ZZ(NTL_NNS ZZ& rop, const fmpz_t op);


/*
   Converts an NTL ZZ_p to an fmpz_t.

   Assumes the fmpz_t has already been allocated to have sufficient space.
*/
FLINT_DLL void fmpz_set_ZZ_p(fmpz_t rop, const NTL_NNS ZZ_p& op);

/*
   Converts an fmpz_t to an NTL ZZ_p. Allocation is automatically handled.
 */
FLINT_DLL void fmpz_get_ZZ_p(NTL_NNS ZZ_p& rop, const fmpz_t op);

/*
   Converts an NTL zz_p to an fmpz_t.
*/
FLINT_DLL void fmpz_set_zz_p(fmpz_t rop, const NTL_NNS zz_p& op);

/*
   Converts an fmpz_t to an NTL zz_p.
 */
FLINT_DLL void fmpz_get_zz_p(NTL_NNS zz_p& rop, const fmpz_t op);

/*
  Converts an fmpz_poly_t to an NTL ZZX.
*/
FLINT_DLL     void fmpz_poly_get_ZZX(NTL_NNS ZZX& rop, const fmpz_poly_t op);

/*
  Converts an NTL ZZX to an fmpz_poly_t.
*/
FLINT_DLL     void fmpz_poly_set_ZZX(fmpz_poly_t rop, const NTL_NNS ZZX& op);

/*
  Converts an fmpz_mod_poly_t to an NTL ZZ_pX.
*/
FLINT_DLL void fmpz_mod_poly_get_ZZ_pX(NTL_NNS ZZ_pX& rop, const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx);

/*
  Converts an NTL ZZ_pX to an fmpz_poly_t.
*/
FLINT_DLL void fmpz_mod_poly_set_ZZ_pX(fmpz_mod_poly_t rop, const NTL_NNS ZZ_pX& op, const fmpz_mod_ctx_t ctx);

/*
  Converts an fq_t to an NTL ZZ_pE.
*/
FLINT_DLL void fq_get_ZZ_pE(NTL_NNS ZZ_pE& rop, const fq_t op, const fq_ctx_t ctx);

/*
  Converts an NTL ZZ_pE to an fq_t.
*/
FLINT_DLL void fq_set_ZZ_pE(fq_t rop, const NTL_NNS ZZ_pE& op, const fq_ctx_t ctx);


/*
  Converts an fq_poly_t to an NTL ZZ_pEX.
*/
FLINT_DLL void fq_poly_get_ZZ_pEX(NTL_NNS ZZ_pEX& rop, const fq_poly_t op, const fq_ctx_t ctx);

/*
  Converts an NTL ZZ_pEX to an fq_poly_t.
*/
FLINT_DLL void fq_poly_set_ZZ_pEX(fq_poly_t rop, const NTL_NNS ZZ_pEX& op, const fq_ctx_t ctx);

/*
  Converts an fmpz_mod_poly_t to an NTL zz_pX.
*/
FLINT_DLL void fmpz_mod_poly_get_zz_pX(NTL_NNS zz_pX& rop, const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx);

/*
  Converts an NTL zz_pX to an fmpz_poly_t.
*/
FLINT_DLL void fmpz_mod_poly_set_zz_pX(fmpz_mod_poly_t rop, const NTL_NNS zz_pX& op, const fmpz_mod_ctx_t ctx);

/*
  Converts an fq_t to an NTL zz_pE.
*/
FLINT_DLL void fq_get_zz_pE(NTL_NNS zz_pE& rop, const fq_t op, const fq_ctx_t ctx);

/*
  Converts an NTL zz_pE to an fq_t.
*/
FLINT_DLL void fq_set_zz_pE(fq_t rop, const NTL_NNS zz_pE& op, const fq_ctx_t ctx);


/*
  Converts an fq_poly_t to an NTL zz_pEX.
*/
FLINT_DLL void fq_poly_get_zz_pEX(NTL_NNS zz_pEX& rop, const fq_poly_t op, const fq_ctx_t ctx);

/*
  Converts an NTL zz_pEX to an fq_poly_t.
*/
FLINT_DLL void fq_poly_set_zz_pEX(fq_poly_t rop, const NTL_NNS zz_pEX& op, const fq_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

