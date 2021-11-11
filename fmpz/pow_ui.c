/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "longlong.h"

#if FLINT_BITS == 64
#define EXPBOUND 63
#else
#define EXPBOUND 39
#endif

static void
mul0(fmpz_t res, slong a)
{ /* a can be whatever, even an mpz */
    fmpz_one(res);
}

static void
mul1(fmpz_t res, slong a)
{ /* a can be whatever, even an mpz */
    fmpz_set(res, &a);
}

static void
mul2(fmpz_t res, slong a)
{ /* Largest slong is a */
    ulong rh, rl;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul3(fmpz_t res, slong a)
{ /* Largest slong is a^2 */
    ulong rh, rl, b;
    b = a * a;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul4(fmpz_t res, slong a)
{ /* Largest slong is a^2 */
    slong rh, rl;
    a = a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul5(fmpz_t res, slong a)
{ /* Largest slong is a^3 */
    ulong rh, rl, b;
    b = a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul6(fmpz_t res, slong a)
{ /* Largest slong is a^3 */
    slong rh, rl;
    a = a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul7(fmpz_t res, slong a)
{ /* Largest slong is a^4 */
    ulong rh, rl, b;
    b = a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul8(fmpz_t res, slong a)
{ /* Largest slong is a^4 */
    slong rh, rl;
    a = a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul9(fmpz_t res, slong a)
{ /* Largest slong is a^5 */
    ulong rh, rl, b;
    b = a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul10(fmpz_t res, slong a)
{ /* Largest slong is a^5 */
    slong rh, rl;
    a = a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul11(fmpz_t res, slong a)
{ /* Largest slong is a^6 */
    ulong rh, rl, b;
    b = a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul12(fmpz_t res, slong a)
{ /* Largest slong is a^6 */
    slong rh, rl;
    a = a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul13(fmpz_t res, slong a)
{ /* Largest slong is a^7 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul14(fmpz_t res, slong a)
{ /* Largest slong is a^7 */
    slong rh, rl;
    a = a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul15(fmpz_t res, slong a)
{ /* Largest slong is a^8 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul16(fmpz_t res, slong a)
{ /* Largest slong is a^8 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul17(fmpz_t res, slong a)
{ /* Largest slong is a^9 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul18(fmpz_t res, slong a)
{ /* Largest slong is a^9 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul19(fmpz_t res, slong a)
{ /* Largest slong is a^10 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul20(fmpz_t res, slong a)
{ /* Largest slong is a^10 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul21(fmpz_t res, slong a)
{ /* Largest slong is a^11 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul22(fmpz_t res, slong a)
{ /* Largest slong is a^11 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul23(fmpz_t res, slong a)
{ /* Largest slong is a^12 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul24(fmpz_t res, slong a)
{ /* Largest slong is a^12 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul25(fmpz_t res, slong a)
{ /* Largest slong is a^13 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul26(fmpz_t res, slong a)
{ /* Largest slong is a^13 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul27(fmpz_t res, slong a)
{ /* Largest slong is a^14 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul28(fmpz_t res, slong a)
{ /* Largest slong is a^14 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul29(fmpz_t res, slong a)
{ /* Largest slong is a^15 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul30(fmpz_t res, slong a)
{ /* Largest slong is a^15 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul31(fmpz_t res, slong a)
{ /* Largest slong is a^16 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul32(fmpz_t res, slong a)
{ /* Largest slong is a^16 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul33(fmpz_t res, slong a)
{ /* Largest slong is a^17 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul34(fmpz_t res, slong a)
{ /* Largest slong is a^17 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul35(fmpz_t res, slong a)
{ /* Largest slong is a^18 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul36(fmpz_t res, slong a)
{ /* Largest slong is a^18 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul37(fmpz_t res, slong a)
{ /* Largest slong is a^19 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul38(fmpz_t res, slong a)
{ /* Largest slong is a^19 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul39(fmpz_t res, slong a)
{ /* Largest slong is a^20 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

#if FLINT_BITS == 64

static void
mul40(fmpz_t res, slong a)
{ /* Largest slong is a^20 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul41(fmpz_t res, slong a)
{ /* Largest slong is a^21 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul42(fmpz_t res, slong a)
{ /* Largest slong is a^21 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul43(fmpz_t res, slong a)
{ /* Largest slong is a^22 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul44(fmpz_t res, slong a)
{ /* Largest slong is a^22 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul45(fmpz_t res, slong a)
{ /* Largest slong is a^23 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul46(fmpz_t res, slong a)
{ /* Largest slong is a^23 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul47(fmpz_t res, slong a)
{ /* Largest slong is a^24 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul48(fmpz_t res, slong a)
{ /* Largest slong is a^24 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul49(fmpz_t res, slong a)
{ /* Largest slong is a^25 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul50(fmpz_t res, slong a)
{ /* Largest slong is a^25 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul51(fmpz_t res, slong a)
{ /* Largest slong is a^26 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul52(fmpz_t res, slong a)
{ /* Largest slong is a^26 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul53(fmpz_t res, slong a)
{ /* Largest slong is a^27 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul54(fmpz_t res, slong a)
{ /* Largest slong is a^27 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul55(fmpz_t res, slong a)
{ /* Largest slong is a^28 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul56(fmpz_t res, slong a)
{ /* Largest slong is a^28 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul57(fmpz_t res, slong a)
{ /* Largest slong is a^29 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul58(fmpz_t res, slong a)
{ /* Largest slong is a^29 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul59(fmpz_t res, slong a)
{ /* Largest slong is a^30 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul60(fmpz_t res, slong a)
{ /* Largest slong is a^30 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul61(fmpz_t res, slong a)
{ /* Largest slong is a^31 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul62(fmpz_t res, slong a)
{ /* Largest slong is a^31 */
    slong rh, rl;
    a = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    smul_ppmm(rh, rl, a, a);
    fmpz_set_signed_uiui(res, rh, rl);
}

static void
mul63(fmpz_t res, slong a)
{ /* Largest slong is a^32 */
    ulong rh, rl, b;
    b = a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a * a;
    a *= b;
    smul_ppmm(rh, rl, a, b);
    fmpz_set_signed_uiui(res, rh, rl);
}

#endif

typedef void (* mulfunc_t)(fmpz_t, slong);

static const mulfunc_t mulfunc[] = {
    mul0, mul1, mul2, mul3, mul4, mul5, mul6, mul7, mul8, mul9,
    mul10,mul11,mul12,mul13,mul14,mul15,mul16,mul17,mul18,mul19,
    mul20,mul21,mul22,mul23,mul24,mul25,mul26,mul27,mul28,mul29,
    mul30,mul31,mul32,mul33,mul34,mul35,mul36,mul37,mul38,mul39,
#if FLINT_BITS == 64
    mul40,mul41,mul42,mul43,mul44,mul45,mul46,mul47,mul48,mul49,
    mul50,mul51,mul52,mul53,mul54,mul55,mul56,mul57,mul58,mul59,
    mul60,mul61,mul62,mul63
#endif
};

/* We do not have a limit for mul0 and mul1. */
static const ulong mulfunc_bound[] = {
#if FLINT_BITS == 64
    0xffffffffffffffff,                         /* mul0            bound */
    0xffffffffffffffff,4611686018427387903,     /* mul1  and mul2  bound */
    3037000499,3037000499,                      /* mul3  and mul4  bound */
    2097151,2097151,                            /* mul5  and mul6  bound */
    55108,55108,                                /* mul7  and mul8  bound */
    6208,6208,                                  /* mul9  and mul10 bound */
    1448,1448,                                  /* mul11 and mul12 bound */
    511,511,                                    /* mul13 and mul14 bound */
    234,234,                                    /* mul15 and mul16 bound */
    127,127,                                    /* mul17 and mul18 bound */
    78,78,                                      /* mul19 and mul20 bound */
    52,52,                                      /* mul21 and mul22 bound */
    38,38,                                      /* mul23 and mul24 bound */
    28,28,                                      /* mul25 and mul26 bound */
    22,22,                                      /* mul27 and mul28 bound */
    18,18,                                      /* mul29 and mul30 bound */
    15,15,                                      /* mul31 and mul32 bound */
    13,13,                                      /* mul33 and mul34 bound */
    11,11,                                      /* mul35 and mul36 bound */
    9,9,                                        /* mul37 and mul38 bound */
    8,8,                                        /* mul39 and mul40 bound */
    7,7,                                        /* mul41 and mul42 bound */
    7,7,                                        /* mul43 and mul44 bound */
    6,6,                                        /* mul45 and mul46 bound */
    6,6,                                        /* mul45 and mul46 bound */
    6,6,                                        /* mul47 and mul48 bound */
    5,5,                                        /* mul49 and mul50 bound */
    5,5,                                        /* mul49 and mul50 bound */
    5,5,                                        /* mul51 and mul52 bound */
    5,5,                                        /* mul53 and mul54 bound */
    4,4,                                        /* mul55 and mul56 bound */
    4,4,                                        /* mul57 and mul58 bound */
    4,4,                                        /* mul59 and mul60 bound */
    4,4,                                        /* mul61 and mul62 bound */
    3                                           /* mul63           bound */
#else
    0xffffffff,                                 /* mul0            bound */
    0xffffffff,1073741823,                      /* mul1  and mul2  bound */
    46340,46340,                                /* mul3  and mul4  bound */
    1290,1290,                                  /* mul5  and mul6  bound */
    215,215,                                    /* mul7  and mul8  bound */
    73,73,                                      /* mul9  and mul10 bound */
    35,35,                                      /* mul11 and mul12 bound */
    21,21,                                      /* mul13 and mul14 bound */
    14,14,                                      /* mul15 and mul16 bound */
    10,10,                                      /* mul17 and mul18 bound */
    8,8,                                        /* mul19 and mul20 bound */
    7,7,                                        /* mul21 and mul22 bound */
    5,5,                                        /* mul23 and mul24 bound */
    5,5,                                        /* mul25 and mul26 bound */
    4,4,                                        /* mul27 and mul28 bound */
    4,4,                                        /* mul29 and mul30 bound */
    3,3,                                        /* mul31 and mul32 bound */
    3,3,                                        /* mul33 and mul34 bound */
    3,3,                                        /* mul35 and mul36 bound */
    3,3,                                        /* mul37 and mul38 bound */
    3                                           /* mul39           bound */
#endif
};


void
fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz s1 = *g;
    __mpz_struct * mptr;

    /* We do not need to check COEFF_IS_MPZ here, since values in mulfunc_bound
     * can be represented as slongs (apart from when exp equals 0 or 1). */
    if (exp <= EXPBOUND
            && s1 <= mulfunc_bound[exp]
            && -s1 <= mulfunc_bound[exp])
    {
        mulfunc[exp](f, s1);
        return;
    }

    if (s1 <= 1 && s1 >= -1)
    {
        fmpz_set_si(f, (s1 == 1)
                     + (s1 == -1 && (exp & 1) == 0)
                     - (s1 == -1 && (exp & 1) == 1));
        return;
    }

    mptr = _fmpz_promote(f);
    if (!COEFF_IS_MPZ(s1))
    {
        /* s1 cannot be zero at this stage */
        __mpz_struct temp[1] = {{ 1, 1, NULL }};
        temp->_mp_d = (mp_ptr) &s1;

        if (s1 < 0)
        {
            temp->_mp_size = -1;
            s1 = -s1;
        }

        flint_mpz_pow_ui(mptr, temp, exp);
    }
    else
    {
        flint_mpz_pow_ui(mptr, COEFF_TO_PTR(s1), exp);
    }
}
