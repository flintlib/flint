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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"
#include "arith.h"


#if FLINT64
#define FLINT_HARMONIC_MAX_TINY 46
#else
#define FLINT_HARMONIC_MAX_TINY 24
#endif

const mp_limb_t FLINT_HARMONIC_TINY_P[] = 
{
  0UL, 1UL, 3UL, 11UL, 25UL, 137UL, 49UL, 363UL, 761UL, 7129UL, 7381UL,
  83711UL, 86021UL, 1145993UL, 1171733UL, 1195757UL, 2436559UL, 42142223UL,
  14274301UL, 275295799UL, 55835135UL, 18858053UL, 19093197UL, 444316699UL,
  1347822955UL,
#if FLINT64
  34052522467UL, 34395742267UL, 312536252003UL, 315404588903UL,
  9227046511387UL, 9304682830147UL, 290774257297357UL, 586061125622639UL,
  53676090078349UL, 54062195834749UL, 54437269998109UL, 54801925434709UL,
  2040798836801833UL, 2053580969474233UL, 2066035355155033UL,
  2078178381193813UL, 85691034670497533UL, 12309312989335019UL,
  532145396070491417UL, 5884182435213075787UL, 5914085889685464427UL,
  5943339269060627227UL,
#endif
};

const mp_limb_t FLINT_HARMONIC_TINY_Q[] =
{
  1UL, 1UL, 2UL, 6UL, 12UL, 60UL, 20UL, 140UL, 280UL, 2520UL, 2520UL,
  27720UL, 27720UL, 360360UL, 360360UL, 360360UL, 720720UL, 12252240UL,
  4084080UL, 77597520UL, 15519504UL, 5173168UL, 5173168UL, 118982864UL,
  356948592UL,
#if FLINT64
  8923714800UL, 8923714800UL, 80313433200UL, 80313433200UL, 2329089562800UL,
  2329089562800UL, 72201776446800UL, 144403552893600UL, 13127595717600UL,
  13127595717600UL, 13127595717600UL, 13127595717600UL, 485721041551200UL,
  485721041551200UL, 485721041551200UL, 485721041551200UL,
  19914562703599200UL, 2844937529085600UL, 122332313750680800UL,
  1345655451257488800UL, 1345655451257488800UL, 1345655451257488800UL,
#endif
};

static void
_mpq_harmonic_odd_balanced(fmpz_t num, fmpz_t den, long n)
{
    mpz_t p, q;

    mp_ptr t, v;
    mp_size_t ts, vs;
    long size;

    if (n <= 0)
    {
        fmpz_zero(num);
        fmpz_one(den);
        return;
    }

    /* TODO: we could avoid the copying/allocation overhead when there
       is guaranteed to be sufficient space in res already */

    size = FLINT_BIT_COUNT(n) * (n+2) + 2*FLINT_BITS;
    mpz_init2(p, size);
    mpz_init2(q, size);
    t = p->_mp_d;
    v = q->_mp_d;

    flint_mpn_harmonic_odd_balanced(t, &ts, v, &vs, 1, n+1, n, 1);
    p->_mp_size = ts;
    q->_mp_size = vs;

    fmpz_set_mpz(num, p);
    fmpz_set_mpz(den, q);

    mpz_clear(p);
    mpz_clear(q);

    _fmpq_canonicalise(num, den);
}

void _arith_harmonic_number(fmpz_t num, fmpz_t den, long n)
{
    n = FLINT_MAX(n, 0);

    if (n <= FLINT_HARMONIC_MAX_TINY)
    {
        fmpz_set_ui(num, FLINT_HARMONIC_TINY_P[n]);
        fmpz_set_ui(den, FLINT_HARMONIC_TINY_Q[n]);
    }
    else
    {
        _mpq_harmonic_odd_balanced(num, den, n);
    }
}

void arith_harmonic_number(fmpq_t x, long n)
{
    _arith_harmonic_number(fmpq_numref(x), fmpq_denref(x), n);
}
