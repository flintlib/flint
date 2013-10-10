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

#include "mpn_extras.h"
#include "arith.h"

#if FLINT64
#define FLINT_HARMONIC_MAX_TINY 46
#else
#define FLINT_HARMONIC_MAX_TINY 24
#endif

const mp_limb_t FLINT_HARMONIC_TINY_P[] = 
{
  UWORD(0), UWORD(1), UWORD(3), UWORD(11), UWORD(25), UWORD(137), UWORD(49), UWORD(363), UWORD(761), UWORD(7129), UWORD(7381),
  UWORD(83711), UWORD(86021), UWORD(1145993), UWORD(1171733), UWORD(1195757), UWORD(2436559), UWORD(42142223),
  UWORD(14274301), UWORD(275295799), UWORD(55835135), UWORD(18858053), UWORD(19093197), UWORD(444316699),
  UWORD(1347822955),
#if FLINT64
  UWORD(34052522467), UWORD(34395742267), UWORD(312536252003), UWORD(315404588903),
  UWORD(9227046511387), UWORD(9304682830147), UWORD(290774257297357), UWORD(586061125622639),
  UWORD(53676090078349), UWORD(54062195834749), UWORD(54437269998109), UWORD(54801925434709),
  UWORD(2040798836801833), UWORD(2053580969474233), UWORD(2066035355155033),
  UWORD(2078178381193813), UWORD(85691034670497533), UWORD(12309312989335019),
  UWORD(532145396070491417), UWORD(5884182435213075787), UWORD(5914085889685464427),
  UWORD(5943339269060627227),
#endif
};

const mp_limb_t FLINT_HARMONIC_TINY_Q[] =
{
  UWORD(1), UWORD(1), UWORD(2), UWORD(6), UWORD(12), UWORD(60), UWORD(20), UWORD(140), UWORD(280), UWORD(2520), UWORD(2520),
  UWORD(27720), UWORD(27720), UWORD(360360), UWORD(360360), UWORD(360360), UWORD(720720), UWORD(12252240),
  UWORD(4084080), UWORD(77597520), UWORD(15519504), UWORD(5173168), UWORD(5173168), UWORD(118982864),
  UWORD(356948592),
#if FLINT64
  UWORD(8923714800), UWORD(8923714800), UWORD(80313433200), UWORD(80313433200), UWORD(2329089562800),
  UWORD(2329089562800), UWORD(72201776446800), UWORD(144403552893600), UWORD(13127595717600),
  UWORD(13127595717600), UWORD(13127595717600), UWORD(13127595717600), UWORD(485721041551200),
  UWORD(485721041551200), UWORD(485721041551200), UWORD(485721041551200),
  UWORD(19914562703599200), UWORD(2844937529085600), UWORD(122332313750680800),
  UWORD(1345655451257488800), UWORD(1345655451257488800), UWORD(1345655451257488800),
#endif
};

static void
_mpq_harmonic_odd_balanced(fmpz_t num, fmpz_t den, slong n)
{
    mpz_t p, q;

    mp_ptr t, v;
    mp_size_t ts, vs;
    slong size;

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

void _arith_harmonic_number(fmpz_t num, fmpz_t den, slong n)
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

void arith_harmonic_number(fmpq_t x, slong n)
{
    _arith_harmonic_number(fmpq_numref(x), fmpq_denref(x), n);
}
