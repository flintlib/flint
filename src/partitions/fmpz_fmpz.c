/*
    Copyright (C) 2013-2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <limits.h>
#include "mpn_extras.h"
#include "fmpz.h"
#include "mp_real.h"
#include "partitions.h"
#include "partitions/impl.h"

/* This nice round number precisely fits on 32 bits */
#define NUMBER_OF_SMALL_PARTITIONS 128

const unsigned int
partitions_lookup[NUMBER_OF_SMALL_PARTITIONS] =
{
    UWORD(1),UWORD(1),UWORD(2),UWORD(3),UWORD(5),UWORD(7),UWORD(11),UWORD(15),UWORD(22),UWORD(30),UWORD(42),UWORD(56),UWORD(77),UWORD(101),UWORD(135),
    UWORD(176),UWORD(231),UWORD(297),UWORD(385),UWORD(490),UWORD(627),UWORD(792),UWORD(1002),UWORD(1255),UWORD(1575),UWORD(1958),
    UWORD(2436),UWORD(3010),UWORD(3718),UWORD(4565),UWORD(5604),UWORD(6842),UWORD(8349),UWORD(10143),UWORD(12310),UWORD(14883),
    UWORD(17977),UWORD(21637),UWORD(26015),UWORD(31185),UWORD(37338),UWORD(44583),UWORD(53174),UWORD(63261),UWORD(75175),
    UWORD(89134),UWORD(105558),UWORD(124754),UWORD(147273),UWORD(173525),UWORD(204226),UWORD(239943),UWORD(281589),
    UWORD(329931),UWORD(386155),UWORD(451276),UWORD(526823),UWORD(614154),UWORD(715220),UWORD(831820),UWORD(966467),
    UWORD(1121505),UWORD(1300156),UWORD(1505499),UWORD(1741630),UWORD(2012558),UWORD(2323520),UWORD(2679689),
    UWORD(3087735),UWORD(3554345),UWORD(4087968),UWORD(4697205),UWORD(5392783),UWORD(6185689),UWORD(7089500),
    UWORD(8118264),UWORD(9289091),UWORD(10619863),UWORD(12132164),UWORD(13848650),UWORD(15796476),UWORD(18004327),
    UWORD(20506255),UWORD(23338469),UWORD(26543660),UWORD(30167357),UWORD(34262962),UWORD(38887673),
    UWORD(44108109),UWORD(49995925),UWORD(56634173),UWORD(64112359),UWORD(72533807),UWORD(82010177),
    UWORD(92669720),UWORD(104651419),UWORD(118114304),UWORD(133230930),UWORD(150198136),UWORD(169229875),
    UWORD(190569292),UWORD(214481126),UWORD(241265379),UWORD(271248950),UWORD(304801365),UWORD(342325709),
    UWORD(384276336),UWORD(431149389),UWORD(483502844),UWORD(541946240),UWORD(607163746),UWORD(679903203),
    UWORD(761002156),UWORD(851376628),UWORD(952050665),UWORD(1064144451),UWORD(1188908248),UWORD(1327710076),
    UWORD(1482074143),UWORD(1653668665),UWORD(1844349560),UWORD(2056148051),UWORD(2291320912),
    UWORD(2552338241),UWORD(2841940500),UWORD(3163127352),UWORD(3519222692),UWORD(3913864295)
};

/* To compute p(n) mod 2^64. */
void
_partitions_vec_ui(nn_ptr v, slong len)
{
    slong i, j, n;
    ulong p;

    for (n = 0; n < FLINT_MIN(len, NUMBER_OF_SMALL_PARTITIONS); n++)
        v[n] = partitions_lookup[n];

    for (n = NUMBER_OF_SMALL_PARTITIONS; n < len; n++)
    {
        p = 0;
        for (i = 1, j = 1; j <= n; j += 3 * i + 1, i++)
            p = v[n - j] - p;
        if ((i & 1) == 0)
            p = -p;
        for (i = 1, j = 2; j <= n; j += 3 * i + 2, i++)
            p = v[n - j] - p;
        if ((i & 1) != 0)
            p = -p;
        v[n] = p;
    }
}

/* x (an exact integer) as an fmpz */
static void
_fmpz_set_mp_real_int(fmpz_t res, const mp_real_t x)
{
    slong sh;
    mpz_ptr z;

    if (x->size == 0)
    {
        fmpz_zero(res);
        return;
    }

    sh = x->exp - x->size;     /* the low zero limbs */
    FLINT_ASSERT(sh >= 0 && x->err == 0);
    if (x->exp > INT_MAX)
        flint_throw(FLINT_ERROR, "partitions: p(n) too large for an fmpz\n");
    z = _fmpz_promote(res);
    if (z->_mp_alloc < x->exp)
        mpz_realloc(z, x->exp);
    flint_mpn_zero(z->_mp_d, sh);
    flint_mpn_copyi(z->_mp_d + sh, x->d, x->size);
    z->_mp_size = x->negative ? -x->exp : x->exp;
    _fmpz_demote_val(res);
}

/* The table, the recurrence while p(n) fits a word (n < 417 on 64-bit),
   then the formula (mp_real_partitions_hrr, which repeats the first two
   steps for its own small arguments). */
void
partitions_fmpz_ui(fmpz_t res, ulong n)
{
    if (n < NUMBER_OF_SMALL_PARTITIONS)
    {
        fmpz_set_ui(res, partitions_lookup[n]);
    }
    else if (FLINT_BITS == 64 && n < 417)
    {
        nn_ptr tmp = flint_malloc((n + 1) * sizeof(ulong));
        _partitions_vec_ui(tmp, n + 1);
        fmpz_set_ui(res, tmp[n]);
        flint_free(tmp);
    }
    else
    {
        mp_real_t x;
        mp_real_init(x);
        mp_real_partitions_hrr(x, 0, n);
        _fmpz_set_mp_real_int(res, x);
        mp_real_clear(x);
    }
}

void
partitions_fmpz_fmpz(fmpz_t res, const fmpz_t n, int FLINT_UNUSED(use_doubles))
{
    if (fmpz_sgn(n) < 0)
    {
        fmpz_zero(res);
    }
    else if (fmpz_abs_fits_ui(n))
    {
        partitions_fmpz_ui(res, fmpz_get_ui(n));
    }
    else if (fmpz_bits(n) <= 2 * FLINT_BITS)
    {
        mp_real_t x;
        ulong hi, lo;
        mp_real_init(x);
        fmpz_get_uiui(&hi, &lo, n);
        mp_real_partitions_hrr(x, hi, lo);
        _fmpz_set_mp_real_int(res, x);
        mp_real_clear(x);
    }
    else
    {
        flint_throw(FLINT_ERROR, "partitions_fmpz_fmpz: n too large\n");
    }
}
