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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"
#include "ulong_extras.h"

#define BERNOULLI_DENOM_MAX_SMALL 178

#if FLINT64
#define __u32 unsigned int
#else
#define __u32 mp_limb_t
#endif

static const __u32 __bernoulli_denom_small[] =
{
    1UL, 6UL, 30UL, 42UL, 30UL, 66UL, 2730UL, 6UL, 510UL, 798UL, 330UL,
    138UL, 2730UL, 6UL, 870UL, 14322UL, 510UL, 6UL, 1919190UL, 6UL, 13530UL,
    1806UL, 690UL, 282UL, 46410UL, 66UL, 1590UL, 798UL, 870UL, 354UL,
    56786730UL, 6UL, 510UL, 64722UL, 30UL, 4686UL, 140100870UL, 6UL, 30UL,
    3318UL, 230010UL, 498UL, 3404310UL, 6UL, 61410UL, 272118UL, 1410UL, 6UL,
    4501770UL, 6UL, 33330UL, 4326UL, 1590UL, 642UL, 209191710UL, 1518UL,
    1671270UL, 42UL, 1770UL, 6UL, 2328255930UL, 6UL, 30UL, 4357878UL, 510UL,
    8646UL, 4206930UL, 6UL, 4110UL, 274386UL, 679470UL, 6UL, 2381714790UL,
    6UL, 4470UL, 2162622UL, 30UL, 138UL, 1794590070UL, 6UL, 230010UL,
    130074UL, 2490UL, 1002UL, 3404310UL, 66UL, 5190UL, 2478UL, 1043970UL,
    1074UL,
};

void arith_bernoulli_number_denom(fmpz_t den, ulong n)
{
    len_t i;
    mp_limb_t p;

    if (n % 2 == 1)
    {
        fmpz_set_ui(den, 1 + (n == 1));
    }
    else if (n <= BERNOULLI_DENOM_MAX_SMALL)
    {
        fmpz_set_ui(den, __bernoulli_denom_small[n / 2]);
    }
    else
    {
        n_prime_pi_bounds(&p, &p, n);
        n_compute_primes(p);

        fmpz_set_ui(den, 6UL);
        for (i = 2; i < n; i++)
        {
            p = flint_primes[i];
            if (p - 1 > n)
                break;
            if (n % (p - 1) == 0)
                fmpz_mul_ui(den, den, p);
        }
    }
}
