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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

#if FLINT64
#define NUM_SMALL_FIB 94
#else
#define NUM_SMALL_FIB 48
#endif

static const mp_limb_t small_fib[NUM_SMALL_FIB] =
{
    0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
    1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393,
    196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887,
    9227465, 14930352, 24157817, 39088169, 63245986, 102334155,
    165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903,
    2971215073UL,
#if FLINT64
    4807526976UL, 7778742049UL, 12586269025UL, 20365011074UL, 32951280099UL,
    53316291173UL, 86267571272UL, 139583862445UL, 225851433717UL,
    365435296162UL, 591286729879UL, 956722026041UL, 1548008755920,
    2504730781961UL, 4052739537881UL, 6557470319842UL, 10610209857723,
    17167680177565UL, 27777890035288UL, 44945570212853UL, 72723460248141,
    117669030460994UL, 190392490709135UL, 308061521170129UL, 498454011879264,
    806515533049393UL, 1304969544928657UL, 2111485077978050UL,
    3416454622906707UL, 5527939700884757UL, 8944394323791464UL,
    14472334024676221UL, 23416728348467685UL, 37889062373143906UL,
    61305790721611591UL, 99194853094755497UL, 160500643816367088UL,
    259695496911122585UL, 420196140727489673UL, 679891637638612258UL,
    1100087778366101931UL, 1779979416004714189UL, 2880067194370816120UL,
    4660046610375530309UL, 7540113804746346429UL, 12200160415121876738UL
#endif
};

void fmpz_fib_ui(fmpz_t f, ulong n)
{
    if (n < NUM_SMALL_FIB)
        fmpz_set_ui(f, small_fib[n]);
    else
        mpz_fib_ui(_fmpz_promote(f), n);
}
