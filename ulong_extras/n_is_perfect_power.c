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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/

#include <stdlib.h>
#include <flint.h>
#include <ulong_extras.h>


/*
check for the bit size of mp_limb_t, which is equal to FLINT_BITS,
q[i] - 1 congruent class modulo q[i] (about q[i] line no. 105)
are stored as follows:

qclass[i][j]th word contain FLINT_BITS Number of bit and
kth bit (from LSB To MSB) represents the,
(say y = FLINT_BITS * max(j - 1, 0) + k), y mod q[i] congruent
class of q[i] and it is set to 1 if y satisfy the property mentioned
in the line no. 106 otherwise it is set to 0.

for total q[i] congruent class modulo q[i],
    q[i] / FLINT_BITS + (q[i] % FLINT_BITS) ? 0 : 1;
    word are used, i.e. the size of qclass[i]th row.
*/

#if FLINT_BITS == 32
    mp_limb_t qclass[18][25] = {{ 0x3 }, { 0x43 }, { 0x403 }, { 0x10021003 }
                               , { 0x400003 },{ 0x40800002 , 0x100001 },
                               { 0x2 , 0x300c000 , 0x0 , 0x41 },
                               { 0x82 , 0x20080 , 0x40000 , 0x2000 ,
                               0x1004000 , 0x41000001 },
                               { 0x2 , 0x4001 }, { 0x2 , 0x4000001 },
                               { 0x42 , 0x100010 , UWORD(0x80000000) , 0x0
                               , 0x0 , 0x0 , 0x1000000 , 0x0 , 0x80008 ,
                               0x420001 }, { 0x2 , 0x1000 , 0x0 , 0x200 ,
                               0x100001 }, { 0x2 , 0x0 , 0x40001 },
                               { 0x2 , 0x0 , 0x20010000 , 0x0 , 0x0 ,
                               0x1001 }, { 0x2 , 0x3000 , 0x0 , 0x0 , 0x0 ,
                               0x0 , 0x0 , 0xc000 , 0x4000001 },
                               { 0x2 , 0x0 , 0x0 , 0x401 }, { 0x2 , 0x0 ,
                               0x8000001 , 0x0 , 0x0 , 0x8000000 , 0x0 ,
                               0x18 , 0x0 , 0x0 , 0x0 , 0x0 , 0x0 , 0x0 ,
                               0x0 , 0x6 , 0x400 , 0x0 , 0x0 , 0x420 , 0x0
                               , 0x0 , 0x11 }, { 0x2 , 0x0 , 0x180000 ,
                               0x0 , 0x0 , 0x0 , 0x0 , 0x0 , 0x18000000 ,
                               0x0 , 0x0 , 0x4001 }
                              };
#else
    mp_limb_t qclass[18][15] = { { 0x3 }, { 0x43 }, { 0x403 },
                                  { 0x10021003 } , { 0x400003 },
                                  { UWORD(0x10000040800003) },
                                  { UWORD(0x300c00000000002) ,
                                  UWORD(0x4000000001) },
                                  { UWORD(0x2008000000082) ,
                                  UWORD(0x200000040000) ,
                                  UWORD(0x4100000001004001) },
                                  { UWORD(0x400000000003) },
                                  { UWORD(0x400000000000003) },
                                  { UWORD(0x10001000000042)
                                  , UWORD(0x80000000) , 0x0 , 0x1000000 ,
                                  UWORD(0x42000000080009) },
                                  { UWORD(0x100000000002) ,
                                  UWORD(0x20000000000) , 0x100001 },{ 0x2 ,
                                  0x40001 }, { 0x2 , 0x20010000 ,
                                  UWORD(0x100000000001) }, {
                                  UWORD(0x300000000002) , 0x0 ,
                                  0x0 , UWORD(0xc00000000000) , 0x4000001 },
                                  { 0x2 , UWORD(0x40000000001) }, { 0x2 ,
                                  UWORD(0x108000000) ,
                                  UWORD(0x800000000000000) ,
                                  UWORD(0x1800000000) , 0x0 , 0x0 , 0x0 ,
                                  UWORD(0x600000000) , 0x400 ,
                                  UWORD(0x42000000000) ,
                                  0x0 , 0x11 },  { 0x2 , 0x180000 , 0x0 ,
                                  0x0 , 0x18000000 , UWORD(0x400000000001) }
                                };
#endif


int n_is_perfect_power(mp_limb_t * r, mp_limb_t x)
{
    /*
       p[i], is prime less than 64.
       q[i], is prime such that q[i] mod p[i] = 1.
       if x = y ^ p[i] , for some y than,
       x ^ ((q[i] - 1) / p[i]) mod q[i] <= 1.

       see this following link for info:

       ("http://citeseerx.ist.psu.edu/viewdoc/
         download?doi=10.1.1.108.458&rep=rep1&type=pdf")

       z[i] = maximum k such that i = j ^ k for some positive integer j.
       v[i] = minimum j such that i = j ^ k for some positive integer k.
    */

    int p[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
               31, 37, 41, 43, 47, 53, 59, 61};
    int q[] = {3, 7, 11, 29, 23, 53, 103, 191, 47, 59,
               311, 149, 83, 173, 283, 107, 709, 367};
    int z[] = {1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
               1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1};
    int v[] = {0, 1, 2, 3, 2, 5, 6, 7, 2, 3,
               10, 11, 12, 13, 14, 15, 2, 17, 18, 19, 20};

    int i, power = 1, pos, k;
    mp_limb_t prev;

    for (i = 0; i < 18; i++)
    {
        if (x <= 20)
        {
            power *= z[x];
            x = v[x];
            break;
        }

        /* k is congruent class of x mod q[i] */
        k = x % q[i];
        /* pos is the position of above class, in the table*/
        pos = k / FLINT_BITS + (k % FLINT_BITS) ? 0 : 1;

        if (pos == 0) pos = 1; /* in case q[i] divides x */
        /* check if bit corresponding to above position is set */
        if ((qclass[i][pos - 1] & (1 << (k % FLINT_BITS))))
        {
            prev = x;
            x = n_rootrem(r, x, p[i]);  /* test if x is p[i]th power */
            if (*r == 0)
            {
                while (*r == 0)    /* test also for subsequent power of p[i] */
                {
                    power = power * p[i];
                    prev = x;
                    x =  n_rootrem(r, x, p[i]);
                }
                x = prev;
            }
            else x = prev;
        }
    }

    *r = x;       /* store root in r */
    return (power == 1 ? 0 : power);  /* return power */
}
