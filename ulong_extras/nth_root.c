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

    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

******************************************************************************/


#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

/* A table of precomputed inverses of values from 1 to 64  
   inv_table[n] = i/n for all n in range[1, 64]  
   inv_table[0] is set to 0, albeit it will never be called */

static const double inv_table[] = { 
                       0.000000000000000, 1.000000000000000, 0.500000000000000, 
                       0.333333333333333, 0.250000000000000, 0.200000000000000, 
                       0.166666666666667, 0.142857142857143, 0.125000000000000, 
                       0.111111111111111, 0.100000000000000, 0.090909090909091, 
                       0.083333333333333, 0.076923076923077, 0.071428571428571, 
                       0.066666666666667, 0.062500000000000, 0.058823529411765, 
                       0.055555555555556, 0.052631578947368, 0.050000000000000, 
                       0.047619047619048, 0.045454545454545, 0.043478260869565, 
                       0.041666666666667, 0.040000000000000, 0.038461538461538, 
                       0.037037037037037, 0.035714285714286, 0.034482758620690, 
                       0.033333333333333, 0.032258064516129, 0.031250000000000, 
                       0.030303030303030, 0.029411764705882, 0.028571428571429, 
                       0.027777777777778, 0.027027027027027, 0.026315789473684, 
                       0.025641025641026, 0.025000000000000, 0.024390243902439, 
                       0.023809523809524, 0.023255813953488, 0.022727272727273, 
                       0.022222222222222, 0.021739130434783, 0.021276595744681, 
                       0.020833333333333, 0.020408163265306, 0.020000000000000, 
                       0.019607843137255, 0.019230769230769, 0.018867924528302, 
                       0.018518518518519, 0.018181818181818, 0.017857142857143, 
                       0.017543859649123, 0.017241379310345, 0.016949152542373, 
                       0.016666666666667, 0.016393442622951, 0.016129032258065, 
                       0.015873015873016, 0.015625000000000 };

                       /* this table consists 65 values */

/* This table has the max possible base for a given root. 
   max_base[n] = UWORD_MAX^(1/n) for n in range [1, FLINT_BITS] 
   max_base[0] is set to 0, although it will never be called */

static const mp_limb_t max_base[] = { 
#if FLINT64
                        0, UWORD_MAX, 4294967296, 2642245 ,65536,
                        7131, 1625, 565, 256, 138, 84, 56, 40, 30,
                        23, 19, 16, 13, 11, 10, 9, 8, 7, 6, 6, 5,
                        5, 5, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 
                        3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 
#else
                        0, 4294967295,  65535, 1625, 255,  84,  40,
                        23, 15, 11, 9, 7, 6, 5, 4, 4, 3, 3, 3, 3, 3,
                        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
#endif
                                                        };

                        /* this table consists of 65 values in case of FLINT64,
                           otherwise 33 */

/* this table contains the value of UWORD_MAX / n, for n in range [1, FLINT_BITS] */

static const mp_limb_t mul_factor[] = {
#if FLINT64    
                        0, 18446744073709551615, 9223372036854775807, 6148914691236517205, 
                        4611686018427387903, 3689348814741910323, 3074457345618258602, 
                        2635249153387078802, 2305843009213693951, 2049638230412172401, 
                        1844674407370955161, 1676976733973595601, 1537228672809129301, 
                        1418980313362273201, 1317624576693539401, 1229782938247303441, 
                        1152921504606846975, 1085102592571150095, 1024819115206086200, 
                        970881267037344821, 922337203685477580, 878416384462359600, 
                        838488366986797800, 802032351030850070, 768614336404564650, 
                        737869762948382064, 709490156681136600, 683212743470724133, 
                        658812288346769700, 636094623231363848, 614891469123651720, 
                        595056260442243600, 576460752303423487, 558992244657865200, 
                        542551296285575047, 527049830677415760, 512409557603043100, 
                        498560650640798692, 485440633518672410, 472993437787424400, 
                        461168601842738790, 449920587163647600, 439208192231179800, 
                        428994048225803525, 419244183493398900, 409927646082434480, 
                        401016175515425035, 392483916461905353, 384307168202282325, 
                        376464164769582686, 368934881474191032, 361700864190383365, 
                        354745078340568300, 348051774975651917, 341606371735362066, 
                        335395346794719120, 329406144173384850, 323627089012448273, 
                        318047311615681924, 312656679215416129, 307445734561825860, 
                        302405640552615600, 297528130221121800, 292805461487453200, 
                        288230376151711743
#else
                        0, 4294967295, 2147483647, 1431655765, 1073741823, 858993459, 
                        715827882, 613566756, 536870911, 477218588, 429496729, 390451572, 
                        357913941, 330382099, 306783378, 286331153, 268435455, 252645135, 
                        238609294, 226050910, 214748364, 204522252, 195225786, 186737708, 
                        178956970, 171798691, 165191049, 159072862, 153391689, 148102320, 
                        143165576, 138547332, 134217727, 
#endif
                                                            };
                        /* this table consists of 65 values in case of FLINT64,
                           otherwise 33 */

/* function to get a good approximation of the cube root */
/* Algorithm for this approximation is mentioned in this article */
/* http://en.wikipedia.org/wiki/Fast_inverse_square_root */
/* Intead of the inverse square root, we calculate the nth root */

double
root_estimate(double a, int n)
{ 
    typedef union { 
        slong      uword_val;
#if FLINT64
        double     double_val;
#else
        float      double_val;
#endif
    } uni;

    uni alias;
    ulong i, hi, lo;
    alias.double_val = a;
    i = alias.uword_val;

#ifdef FLINT64
    slong s = ((1 << 10) - 1);
    s <<= 52;
#else
    slong s = ((1 << 7) - 1);
    s <<= 23;
#endif

    i -= s;
    umul_ppmm(hi, lo, i, mul_factor[n]);
    i = hi;
    i += s;
    alias.uword_val = i;
    return alias.double_val;
}

mp_limb_t
n_root(mp_limb_t n, mp_limb_t root)
{
    mp_limb_t x, currval, base; 
    double dx;

    if (!n || !root)
        return 0;

    if (root == 1)
        return n;

    if (root == 2)
        return n_sqrt(n);

    if (root == 3)
        return n_cbrt(n);

    if (root >= FLINT_BITS || (UWORD(1) << root) > n)
        return 1;

    const mp_limb_t upper_limit = max_base[root];      /* n <= upper_limit^root */

    x =  (mp_limb_t) root_estimate((double)n, root);

    /* one round of newton iteration */

    currval = n_pow(x, root-1);
    dx = n / currval;
    dx -= x;
    dx *= inv_table[root];
    dx = floor(dx);
    x += dx;

    base = x;

    if (base > upper_limit)
        base = upper_limit - 1;

    currval = n_pow(base, root);
    if (currval == n)
        goto final;

    while (currval <= n)
    {
        (base) += 1;
        if (base == upper_limit)
            goto final;
        currval = n_pow(base, root);
    }

    while (currval > n)
    {
        (base) -= 1;
        currval = n_pow(base, root);
    }

    final:
    return base;
}
