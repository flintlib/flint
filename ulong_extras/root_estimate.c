/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include "flint.h"
#include "ulong_extras.h"
#include "longlong.h"


/* this table contains the value of UWORD_MAX / n, for n in range [1, FLINT_BITS] */

static const mp_limb_t mul_factor[] = {
#ifdef FLINT64    
                        UWORD(0), UWORD_MAX,        UWORD(9223372036854775807),
                        UWORD(6148914691236517205), UWORD(4611686018427387903), 
                        UWORD(3689348814741910323), UWORD(3074457345618258602), 
                        UWORD(2635249153387078802), UWORD(2305843009213693951), 
                        UWORD(2049638230412172401), UWORD(1844674407370955161), 
                        UWORD(1676976733973595601), UWORD(1537228672809129301), 
                        UWORD(1418980313362273201), UWORD(1317624576693539401), 
                        UWORD(1229782938247303441), UWORD(1152921504606846975), 
                        UWORD(1085102592571150095), UWORD(1024819115206086200), 
                        UWORD(970881267037344821), UWORD(922337203685477580), 
                        UWORD(878416384462359600), UWORD(838488366986797800), 
                        UWORD(802032351030850070), UWORD(768614336404564650), 
                        UWORD(737869762948382064), UWORD(709490156681136600), 
                        UWORD(683212743470724133), UWORD(658812288346769700), 
                        UWORD(636094623231363848), UWORD(614891469123651720), 
                        UWORD(595056260442243600), UWORD(576460752303423487), 
                        UWORD(558992244657865200), UWORD(542551296285575047), 
                        UWORD(527049830677415760), UWORD(512409557603043100), 
                        UWORD(498560650640798692), UWORD(485440633518672410), 
                        UWORD(472993437787424400), UWORD(461168601842738790), 
                        UWORD(449920587163647600), UWORD(439208192231179800), 
                        UWORD(428994048225803525), UWORD(419244183493398900), 
                        UWORD(409927646082434480), UWORD(401016175515425035), 
                        UWORD(392483916461905353), UWORD(384307168202282325), 
                        UWORD(376464164769582686), UWORD(368934881474191032), 
                        UWORD(361700864190383365), UWORD(354745078340568300), 
                        UWORD(348051774975651917), UWORD(341606371735362066), 
                        UWORD(335395346794719120), UWORD(329406144173384850), 
                        UWORD(323627089012448273), UWORD(318047311615681924), 
                        UWORD(312656679215416129), UWORD(307445734561825860), 
                        UWORD(302405640552615600), UWORD(297528130221121800), 
                        UWORD(292805461487453200), UWORD(288230376151711743)
#else
                        UWORD(0), UWORD_MAX, UWORD(2147483647), UWORD(1431655765), 
                        UWORD(1073741823), UWORD(858993459), UWORD(715827882), 
                        UWORD(613566756), UWORD(536870911), UWORD(477218588), 
                        UWORD(429496729), UWORD(390451572), UWORD(357913941), 
                        UWORD(330382099), UWORD(306783378), UWORD(286331153), 
                        UWORD(268435455), UWORD(252645135), UWORD(238609294), 
                        UWORD(226050910), UWORD(214748364), UWORD(204522252), 
                        UWORD(195225786), UWORD(186737708), UWORD(178956970), 
                        UWORD(171798691), UWORD(165191049), UWORD(159072862), 
                        UWORD(153391689), UWORD(148102320), UWORD(143165576), 
                        UWORD(138547332), UWORD(134217727), 
#endif
                                                            };
                        /* this table consists of 65 values in case of FLINT64,
                           otherwise 33 */

/* function to get a good approximation of the cube root */
/* Algorithm for this approximation is mentioned in this article */
/* https://en.wikipedia.org/wiki/Fast_inverse_square_root */
/* Intead of the inverse square root, we calculate the nth root */

mp_limb_t
n_root_estimate(double a, int n)
{ 
    typedef union { 
        slong      uword_val;
#ifdef FLINT64
        double     double_val;
#else
        float      double_val;
#endif
    } uni;

    uni alias;
    ulong i, hi, lo, s;
    alias.double_val = a;

#ifdef FLINT64
    s = ((1 << 10) - 1);
    s <<= 52;
#else
    s = ((1 << 7) - 1);
    s <<= 23;
#endif

    i = alias.uword_val;
    i -= s;
    umul_ppmm(hi, lo, i, mul_factor[n]);
    i = hi;
    i += s;
    alias.uword_val = i;
    return (mp_limb_t)alias.double_val;
}
