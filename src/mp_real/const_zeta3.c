/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "flint.h"
#include "mp_real.h"
#include "impl.h"

/* zeta(3) by Zuniga's series (2023-vi),

    zeta(3) = (1/48) sum_{k>=1} (-1)^(k-1) P(k) / ( k^5 (2k-1)^3
              (3k-1)(3k-2)(4k-1)(4k-3)(6k-1)(6k-5)
              C(5k,k) C(5k,2k) C(9k,4k) C(10k,5k) C(12k,6k) ),

   P of degree 11, in the hypergeometric form of y-cruncher's formula
   file (https://mathoverflow.net/questions/454929).  Of the twelve
   files for zeta(3) this is the fastest (cost 2.05 bits per term
   against 2.31 for Zuniga 2023-v; measured 20% faster than the
   Amdeberhan-Zeilberger series used by the old arb code). */
void
_mp_real_const_zeta3_compute(mp_real_t res, slong n)
{
    static const int64_t P[] = {
        INT64_C(-3143448000), INT64_C(156286859400), INT64_C(-3292502315430),
        INT64_C(38721705264979), INT64_C(-282805786014979),
        INT64_C(1352700034136826), INT64_C(-4348596587040104),
        INT64_C(9451223531851808), INT64_C(-13684352515879536),
        INT64_C(12632254526031264), INT64_C(-6719460725627136),
        INT64_C(1565994397644288) };
    static const int64_t Q[] = {
        INT64_C(44008272000), INT64_C(-2334151436400),
        INT64_C(53522442803340), INT64_C(-703273183134030),
        INT64_C(5931859745397870), INT64_C(-34140867105175650),
        INT64_C(139058868850409430), INT64_C(-409481300311614720),
        INT64_C(880500176512163280), INT64_C(-1382139595517666400),
        INT64_C(1565294958171053280), INT64_C(-1244539247650560000),
        INT64_C(658690593528960000), INT64_C(-208277254886400000),
        INT64_C(29753893555200000) };
    static const int64_t R[] = {
        INT64_C(0), INT64_C(0), INT64_C(0), INT64_C(0), INT64_C(0),
        INT64_C(30), INT64_C(-691), INT64_C(6781), INT64_C(-37374),
        INT64_C(127976), INT64_C(-283232), INT64_C(406224),
        INT64_C(-364896), INT64_C(186624), INT64_C(-41472) };

    mp_real_hypgeom_series_int64(res, 1, 1, 0, 48, P, 12, Q, 15, R, 15, n);
}
