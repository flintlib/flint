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

    Copyright (C) 2009 Tom Boothby
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

#if FLINT64

mp_limb_t FLINT_ODD_PRIME_LOOKUP[] =
{
    0x816d129a64b4cb6eUL, 0x2196820d864a4c32UL, 0xa48961205a0434c9UL,
    0x4a2882d129861144UL, 0x834992132424030UL,  0x148a48844225064bUL,
    0xb40b4086c304205UL,  0x65048928125108a0UL, 0x80124496804c3098UL,
    0xc02104c941124221UL, 0x804490000982d32UL,  0x220825b082689681UL,
    0x9004265940a28948UL, 0x6900924430434006UL, 0x12410da408088210UL,
    0x86122d22400c060UL,  0x110d301821b0484UL,  0x14916022c044a002UL,
    0x92094d204a6400cUL,  0x4ca2100800522094UL, 0xa48b081051018200UL,
    0x34c108144309a25UL,  0x2084490880522502UL, 0x241140a218003250UL,
    0xa41a00101840128UL,  0x2926000836004512UL, 0x10100480c0618283UL,
    0xc20c26584822006dUL, 0x4520582024894810UL, 0x10c0250219002488UL,
    0x802832ca01140868UL, 0x60901300264b0400UL
};
 
#else

mp_limb_t FLINT_ODD_PRIME_LOOKUP[] =
{
    0x64b4cb6eUL, 0x816d129aUL, 0x864a4c32UL, 0x2196820dUL, 
    0x5a0434c9UL, 0xa4896120UL, 0x29861144UL, 0x4a2882d1UL, 
    0x32424030UL, 0x8349921UL,  0x4225064bUL, 0x148a4884UL,
    0x6c304205UL, 0xb40b408UL,  0x125108a0UL, 0x65048928UL, 
    0x804c3098UL, 0x80124496UL, 0x41124221UL, 0xc02104c9UL, 
    0x982d32UL,   0x8044900UL,  0x82689681UL, 0x220825b0UL,
    0x40a28948UL, 0x90042659UL, 0x30434006UL, 0x69009244UL, 
    0x8088210UL,  0x12410da4UL, 0x2400c060UL, 0x86122d2UL,  
    0x821b0484UL, 0x110d301UL,  0xc044a002UL, 0x14916022UL,
    0x4a6400cUL,  0x92094d2UL,  0x522094UL,   0x4ca21008UL, 
    0x51018200UL, 0xa48b0810UL, 0x44309a25UL, 0x34c1081UL,  
    0x80522502UL, 0x20844908UL, 0x18003250UL, 0x241140a2UL,
    0x1840128UL,  0xa41a001UL,  0x36004512UL, 0x29260008UL, 
    0xc0618283UL, 0x10100480UL, 0x4822006dUL, 0xc20c2658UL, 
    0x24894810UL, 0x45205820UL, 0x19002488UL, 0x10c02502UL,
    0x1140868UL,  0x802832caUL, 0x264b0400UL, 0x60901300UL
};

#endif

int n_is_oddprime_small(mp_limb_t n) 
{
    mp_limb_t q = n / 2;
    mp_limb_t x = (q & (FLINT_BITS - 1UL));
    return (FLINT_ODD_PRIME_LOOKUP[q / FLINT_BITS] & (1UL << x)) >> x;
}
