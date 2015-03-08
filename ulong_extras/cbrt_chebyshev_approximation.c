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
#include "flint.h"
#include "ulong_extras.h"

/* Values of coefficient Chebyshev's approximation polynomial calculated
   from the python mpmath module : http://mpmath.org/   */

/* Coefficients of Chebyshev's approximation polynomial (deg 4) [c0, c1, c2, c3, c4, c5] 
   splitting 0.5 to 1 into 8 equal intervals */

float poly1[] = { 0.36643453, 1.38144515, -1.56172240, 1.22552838, -0.41947252}; /* range : [0.5000, 0.5625] */
float poly2[] = { 0.38031415, 1.28249006, -1.29698844, 0.91055277, -0.27885102}; /* range : [0.5625, 0.6250] */
float poly3[] = { 0.38031415, 1.28249006, -1.29698844, 0.91055277, -0.27885102}; /* range : [0.6250, 0.6875] */
float poly4[] = { 0.40537531, 1.12886076, -0.94284858, 0.54673681, -0.13831223}; /* range : [0.6875, 0.7500] */
float poly5[] = { 0.41681825, 1.06774310, -0.82039142, 0.43764908, -0.10185740}; /* range : [0.7500, 0.8125] */
float poly6[] = { 0.42766399, 1.01428265, -0.72154266, 0.35639206, -0.07680121}; /* range : [0.8125, 0.8750] */
float poly7[] = { 0.43798496, 0.96705079, -0.64046590, 0.29452059, -0.05909072}; /* range : [0.8750, 0.9375] */
float poly8[] = { 0.44784046, 0.92496163, -0.57304513, 0.24651021, -0.04626718}; /* range : [0.9375, 1.0000] */

float* poly_table[] = { poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8};

float factor_table[] = {1.000000, 1.259921, 1.587401};

mp_limb_t
n_cbrt_chevyshef_poly(mp_limb_t n)
{
    int rem, mul;
    double val, factor, root, dec, dec2, dec4;
    float* table;
    mp_limb_t ret, upper_limit, expo, table_index;
    mp_limb_t* n_ptr;
    
    /* upper_limit is the max cube root possible for one word */

    upper_limit = 1626;         /* 1626 < (2^32)^(1/3) */
#if FLINT64
    upper_limit = 2642245;      /* 2642245 < (2^64)^(1/3) */
#endif

    val = (double)n;

    /* aliasing val as an mp_limb_t to get binary representation of double */
    n_ptr = (mp_limb_t*)&val;
    expo = *n_ptr & 0x7FF0000000000000; /* extracting exponent */
    expo>>=52;
    expo-=1022;     /* Subtracting bias */
    
    /* extracting first 3 bits of mantissa, this will help select correct poly */
    /* note mantissa of 0.5 is 0x0000000000000 not 0x1000000000000 */
    table_index = *n_ptr & 0x000F000000000000;  
    table_index>>=49;

    /* extracting decimal part, 0.5 <= dec <= 1 */
    ret = (*(mp_limb_t*)&val) & 0x000FFFFFFFFFFFFF;
    ret |= 0x3FE0000000000000;
    dec = *(double*)&ret;

    rem = expo%3;
    expo/=3;                            /* cube root of 2^expo */
    factor = factor_table[rem];         /* select factor */
    table = poly_table[table_index];    /* select poly */

    /* Calculating cube root of dec using chebyshev approximation polynomial */
    /* Evaluating approx polynomial at (dec) by Estrin's scheme */
    
    dec2 = dec*dec;
    dec4 = dec2*dec2;
    root = (table[0] + table[1]*dec);
    root+= ((table[2] + table[3]*dec)*dec2);
    root+= (table[4]*dec4);

    mul = UWORD(1)<<expo;   /* mul = 2^expo */
    root*=mul;              /* dec^(1/3) * 2^(expo/3) */
    root*=factor;           /* root*= (expo%3)^(1/3) */
    ret = root;

    /* In case ret^3 or (ret+1)^3 will cause overflow */

    if (ret>= upper_limit)      
    {
        if (n >= upper_limit*upper_limit*upper_limit)
            return upper_limit;
        ret = upper_limit - 1;
    }
    while (ret*ret*ret <= n)
    {
        (ret) += 1;
        if (ret == upper_limit)
            break;
    }
    while (ret*ret*ret > n)
        (ret) -= 1;

    return ret;
}