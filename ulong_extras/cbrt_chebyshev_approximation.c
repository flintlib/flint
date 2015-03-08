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

/* Coefficients of Chebyshev's approximation polynomial (deg 4) {c0, c1, c2, c3, c4, c5}
   splitting 0.5 to 1 into 8 equal intervals */

float poly1[] = { 0.366434535, 1.381445155, -1.561722404, 1.225528386, -0.4194725256}; /* range : [0.5000, 0.5625] */
float poly2[] = { 0.380314156, 1.282490068, -1.296988440, 0.910552775, -0.2788510285}; /* range : [0.5625, 0.6250] */
float poly3[] = { 0.393245231, 1.199557840, -1.097430124, 0.697022553, -0.1931266447}; /* range : [0.6250, 0.6875] */
float poly4[] = { 0.405375312, 1.128860765, -0.942848589, 0.546736819, -0.1383122321}; /* range : [0.6875, 0.7500] */
float poly5[] = { 0.416818256, 1.067743101, -0.820391423, 0.437649088, -0.1018574075}; /* range : [0.7500, 0.8125] */
float poly6[] = { 0.427663998, 1.014282651, -0.721542665, 0.356392064, -0.0768012174}; /* range : [0.8125, 0.8750] */
float poly7[] = { 0.437984966, 0.967050791, -0.640465903, 0.294520591, -0.0590907229}; /* range : [0.8750, 0.9375] */
float poly8[] = { 0.447840465, 0.924961634, -0.573045133, 0.246510218, -0.0462671851}; /* range : [0.9375, 1.0000] */

float* poly_table[] = { poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8};
float factor_table[] = {1.000000, 1.259921, 1.587401};  /* {1^(1/3), 2^(1/3), 4^(1/3)} */

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