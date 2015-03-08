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

/* Coefficients of Chebyshev's approximation polynomial (deg 4) {c0, c1, c2, c3, c4, c5}
   splitting 0.5 to 1 into 8 equal intervals 

   Values of these coefficients of Chebyshev's approximation polynomial have been
   calculated from the python module, "mpmath" - http://mpmath.org/   
   function call:
   mpmath.chebyfit(lambda x: mpmath.root(x,3), [i, j], 5, error=True)
   where (i, j) is the  range                                                          */

/*                                      c0           c1            c2           c3             c4               range */
static const float coeff[8][5] = {{ 0.366434535, 1.381445155, -1.561722404, 1.225528386, -0.4194725256}, /*[0.5000, 0.5625]*/
                                  { 0.380314156, 1.282490068, -1.296988440, 0.910552775, -0.2788510285}, /*[0.5625, 0.6250]*/
                                  { 0.393245231, 1.199557840, -1.097430124, 0.697022553, -0.1931266447}, /*[0.6250, 0.6875]*/
                                  { 0.405375312, 1.128860765, -0.942848589, 0.546736819, -0.1383122321}, /*[0.6875, 0.7500]*/
                                  { 0.416818256, 1.067743101, -0.820391423, 0.437649088, -0.1018574075}, /*[0.7500, 0.8125]*/
                                  { 0.427663998, 1.014282651, -0.721542665, 0.356392064, -0.0768012174}, /*[0.8125, 0.8750]*/
                                  { 0.437984966, 0.967050791, -0.640465903, 0.294520591, -0.0590907229}, /*[0.8750, 0.9375]*/
                                  { 0.447840465, 0.924961634, -0.573045133, 0.246510218, -0.0462671851}};/*[0.9375, 1.0000]*/

static const float factor_table[] = {1.000000, 1.259921, 1.587401};  /* {1^(1/3), 2^(1/3), 4^(1/3)} */

mp_limb_t
n_cbrt_chebyshef_poly(mp_limb_t n)
{
    typedef union { 
        mp_limb_t  uword_val;
        double     double_val;
    } uni;

    int rem, mul;
    double factor, root, dec, dec2;
    mp_limb_t ret, upper_limit, expo, table_index;
    mp_limb_t* n_ptr;
    uni alias;

    /* upper_limit is the max cube root possible for one word */

    upper_limit = 1626;         /* 1626 < (2^32)^(1/3) */
#if FLINT64
    upper_limit = 2642245;      /* 2642245 < (2^64)^(1/3) */
#endif

    alias.double_val = (double)n;

    expo = alias.uword_val & 0x7FF0000000000000;  /* extracting exponent */
    expo >>= 52;
    expo -= 1022;                                 /* Subtracting bias */

    /* extracting first 3 bits of mantissa, this will help select correct poly */
    /* note mantissa of 0.5 is 0x0000000000000 not 0x1000000000000 */

    table_index = alias.uword_val & 0x000F000000000000;
    table_index >>= 49;

    /* extracting decimal part, 0.5 <= dec <= 1 */
    ret = alias.uword_val & 0x000FFFFFFFFFFFFF;
    ret |= 0x3FE0000000000000;
    alias.uword_val = ret;
    dec = alias.double_val;

    rem = expo%3;
    expo/=3;                            /* cube root of 2^expo */
    factor = factor_table[rem];         /* select factor */

    /* Calculating cube root of dec using chebyshev approximation polynomial */
    /* Evaluating approx polynomial at (dec) by Estrin's scheme */
    
    dec2 = dec*dec;
    root = (coeff[table_index][0] + coeff[table_index][1]*dec);
    root += ((coeff[table_index][2] + coeff[table_index][3]*dec + coeff[table_index][4]*dec2)*dec2);

    mul = UWORD(1) << expo;       /* mul = 2^expo */
    root *= mul;                  /* dec^(1/3) * 2^(expo/3) */
    root *= factor;               /* root*= (expo%3)^(1/3) */
    ret = root;

    /* In case ret^3 or (ret+1)^3 will cause overflow */

    if (ret >= upper_limit)      
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