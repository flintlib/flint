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
#include "flint.h"
#include "ulong_extras.h"

/* Coefficients of Chebyshev's approximation polynomial (deg 2) {c0, c1, c2}
   splitting 0.5 to 1 into 8 equal intervals 

   Values of these coefficients of Chebyshev's approximation polynomial have been
   calculated from the python module, "mpmath" - http://mpmath.org/   
   function call:
   mpmath.chebyfit(lambda x: mpmath.root(x,3), [i, j], 3, error=True)
   where (i, j) is the  range                                                           */

static const float factor_table[] = {1.000000, 1.259921, 1.587401};  /* {1^(1/3), 2^(1/3), 4^(1/3)}  */

/*                                       c0          c1             c2               range             */
static const float coeff[16][3] = {{0.445434042, 0.864136635, -0.335205926},    /* [0.50000, 0.53125]  */
                                   {0.454263239, 0.830878907, -0.303884962},    /* [0.53125, 0.56250]  */
                                   {0.462761624, 0.800647514, -0.276997626},    /* [0.56250, 0.59375]  */
                                   {0.470958569, 0.773024522, -0.253724515},    /* [0.59375, 0.62500]  */
                                   {0.478879482, 0.747667468, -0.233429710},    /* [0.62500, 0.65625]  */
                                   {0.486546506, 0.724292830, -0.215613166},    /* [0.65625, 0.68750]  */
                                   {0.493979069, 0.702663686, -0.199877008},    /* [0.68750, 0.71875]  */
                                   {0.501194325, 0.682580388, -0.185901247},    /* [0.71875, 0.75000]  */
                                   {0.508207500, 0.663873398, -0.173426009},    /* [0.75000, 0.78125]  */
                                   {0.515032183, 0.646397742, -0.162238357},    /* [0.78125, 0.81250]  */
                                   {0.521680556, 0.630028647, -0.152162376},    /* [0.81250, 0.84375]  */
                                   {0.528163588, 0.614658092, -0.143051642},    /* [0.84375, 0.87500]  */
                                   {0.534491194, 0.600192044, -0.134783425},    /* [0.87500, 0.90625]  */
                                   {0.540672371, 0.586548233, -0.127254189},    /* [0.90625, 0.93750]  */
                                   {0.546715310, 0.573654340, -0.120376066},    /* [0.93750, 0.96875]  */
                                   {0.552627494, 0.561446514, -0.114074068}};   /* [0.96875, 1.00000]  */
mp_limb_t
n_cbrt_chebyshev_approx(mp_limb_t n)
{
    typedef union { 
        mp_limb_t  uword_val;
#ifdef FLINT64
        double     double_val;
#else
        float      double_val;
#endif
    } uni;

    int rem, mul;
    double factor, root, dec, dec2;
    mp_limb_t ret, expo, table_index;
    uni alias;

    /* upper_limit is the max cube root possible for one word */

#ifdef FLINT64
    const mp_limb_t upper_limit = 2642245;              /* 2642245 < (2^64)^(1/3) */
    const mp_limb_t expo_mask = 0x7FF0000000000000;     /* exponent bits in double */
    const mp_limb_t mantissa_mask = 0x000FFFFFFFFFFFFF; /* mantissa bits in float */
    const mp_limb_t table_mask = 0x000F000000000000;    /* first 4 bits of mantissa */
    const int mantissa_bits = 52;
    const mp_limb_t bias_hex = 0x3FE0000000000000;      
    const int bias = 1022;
    alias.double_val = (double)n;
#else
    const mp_limb_t upper_limit = 1625;         /* 1625 < (2^32)^(1/3) */
    const mp_limb_t expo_mask = 0x7F800000;     /* exponent bits in float */
    const mp_limb_t mantissa_mask = 0x007FFFFF; /* mantissa bits in float */
    const mp_limb_t table_mask = 0x00780000;    /* first 4 bits of mantissa */
    const int mantissa_bits = 23;
    const mp_limb_t bias_hex = 0x3F000000;      
    const int bias = 126;
    alias.double_val = (float)n;
#endif


    expo = alias.uword_val & expo_mask;  /* extracting exponent */
    expo >>= mantissa_bits;
    expo -= bias;                        /* Subtracting bias */

    /* extracting first 4 bits of mantissa, this will help select correct poly */
    /* note mantissa of 0.5 is 0x0000000000000 not 0x1000000000000 */

    table_index = alias.uword_val & table_mask;
    table_index >>= (mantissa_bits - 4);

    /* extracting decimal part, 0.5 <= dec <= 1 */
    ret = alias.uword_val & mantissa_mask;
    ret |= bias_hex;
    alias.uword_val = ret;
    dec = alias.double_val;

    rem = expo % 3;
    expo /= 3;                          /* cube root of 2^expo */
    factor = factor_table[rem];         /* select factor */

    /* Calculating cube root of dec using chebyshev approximation polynomial */
    /* Evaluating approx polynomial at (dec) by Estrin's scheme */
    
    dec2 = dec*dec;
    root = (coeff[table_index][0] + coeff[table_index][1] * dec);
    root += (coeff[table_index][2] * dec2);

    mul = UWORD(1) << expo;       /* mul = 2^expo */
    root *= mul;                  /* dec^(1/3) * 2^(expo/3) */
    root *= factor;               /* root*= (expo%3)^(1/3) */
    ret = root;

    /* In case ret^3 or (ret+1)^3 will cause overflow */

    if (ret >= upper_limit)      
    {
        if (n >= upper_limit * upper_limit * upper_limit)
            return upper_limit;
        ret = upper_limit - 1;
    }
    while (ret * ret * ret <= n)
    {
        (ret) += 1;
        if (ret == upper_limit)
            break;
    }
    while (ret * ret * ret > n)
        (ret) -= 1;

    return ret;
}
