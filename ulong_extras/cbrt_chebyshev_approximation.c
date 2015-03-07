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

double poly1[] = { 0.366434, 1.381445, -1.561725, 1.225528, -0.419472}; /* range : [0.5000, 0.5625] */
double poly2[] = { 0.380314, 1.282490, -1.296988, 0.910552, -0.278851}; /* range : [0.5625, 0.6250] */
double poly3[] = { 0.393245, 1.199557, -1.097430, 0.697022, -0.193126}; /* range : [0.6250, 0.6875] */
double poly4[] = { 0.405375, 1.128860, -0.942848, 0.546736, -0.138312}; /* range : [0.6875, 0.7500] */
double poly5[] = { 0.416818, 1.067743, -0.820391, 0.437649, -0.101857}; /* range : [0.7500, 0.8125] */
double poly6[] = { 0.427663, 1.014282, -0.721542, 0.356392, -0.076801}; /* range : [0.8125, 0.8750] */
double poly7[] = { 0.437984, 0.967050, -0.640465, 0.294520, -0.059090}; /* range : [0.8750, 0.9375] */
double poly8[] = { 0.447840, 0.924961, -0.573045, 0.246510, -0.046267}; /* range : [0.9375, 1.0000] */

/* Coefficients of Chebyshev's approximation polynomial (deg 2) {c0, c1, c2} 
   splitting 0.5 to 1 into 8 equal intervals 

   double poly1[] = {-0.319161, 0.847528, 0.449730}; range : [0.5000, 0.5625] 
   double poly2[] = {-0.265106, 0.786851, 0.466761}; range : [0.5625, 0.6250] 
   double poly3[] = {-0.224345, 0.735992, 0.482629}; range : [0.6250, 0.6875] 
   double poly4[] = {-0.192763, 0.692631, 0.497515}; range : [0.6875, 0.7500] 
   double poly5[] = {-0.167739, 0.655143, 0.511557}; range : [0.7500, 0.8125] 
   double poly6[] = {-0.147537, 0.622349, 0.524867}; range : [0.8125, 0.8750] 
   double poly7[] = {-0.130965, 0.593375, 0.537533}; range : [0.8750, 0.9375] 
   double poly8[] = {-0.117183, 0.567555, 0.549628}; range : [0.9375, 1.0000] */

/* Coefficients of Chebyshev's approximation polynomial (deg 4) [c0, c1, c2, c3, c4, c5] 
   splitting 0.5 to 1 into 4 equal intervals 

   double poly1[] = { 0.372988, 1.332826, -1.42658, 1.058728, -0.342333}; range : [0.500, 0.625]
   double poly2[] = { 0.399034, 1.164711, -1.01881, 0.618246, -0.163540}; range : [0.625, 0.750]
   double poly3[] = { 0.422032, 1.041334, -0.77025, 0.395358, -0.088486}; range : [0.750, 0.875]
   double poly4[] = { 0.442748, 0.946225, -0.60633, 0.269664, -0.052305}; range : [0.875, 1.000] */

mp_limb_t
n_cbrt_chevyshef_poly(mp_limb_t n)
{
    int expo, rem, i, mul;
    double val, factor, root, dec;
    double* table;
    mp_limb_t ret, upper_limit;
    
    /* upper_limit is the max cube root possible for one word */

    upper_limit = 1626;         /* 1626 < (2^32)^(1/3) */
#if FLINT64
    upper_limit = 2642245;      /* 2642245 < (2^64)^(1/3) */
#endif

    val = (double)n;
    dec = frexp(val, &expo);    /* val = dec * 2^expo */
    rem = expo%3;
    expo/=3;                    /* cube root of 2^expo */

    factor = 1;
    if (rem == 1)
        factor = 1.259921;  /* 2^(1/3) */
    if (rem == 2)
        factor = 1.587401;  /* 4^(1/3) */

    /* Selecting correct polynomial based on value of dec */

    table = poly8;

    if (dec < 0.75)
        if (dec < 0.5625)
            table = poly1;
        else if (dec < 0.625)
            table = poly2;
        else if (dec < 0.6875)
            table = poly3;
        else
            table = poly4;
    else 
        if (dec < 0.8125)
            table = poly5;
        else if (dec < 0.8750)
            table = poly6;
        else if (dec < 0.9375)
            table = poly7;
    
    /* Calculating cube root of dec using chebyshev approximation polynomial */
    /* Evaluating approx polynomial at (dec) by Horner's method */

    root = table[4];
    for (i = 3; i > -1; i--)
        root = root * dec + table[i];

    mul = UWORD(1)<<expo;
    root*=mul;           /* dec^(1/3) * 2^(expo/3) */
    root*=factor;        /* root*= (expo%3)^(1/3) */
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
