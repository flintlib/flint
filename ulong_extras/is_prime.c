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

    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Dana Jacobsen

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int n_is_prime(mp_limb_t n)
{
    /* flint's "BPSW" checked against Feitsma and Galway's database [1, 2] 
       up to 2^64 by Dana Jacobsen.
       [1]  http://www.janfeitsma.nl/math/psp2/database
       [2]  http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html
    */

    if (n < 11) {
        if (n == 2 || n == 3 || n == 5 || n == 7)   return 1;
        else                                        return 0;
    }
    if (!(n%2) || !(n%3) || !(n%5) || !(n%7))       return 0;
    if (n <  121) /* 11*11 */                       return 1;
    if (!(n%11) || !(n%13) || !(n%17) || !(n%19) ||
        !(n%23) || !(n%29) || !(n%31) || !(n%37) ||
        !(n%41) || !(n%43) || !(n%47) || !(n%53))   return 0;
    if (n < 3481) /* 59*59 */                       return 1;

    return n_is_probabprime(n);
}
