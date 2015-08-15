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

    Built upon existing FLINT siqs
    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpqs.h"

void
mpqs_sqrtmod_psq(fmpz_t root, fmpz_t n, mp_limb_t prime)
{
    fmpz_t p, psq, a, b, k;
    
    fmpz_init(psq);
    fmpz_init(b);
    fmpz_init(k);
    fmpz_init_set(a, n);
    fmpz_init_set_ui(p, prime);

    fmpz_mod(b, a, p);
    fmpz_sqrtmod(root, a, p);
    fmpz_mul(psq, p, p);

    /* b = root^2 mod p^2 */
    fmpz_mul(b, root, root);
    fmpz_mod(b, b, psq);

    if (fmpz_cmp(a, b) > 0)
        fmpz_sub(k, a, b);
    else
        fmpz_sub(k, b, a);

    fmpz_mod(k, k, psq);
    fmpz_divexact_ui(k, k, prime);

    if (fmpz_cmp(a, b) < 0)
        fmpz_negmod(k, k, p);

    fmpz_mul_2exp(b, root, 1); /* b = root * 2 */
    fmpz_invmod(b, b, p);
    fmpz_mul(k, k, b);
    fmpz_mod(k, k, p);
    fmpz_mul(k, k, p);
    fmpz_add(root, root, k);

    /* if first root is even, second root is odd */
    
    if (fmpz_is_even(root) == 1)
        fmpz_sub(root, psq, root);

    fmpz_clear(psq);
    fmpz_clear(b);
    fmpz_clear(k);
    fmpz_clear(p);
}