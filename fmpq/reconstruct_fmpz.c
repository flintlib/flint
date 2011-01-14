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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"


int
_fmpq_reconstruct_fmpz(fmpz_t num, fmpz_t den,
    const fmpz_t a, const fmpz_t m)
{
    int success;

    fmpz_t tmp, q, T0, T1, T2, U0, U1, U2, V0, V1, V2;
    fmpz s;

    if (fmpz_is_zero(a) || fmpz_is_one(a))
    {
        fmpz_set(num, a);
        fmpz_set_ui(den, 1UL);
        return 1;
    }

    fmpz_init(U0); fmpz_init(U1); fmpz_init(U2);
    fmpz_init(V0); fmpz_init(V1); fmpz_init(V2);
    fmpz_init(T0); fmpz_init(T1); fmpz_init(T2);
    fmpz_init(tmp); fmpz_init(q);

    fmpz_set_ui(U0, 1UL);
    fmpz_set_ui(U1, 0UL);
    fmpz_set(U2, m);

    fmpz_set_ui(V0, 0UL);
    fmpz_set_ui(V1, 1UL);
    fmpz_set(V2, a);

    while (1)
    {
        fmpz_mul(tmp, V2, V2);
        fmpz_mul_2exp(tmp, tmp, 1);

        if (fmpz_cmp(tmp, m) <= 0)
            break;

        fmpz_fdiv_q(q, U2, V2);

        fmpz_mul(tmp, q, V0);
        fmpz_sub(T0, U0, tmp);

        fmpz_mul(tmp, q, V1);
        fmpz_sub(T1, U1, tmp);

        fmpz_mul(tmp, q, V2);
        fmpz_sub(T2, U2, tmp);

        s = *U0; *U0 = *V0; *V0 = *T0; *T0 = s;
        s = *U1; *U1 = *V1; *V1 = *T1; *T1 = s;
        s = *U2; *U2 = *V2; *V2 = *T2; *T2 = s;
    }

    fmpz_abs(den, V1);
    fmpz_set(num, V2);

    fmpz_mul(tmp, den, den);
    fmpz_mul_2exp(tmp, tmp, 1);

    if (fmpz_cmp(tmp, m) <= 0)
    {
        fmpz_gcd(tmp, num, den);
        success = (fmpz_cmp_ui(tmp, 1UL) == 0) ? 1 : 0;
    }
    else
    {
        success = 0;
    }

    if (fmpz_sgn(V1) < 0)
        fmpz_neg(num, num);

    fmpz_clear(U0); fmpz_clear(U1); fmpz_clear(U2);
    fmpz_clear(V0); fmpz_clear(V1); fmpz_clear(V2);
    fmpz_clear(T0); fmpz_clear(T1); fmpz_clear(T2);
    fmpz_clear(tmp); fmpz_clear(q);

    return success;
}

int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m)
{
    return _fmpq_reconstruct_fmpz(&res->num, &res->den, a, m);
}
