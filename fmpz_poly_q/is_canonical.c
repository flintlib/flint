/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include "fmpz_poly_q.h"

int fmpz_poly_q_is_canonical(const fmpz_poly_q_t op)
{
    int ans;
    fmpz_poly_t t;

    if (fmpz_poly_is_zero(op->den))
        return 0;

    if (fmpz_sgn(fmpz_poly_lead(op->den)) < 0)
        return 0;

    fmpz_poly_init(t);
    fmpz_poly_gcd(t, op->num, op->den);
    ans = fmpz_poly_is_one(t);
    fmpz_poly_clear(t);
    return ans;
}

