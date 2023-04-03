/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void _fmpq_ball_init(_fmpq_ball_t x)
{
    fmpz_init(x->left_num);
    fmpz_init(x->left_den);
    fmpz_init(x->right_num);
    fmpz_init(x->right_den);
    x->exact = 0;
}

void _fmpq_ball_clear(_fmpq_ball_t x)
{
    fmpz_clear(x->left_num);
    fmpz_clear(x->left_den);
    fmpz_clear(x->right_num);
    fmpz_clear(x->right_den);
}

