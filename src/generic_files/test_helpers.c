/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"

double _flint_test_multiplier = -1.0;

double flint_test_multiplier(void)
{
    if (_flint_test_multiplier == -1.0)
    {
        const char * s = getenv("FLINT_TEST_MULTIPLIER");

        if (s == NULL)
        {
            _flint_test_multiplier = 1.0;
        }
        else
        {
            _flint_test_multiplier = strtod(s, NULL);

            if (!(_flint_test_multiplier >= 0.0 && _flint_test_multiplier <= 1000.0))
                _flint_test_multiplier = 1.0;
        }
    }

    return _flint_test_multiplier;
}
