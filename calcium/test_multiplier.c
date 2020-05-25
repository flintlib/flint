/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "calcium.h"

double _calcium_test_multiplier = -1.0;

double calcium_test_multiplier()
{
    if (_calcium_test_multiplier == -1.0)
    {
        const char * s = getenv("CALCIUM_TEST_MULTIPLIER");

        if (s == NULL)
        {
            _calcium_test_multiplier = 1.0;
        }
        else
        {
            _calcium_test_multiplier = strtod(s, NULL);

            if (!(_calcium_test_multiplier >= 0.0 && _calcium_test_multiplier <= 1000.0))
                _calcium_test_multiplier = 1.0;
        }
    }

    return _calcium_test_multiplier;
}

