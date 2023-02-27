/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "arb.h"

double _arb_test_multiplier = -1.0;

double arb_test_multiplier()
{
    if (_arb_test_multiplier == -1.0)
    {
        const char * s = getenv("ARB_TEST_MULTIPLIER");

        if (s == NULL)
        {
            _arb_test_multiplier = 1.0;
        }
        else
        {
            _arb_test_multiplier = strtod(s, NULL);

            if (!(_arb_test_multiplier >= 0.0 && _arb_test_multiplier <= 1000.0))
                _arb_test_multiplier = 1.0;
        }
    }

    return _arb_test_multiplier;
}

