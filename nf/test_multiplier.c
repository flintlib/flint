/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nf.h"

long int _antic_test_multiplier = -1;

long int antic_test_multiplier()
{
    if (_antic_test_multiplier == -1)
    {
        const char * s = getenv("ANTIC_TEST_MULTIPLIER");

        if (s == NULL)
        {
            _antic_test_multiplier = 10;
        }
        else
        {
            _antic_test_multiplier = strtol(s, NULL, 10);

            if (!(_antic_test_multiplier >= 0.0 && _antic_test_multiplier <= 1000))
                _antic_test_multiplier = 1;
        }
    }

    return _antic_test_multiplier;
}
