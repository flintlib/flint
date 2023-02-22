/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program to demonstrate some use of the padic module.
*/

#include <iostream>
#include "padicxx.h"

using namespace flint;
using namespace std;

int main()
{
    std::cout << "Output:\n\n";

    // Case 1
    {
        std::cout << "Positive integer:  x = 127 mod 7^10\n";
        fmpzxx p(7);
        padicxx_ctx ctx(p, 8, 12, PADIC_TERSE);
        padicxx x = padicxx::from_QQ(127, ctx, 10);

        ctx.mode() = PADIC_TERSE;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';

        ctx.mode() = PADIC_SERIES;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';

        ctx.mode() = PADIC_VAL_UNIT;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';
    }

    // Case 2
    {
        std::cout << "Positive integer larger than p^N:  x = 1057 mod 2^10\n";
        fmpzxx p(2);
        padicxx_ctx ctx(p, 10, 12, PADIC_TERSE);
        padicxx x = padicxx::from_QQ(1057, ctx, 10);

        ctx.mode() = PADIC_TERSE;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';

        ctx.mode() = PADIC_SERIES;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';

        ctx.mode() = PADIC_VAL_UNIT;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';
    }

    // Case 3
    {
        std::cout << "Negative integer:  x = -127 mod 3^10\n";
        fmpzxx p(3);
        padicxx_ctx ctx(p, 10, 12, PADIC_TERSE);
        padicxx x = padicxx::from_QQ(-127, ctx, 10);

        ctx.mode() = PADIC_TERSE;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';

        ctx.mode() = PADIC_SERIES;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';

        ctx.mode() = PADIC_VAL_UNIT;
        std::cout << "print:   ";print(x);std::cout << '\n';
        std::cout << "get_str: " << x.to_string() << '\n';
    }

    // Log
    {
        std::cout << "Log of 7380996 mod 5^20\n";
        fmpzxx p(5);
        padicxx_ctx ctx(p, 10, 25, PADIC_SERIES);

        padicxx x = padicxx::from_QQ(7380996, ctx);
        padicxx y(log(x));

        std::cout << "x = " << x << '\n';
        std::cout << "y = " << y << '\n';
    }

    return 0;
}

