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

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/
#include <iostream>
#include <string>
#include <stdexcept>
#include <stdio.h>

#include "cxx.h"
#include "cxx/test/helpers.h"

class skippable_exception
    : public std::runtime_error
{
public:
    skippable_exception(const std::string& n) : std::runtime_error(n) {}
};

std::string exec(const std::string& cmd)
{
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe)
        throw skippable_exception("cannot execute command" + cmd);
    char buffer[128];
    std::string result = "";
    while(!feof(pipe))
    {
        if(fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;
}

int count(const std::string& str, const std::string& sub)
{
    if(sub.length() == 0)
        return 0;
    int count = 0;
    for(size_t offset = str.find(sub); offset != std::string::npos;
            offset = str.find(sub, offset + sub.length()))
        ++count;
    return count;
}

template<class T>
bool fuzzy_equals(T v1, T v2, double perc)
{
    double d1 = double(v1);
    double d2 = double(v2);
    return d1*(1-perc) <= d2 && d2 <= d1*(1+perc);
}

std::string disass(const char* program, const std::string& function)
{
    std::string all = exec(std::string("objdump -d ") + program);
    size_t start = all.find("<" + function + ">:");
    if(start == std::string::npos)
        throw skippable_exception("cannot find function " + function);
    all = all.substr(start);
    size_t len = all.find("\n\n");
    if(len == std::string::npos)
        throw skippable_exception("cannot find end of function " + function);
    return all.substr(0, len);
}


extern "C" {
void
test1(flint::mpz& out, const flint::mpz& a, const flint::mpz& b,
     const flint::mpz& c, const flint::mpz& d) __attribute__((__optimize__("no-exceptions")));
void
test1(flint::mpz& out, const flint::mpz& a, const flint::mpz& b,
     const flint::mpz& c, const flint::mpz& d)
{
    out = ((a + b) + (c + d)) + ((a + c) + (b + d));
}

void
test2(flint::mpz& out, const flint::mpz& a, const flint::mpz& b,
     const flint::mpz& c, const flint::mpz& d) __attribute__((__optimize__("no-exceptions")));
void
test2(flint::mpz& out, const flint::mpz& a, const flint::mpz& b,
     const flint::mpz& c, const flint::mpz& d)
{
    fmpz_t tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
    fmpz_init(tmp1);
    fmpz_init(tmp2);
    fmpz_init(tmp3);
    fmpz_init(tmp4);
    fmpz_init(tmp5);
    fmpz_init(tmp6);

    fmpz_add(tmp1, a._data(), b._data());
    fmpz_add(tmp2, c._data(), d._data());
    fmpz_add(tmp3, tmp1, tmp2);
    fmpz_add(tmp4, a._data(), c._data());
    fmpz_add(tmp5, b._data(), d._data());
    fmpz_add(tmp6, tmp1, tmp2);
    fmpz_add(out._data(), tmp3, tmp6);

    fmpz_clear(tmp1);
    fmpz_clear(tmp2);
    fmpz_clear(tmp3);
    fmpz_clear(tmp4);
    fmpz_clear(tmp5);
    fmpz_clear(tmp6);
}
}

// TODO move all codegen stuff to t-codegen
int
main(int argc, char** argv)
{
    std::cout << "mpz....";
    //flint::mpz a('a');
    flint::mpz a((short)-4);
    flint::mpz b; b = "7";
    //std::cout << a + b << std::endl;
    std::string ass1 = disass(argv[0], "test1");
    std::string ass2 = disass(argv[0], "test2");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.1));
    std::cout << "PASS" << std::endl;
}
