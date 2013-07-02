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
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>

#include "cxx/test/helpers.h"
#include "cxx/tuple.h"

#include "cxx/prototype.h"

// Exception class to indicate that this test cannot proceed, e.g. because
// binutils is not installed or because we cannot interpret the disassembled
// code.
// This should not cause the test to fail, since it is not portable.
class skippable_exception
    : public std::runtime_error
{
public:
    skippable_exception(const std::string& n) : std::runtime_error(n) {}
};

// Run a command and recver the output.
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

// Disassembled a function contained in the binary "program".
std::string disass(const char* program, const std::string& function)
{
    std::string all = exec(std::string("objdump -d ") + program);
    std::string marker = "<" + function + ">:";
    std::size_t start = all.find(marker);
    if(start == std::string::npos)
        throw skippable_exception("cannot find function " + function);
    all = all.substr(start + marker.size());
    std::size_t len = all.find("\n\n");
    if(len == std::string::npos)
        throw skippable_exception("cannot find end of function " + function);
    return all.substr(0, len);
}

std::vector<std::string>&
split(const std::string &s, char delim, std::vector<std::string>& elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}
std::vector<std::string> lines(const std::string& input)
{
    std::vector<std::string> result;
    return split(input, '\n', result);
}

// Strip addresses from disassembly, to allow comparison. E.g.
//
//  401790:       45 0f be c9             movsbl %r9b,%r9d
//  401794:       44 89 07                mov    %r8d,(%rdi)
//
// -->
//
//  45 0f be c9             movsbl %r9b,%r9d
//  44 89 07                mov    %r8d,(%rdi)
// 
std::string stripaddr(const std::string& input)
{
    std::string result;
    std::vector<std::string> ls = lines(input);
    for(unsigned i = 0;i < ls.size();++i)
    {
        if(ls[i] == "")
            continue;
        std::size_t pos = ls[i].find(":\t");
        if(pos == std::string::npos)
            result += ls[i];
        else
            result += ls[i].substr(pos+2);
        result += '\n';
    }
    return result;
}

std::string stripnop(const std::string& input)
{
    std::string result;
    std::vector<std::string> ls = lines(input);
    for(unsigned i = 0;i < ls.size();++i)
    {
        if(ls[i].find("\tnop") != std::string::npos)
            continue;
        result += ls[i];
        result += '\n';
    }
    return result;
}

// Count the number of occurrences of a substring.
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

// Check if two quantities are equal up to a certain percentage of error.
template<class T>
bool fuzzy_equals(T v1, T v2, double perc)
{
    double d1 = double(v1);
    double d2 = double(v2);
    return d1*(1-perc) <= d2 && d2 <= d1*(1+perc);
}

using namespace flint;
using namespace mp;

extern "C" {
#define DEFINE_FUNC(name, args) \
void name args __attribute__((__optimize__("no-exceptions"))); \
void name args 

DEFINE_FUNC(test_tuple_merge_1,
        (int& o1, int& o2, int& o3, int i1, long i2, char i3, int i4))
{
    typedef make_tuple<long, int, char> make3;
    typedef make_tuple<int, int> make2;
    typedef merge_tuple<make3::type, make2::type> merge;
    merge::type merged; // <char, int, long, int>
    merged.head = i3;
    merged.tail.head = i1;
    merged.tail.tail.head = i2;
    merged.tail.tail.tail.head = i4;
    make3::type extract = merge::get_first(merged);
    o1 = extract.head; // i2
    o2 = extract.tail.head; // i1
    o3 = extract.tail.tail.head; // i3
}

DEFINE_FUNC(test_tuple_merge_2,
        (int& o1, int& o2, int& o3, int i1, long i2, char i3, int i4))
{
    o1 = i2;
    o2 = i1;
    o3 = i3;
}

DEFINE_FUNC(test_tuple_concat_1,
        (tuple<int, make_tuple<int, int, int>::type>& out,
         const make_tuple<int, int>::type& in1,
         const make_tuple<int, int>::type& in2))
{
    out = concat_tuple<
        make_tuple<int, int>::type,
        make_tuple<int, int>::type
      >::doit(in1, in2);
}

DEFINE_FUNC(test_tuple_concat_2,
        (tuple<int, make_tuple<int, int, int>::type>& out,
         const make_tuple<int, int>::type& in1,
         const make_tuple<int, int>::type& in2))
{
    out.head = in1.head;
    out.tail.head = in1.tail.head;
    out.tail.tail.head = in2.head;
    out.tail.tail.tail.head = in2.tail.head;
}

DEFINE_FUNC(test_tuple_back_1,
        (make_tuple<int, long, char>::type& out, int i1, long i2, char i3))
{
    typedef make_tuple<int*, long*, char*>::type type;
    type pointers;
    back_tuple<type>::init(pointers, out, 0);
    *pointers.head = i1;
    *pointers.tail.head = i2;
    *pointers.tail.tail.head = i3;
}

DEFINE_FUNC(test_tuple_back_2,
        (make_tuple<int, long, char>::type& out, int i1, long i2, char i3))
{
    out.head = i1;
    out.tail.head = i2;
    out.tail.tail.head = i3;
}

DEFINE_FUNC(test_tuple_extract_1,
        (make_homogeneous_tuple<int, 4>::type& out,
         const make_homogeneous_tuple<int, 7>::type& in))
{
    out = htuples::extract<4>(in);
}

DEFINE_FUNC(test_tuple_extract_2,
        (make_homogeneous_tuple<int, 4>::type& out,
         const make_homogeneous_tuple<int, 7>::type& in))
{
    out.head = in.head;
    out.tail.head = in.tail.head;
    out.tail.tail.head = in.tail.tail.head;
    out.tail.tail.tail.head = in.tail.tail.tail.head;
}

DEFINE_FUNC(test_tuple_removeres_1,
        (make_homogeneous_tuple<int, 2>::type& out, int in1, int in2))
{
    make_homogeneous_tuple<int, 3>::type tmp =
        make_tuple<int, int, int>::make(1, in1, in2);
    out = htuples::removeres(tmp, 1);
}

DEFINE_FUNC(test_tuple_removeres_2,
        (make_homogeneous_tuple<int, 2>::type& out, int in1, int in2))
{
    out.head = in1;
    out.tail.head = in2;
}


DEFINE_FUNC(test_mpz_symadd_1,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    out = ((a + b) + (c + d)) + ((a + c) + (b + d));
}

DEFINE_FUNC(test_mpz_symadd_2,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    fmpz_t tmp1, tmp2, tmp3;
    fmpz_init(tmp1);
    fmpz_init(tmp2);
    fmpz_init(tmp3);

    fmpz_add(tmp1, a._data(), b._data());
    fmpz_add(tmp2, c._data(), d._data());
    fmpz_add(tmp3, tmp1, tmp2);
    fmpz_add(tmp1, a._data(), c._data());
    fmpz_add(tmp2, b._data(), d._data());
    fmpz_add(tmp1, tmp1, tmp2);
    fmpz_add(out._data(), tmp1, tmp3);

    fmpz_clear(tmp1);
    fmpz_clear(tmp2);
    fmpz_clear(tmp3);
}

DEFINE_FUNC(test_mpz_asymadd_1,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    out = (a + (((b + (c + (a + b))) + c) + d));
}

DEFINE_FUNC(test_mpz_asymadd_2,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    fmpz_t tmp;
    fmpz_init(tmp);

    fmpz_add(tmp, a._data(), b._data());
    fmpz_add(tmp, c._data(), tmp);
    fmpz_add(tmp, b._data(), tmp);
    fmpz_add(tmp, tmp, c._data());
    fmpz_add(tmp, tmp, d._data());
    fmpz_add(out._data(), a._data(), tmp);

    fmpz_clear(tmp);
}

DEFINE_FUNC(test_mpz_asymadd_3,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    out = a + b + c + a + b + c + d;
}

DEFINE_FUNC(test_mpz_ternary_1,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    out = a + ((a+b) + c*d);
}

DEFINE_FUNC(test_mpz_ternary_2,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    fmpz_t tmp;
    fmpz_init(tmp);

    fmpz_add(tmp, a._data(), b._data());
    fmpz_addmul(tmp, c._data(), d._data());
    fmpz_add(out._data(), a._data(), tmp);

    fmpz_clear(tmp);
}

DEFINE_FUNC(test_mpz_ternary_3,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    out = a + ((a+b) + ((c+a) + (a+d))*d);
}

DEFINE_FUNC(test_mpz_ternary_4,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    fmpz_t tmp1, tmp2;
    fmpz_init(tmp1);fmpz_init(tmp2);

    fmpz_add(tmp1, c._data(), a._data());
    fmpz_add(tmp2, a._data(), d._data());
    fmpz_add(tmp1, tmp1, tmp2);
    fmpz_add(tmp2, a._data(), b._data());
    fmpz_addmul(tmp2, tmp1, d._data());
    fmpz_add(out._data(), a._data(), tmp2);

    fmpz_clear(tmp1);fmpz_clear(tmp2);
}

DEFINE_FUNC(test_mpz_ternary_5,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    out += a*b;
}

DEFINE_FUNC(test_mpz_ternary_6,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    fmpz_addmul(out._data(), a._data(), b._data());
}

DEFINE_FUNC(test_mpz_ternary_7,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    out += a*(b+c);
}

DEFINE_FUNC(test_mpz_ternary_8,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    fmpz_t tmp; fmpz_init(tmp);
    fmpz_add(tmp, b._data(), c._data());
    fmpz_addmul(out._data(), a._data(), tmp);
    fmpz_clear(tmp);
}

DEFINE_FUNC(test_mpz_ternary_9,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    out += (a+d)*(b+c);
}

DEFINE_FUNC(test_mpz_ternary_10,
        (mpz& out, const mpz& a, const mpz& b, const mpz& c, const mpz& d))
{
    fmpz_t tmp1, tmp2; fmpz_init(tmp1); fmpz_init(tmp2);
    fmpz_add(tmp1, b._data(), c._data());
    fmpz_add(tmp2, a._data(), d._data());
    fmpz_addmul(out._data(), tmp1, tmp2);
    fmpz_clear(tmp1); fmpz_clear(tmp2);
}
} // extern "C"

// Global variable, initialized by main.
const char* program = 0;

void
test_tuple()
{
    std::string ass1 = disass(program, "test_tuple_merge_1");
    std::string ass2 = disass(program, "test_tuple_merge_2");
    // XXX is this deterministic?
    tassert(stripaddr(ass1) == stripaddr(ass2));

    ass1 = disass(program, "test_tuple_concat_1");
    ass2 = disass(program, "test_tuple_concat_2");
    tassert(count(ass1, "\n") == count(ass2, "\n"));
    tassert(count(ass1, "mov") == count(ass2, "mov"));
    tassert(count(ass1, "call") == count(ass2, "call")); // 0

    ass1 = disass(program, "test_tuple_back_1");
    ass2 = disass(program, "test_tuple_back_2");
    // XXX is this deterministic?
    tassert(stripaddr(ass1) == stripaddr(ass2));

    ass1 = disass(program, "test_tuple_extract_1");
    ass2 = disass(program, "test_tuple_extract_2");
    tassert(count(ass1, "\n") == count(ass2, "\n"));
    tassert(count(ass1, "mov") == count(ass2, "mov"));
    tassert(count(ass1, "call") == count(ass2, "call")); // 0

    ass1 = disass(program, "test_tuple_removeres_1");
    ass2 = disass(program, "test_tuple_removeres_2");
    // XXX is this deterministic?
    tassert(stripaddr(ass1) == stripaddr(ass2));
}

void
test_mpz()
{
    std::string ass1 = disass(program, "test_mpz_symadd_1");
    std::string ass2 = disass(program, "test_mpz_symadd_2");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(fuzzy_equals(count(stripnop(ass1), "\n"),
                count(stripnop(ass2), "\n"), 0.1));

    ass1 = disass(program, "test_mpz_asymadd_1");
    ass2 = disass(program, "test_mpz_asymadd_2");
    std::string ass3 = disass(program, "test_mpz_asymadd_3");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(count(ass1, "mov") == count(ass2, "mov"));
    tassert(count(ass3, "call") == count(ass2, "call"));
    tassert(count(ass3, "mov") == count(ass2, "mov"));

    // This is probably not the best idea.
    // (Actually the same code is generated, up to jump targets, register names
    //  and addresses.)
    tassert(count(ass1, "\n") == count(ass2, "\n"));
    tassert(count(ass3, "\n") == count(ass2, "\n"));

    ass1 = disass(program, "test_mpz_ternary_1");
    ass2 = disass(program, "test_mpz_ternary_2");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(count(ass1, "mov") == count(ass2, "mov"));
    tassert(count(ass1, "\n") == count(ass2, "\n")); // XXX

    ass1 = disass(program, "test_mpz_ternary_3");
    ass2 = disass(program, "test_mpz_ternary_4");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(count(ass1, "je") == count(ass2, "je"));
    tassert(count(ass1, "jmp") == count(ass2, "jmp"));
    tassert(fuzzy_equals(count(ass1, "mov") + count(ass1, "lea"),
                         count(ass2, "mov") + count(ass2, "lea"), 0.1));
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.1));

    ass1 = disass(program, "test_mpz_ternary_6");
    ass2 = disass(program, "test_mpz_ternary_5");
    tassert(count(ass1, "\n") == count(ass2, "\n")); // XXX

    ass1 = disass(program, "test_mpz_ternary_7");
    ass2 = disass(program, "test_mpz_ternary_8");
    tassert(count(ass1, "call") == count(ass2, "call"));
    // NB: ass1 is actually shorter?!?
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.1));

    ass1 = disass(program, "test_mpz_ternary_9");
    ass2 = disass(program, "test_mpz_ternary_10");
    tassert(count(ass1, "call") == count(ass2, "call"));
    // NB: ass1 is actually shorter?!?
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.1));
}

int
main(int argc, char** argv)
{
    std::cout << "codegen....";
    program = argv[0];
    try
    {
        test_tuple();
        test_mpz();
    }
    catch(skippable_exception e)
    {
        std::cout << "SKIP (" << e.what() << ")\n";
        return 0;
    }
    std::cout << "PASS" << std::endl;
    return 0;
}

