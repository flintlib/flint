/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>


// whether or not the compiler suppors __attribute__(__optimize__(..))
#ifndef HAVE_OPTIMIZE_ATTRIBUTE
#ifdef __GNUC__
#define HAVE_OPTIMIZE_ATTRIBUTE (__GNUC__ >= 4 && __GNUC_MINOR__ >= 4)
#else
#define HAVE_OPTIMIZE_ATTRIBUTE 0
#endif
#endif

#ifndef DO_CODEGEN_CHECKS
#if !HAVE_OPTIMIZE_ATTRIBUTE
#define DO_CODEGEN_CHECKS 0
#elif defined(__GNUC__)
// XXX this will need tweaking - look here if things break
//     b/c of compiler versions
#define DO_CODEGEN_CHECKS \
    (__GNUC__ == 4 && __GNUC_MINOR__ == 7 && __GNUC_PATCHLEVEL__ == 3 \
     && __x86_64__)
#else
#define DO_CODEGEN_CHECKS 0
#endif
#endif

#if !DO_CODEGEN_CHECKS
#define EXIT_STATEMENT throw skippable_exception("did not expect pass anyway")
#endif

#include <string>
#include <stdexcept>

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

#include "flintxx/test/helpers.h"

#include "fmpz_matxx.h"
#include "fmpz_polyxx.h"
#include "fmpzxx.h"
#include "nmod_matxx.h"
#include "nmod_polyxx.h"
#include "nmod_vecxx.h"

#include "flintxx/tuple.h"
#include "flintxx/vector.h"

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
//  movsbl %r9b,%r9d
//  mov    %r8d,(%rdi)
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
        {
            std::string tmp = ls[i].substr(pos+2);
            pos = tmp.find("\t");
            if(pos == std::string::npos)
                result += tmp;
            result += tmp.substr(pos+1);
        }
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

typedef make_vector_n<fmpzxx, 10>::type fmpzxx_vector;

extern "C" {
#if HAVE_OPTIMIZE_ATTRIBUTE
#define DEFINE_FUNC(name, args) \
void name args __attribute__((__optimize__("no-exceptions"))); \
void name args 
#else
#define DEFINE_FUNC(name, args) \
void name args 
#endif

DEFINE_FUNC(test_tuple_merge_1,
        (int& o1, int& o2, int& o3, int i1, slong i2, char i3, int i4))
{
    typedef make_tuple<slong, int, char> make3;
    typedef make_tuple<int, int> make2;
    typedef merge_tuple<make3::type, make2::type> merge;
    merge::type merged; // <char, int, slong, int>
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
        (int& o1, int& o2, int& o3, int i1, slong i2, char i3, int i4))
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
        (make_tuple<int, slong, char>::type& out, int i1, slong i2, char i3))
{
    typedef make_tuple<int*, slong*, char*>::type type;
    type pointers;
    back_tuple<type>::init(pointers, out, 0);
    *pointers.head = i1;
    *pointers.tail.head = i2;
    *pointers.tail.tail.head = i3;
}

DEFINE_FUNC(test_tuple_back_2,
        (make_tuple<int, slong, char>::type& out, int i1, slong i2, char i3))
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


DEFINE_FUNC(test_fmpzxx_symadd_1,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    out = ((a + b) + (c + d)) + ((a + c) + (b + d));
}

DEFINE_FUNC(test_fmpzxx_symadd_2,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    fmpz_t tmp1, tmp2, tmp3;
    fmpz_init(tmp1);
    fmpz_init(tmp2);
    fmpz_init(tmp3);

    fmpz_add(tmp1, a._fmpz(), b._fmpz());
    fmpz_add(tmp2, c._fmpz(), d._fmpz());
    fmpz_add(tmp3, tmp1, tmp2);
    fmpz_add(tmp1, a._fmpz(), c._fmpz());
    fmpz_add(tmp2, b._fmpz(), d._fmpz());
    fmpz_add(tmp1, tmp1, tmp2);
    fmpz_add(out._fmpz(), tmp1, tmp3);

    fmpz_clear(tmp1);
    fmpz_clear(tmp2);
    fmpz_clear(tmp3);
}

DEFINE_FUNC(test_fmpzxx_asymadd_1,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    out = (a + (((b + (c + (a + b))) + c) + d));
}

DEFINE_FUNC(test_fmpzxx_asymadd_2,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    fmpz_t tmp;
    fmpz_init(tmp);

    fmpz_add(tmp, a._fmpz(), b._fmpz());
    fmpz_add(tmp, c._fmpz(), tmp);
    fmpz_add(tmp, b._fmpz(), tmp);
    fmpz_add(tmp, tmp, c._fmpz());
    fmpz_add(tmp, tmp, d._fmpz());
    fmpz_add(out._fmpz(), a._fmpz(), tmp);

    fmpz_clear(tmp);
}

DEFINE_FUNC(test_fmpzxx_asymadd_3,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    out = a + b + c + a + b + c + d;
}

DEFINE_FUNC(test_fmpzxx_ternary_1,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    out = a + ((a+b) + c*d);
}

DEFINE_FUNC(test_fmpzxx_ternary_2,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    fmpz_t tmp;
    fmpz_init(tmp);

    fmpz_add(tmp, a._fmpz(), b._fmpz());
    fmpz_addmul(tmp, c._fmpz(), d._fmpz());
    fmpz_add(out._fmpz(), a._fmpz(), tmp);

    fmpz_clear(tmp);
}

DEFINE_FUNC(test_fmpzxx_ternary_3,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    out = a + ((a+b) + ((c+a) + (a+d))*d);
}

DEFINE_FUNC(test_fmpzxx_ternary_4,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    fmpz_t tmp1, tmp2;
    fmpz_init(tmp1);fmpz_init(tmp2);

    fmpz_add(tmp1, c._fmpz(), a._fmpz());
    fmpz_add(tmp2, a._fmpz(), d._fmpz());
    fmpz_add(tmp1, tmp1, tmp2);
    fmpz_add(tmp2, a._fmpz(), b._fmpz());
    fmpz_addmul(tmp2, tmp1, d._fmpz());
    fmpz_add(out._fmpz(), a._fmpz(), tmp2);

    fmpz_clear(tmp1);fmpz_clear(tmp2);
}

DEFINE_FUNC(test_fmpzxx_ternary_5,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    out += a*b;
}

DEFINE_FUNC(test_fmpzxx_ternary_6,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    fmpz_addmul(out._fmpz(), a._fmpz(), b._fmpz());
}

DEFINE_FUNC(test_fmpzxx_ternary_7,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    out += a*(b+c);
}

DEFINE_FUNC(test_fmpzxx_ternary_8,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    fmpz_t tmp; fmpz_init(tmp);
    fmpz_add(tmp, b._fmpz(), c._fmpz());
    fmpz_addmul(out._fmpz(), a._fmpz(), tmp);
    fmpz_clear(tmp);
}

DEFINE_FUNC(test_fmpzxx_ternary_9,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    out += (a+d)*(b+c);
}

DEFINE_FUNC(test_fmpzxx_ternary_10,
        (fmpzxx& out, const fmpzxx& a, const fmpzxx& b, const fmpzxx& c, const fmpzxx& d))
{
    fmpz_t tmp1, tmp2; fmpz_init(tmp1); fmpz_init(tmp2);
    fmpz_add(tmp1, b._fmpz(), c._fmpz());
    fmpz_add(tmp2, a._fmpz(), d._fmpz());
    fmpz_addmul(out._fmpz(), tmp1, tmp2);
    fmpz_clear(tmp1); fmpz_clear(tmp2);
}

DEFINE_FUNC(test_fmpzxx_ref_1,
        (fmpzxx_ref out, fmpzxx_srcref in1, fmpzxx_srcref in2))
{
    out = in1 + in2;
}

DEFINE_FUNC(test_fmpzxx_ref_2,
        (fmpzxx_ref out, fmpzxx_srcref in1, fmpzxx_srcref in2))
{
    fmpz_add(out._fmpz(), in1._fmpz(), in2._fmpz());
}

DEFINE_FUNC(test_fmpzxx_ltuple_comp_1,
        (bool& out, const fmpzxx& in1, const fmpzxx& in2))
{
    out = (ltupleref(in1, in2) == ltuple(1, 2));
}

DEFINE_FUNC(test_fmpzxx_ltuple_comp_2,
        (bool& out, const fmpzxx& in1, const fmpzxx& in2))
{
    out = ((fmpz_cmp_si(in1._fmpz(), 1) == 0)
           && (fmpz_cmp_si(in2._fmpz(), 2) == 0));
}

DEFINE_FUNC(test_fmpzxx_vec_1,
        (fmpzxx_vector& out, const fmpzxx_vector& in1, const fmpzxx_vector& in2))
{
    out = in1 + in2;
}

DEFINE_FUNC(test_fmpzxx_vec_2,
        (fmpz* out, const fmpz* in1, const fmpz* in2))
{
    fmpz_add(out + 0, in1 + 0, in2 + 0);
    fmpz_add(out + 1, in1 + 1, in2 + 1);
    fmpz_add(out + 2, in1 + 2, in2 + 2);
    fmpz_add(out + 3, in1 + 3, in2 + 3);
    fmpz_add(out + 4, in1 + 4, in2 + 4);
    fmpz_add(out + 5, in1 + 5, in2 + 5);
    fmpz_add(out + 6, in1 + 6, in2 + 6);
    fmpz_add(out + 7, in1 + 7, in2 + 7);
    fmpz_add(out + 8, in1 + 8, in2 + 8);
    fmpz_add(out + 9, in1 + 9, in2 + 9);
}


DEFINE_FUNC(test_fmpz_matxx_1,
        (fmpzxx& to, const fmpz_matxx& A))
{
    to = trace(transpose(A)*A);
}
DEFINE_FUNC(test_fmpz_matxx_2,
        (fmpzxx& to, const fmpz_matxx& A))
{
    fmpz_mat_t tmp1, tmp2;
    fmpz_mat_init(tmp1, fmpz_mat_nrows(A._mat()), fmpz_mat_ncols(A._mat()));
    fmpz_mat_init(tmp2, fmpz_mat_ncols(A._mat()), fmpz_mat_ncols(A._mat()));
    fmpz_mat_transpose(tmp1, A._mat());
    fmpz_mat_mul(tmp1, tmp2, A._mat());
    fmpz_mat_trace(to._fmpz(), tmp2);
    fmpz_mat_clear(tmp1);
    fmpz_mat_clear(tmp2);
}

DEFINE_FUNC(test_fmpz_matxx_manip_1,
        (fmpz_matxx& A))
{
    A.at(0, 0) += 2u;
}

DEFINE_FUNC(test_fmpz_matxx_manip_2,
        (fmpz_matxx& A))
{
    fmpz_add_ui(fmpz_mat_entry(A._mat(), 0, 0),
            fmpz_mat_entry(A._mat(), 0, 0), 2u);
}


DEFINE_FUNC(test_fmpz_polyxx_divrem_1,
        (fmpz_polyxx& A, fmpz_polyxx& B,
         const fmpz_polyxx& f, const fmpz_polyxx& g))
{
    ltupleref(A, B) = divrem(f, g);
}
DEFINE_FUNC(test_fmpz_polyxx_divrem_2,
        (fmpz_polyxx& A, fmpz_polyxx& B,
         const fmpz_polyxx& f, const fmpz_polyxx& g))
{
    fmpz_poly_divrem(A._poly(), B._poly(), f._poly(), g._poly());
}

DEFINE_FUNC(test_fmpz_polyxx_divrem_3,
        (fmpz_polyxx& A,
         const fmpz_polyxx& f, const fmpz_polyxx& g))
{
    ltupleref(A, _) = divrem(f, g);
}
DEFINE_FUNC(test_fmpz_polyxx_divrem_4,
        (fmpz_polyxx& A,
         const fmpz_polyxx& f, const fmpz_polyxx& g))
{
    fmpz_poly_t tmp;
    fmpz_poly_init(tmp);
    fmpz_poly_divrem(A._poly(), tmp, f._poly(), g._poly());
    fmpz_poly_clear(tmp);
}

DEFINE_FUNC(test_nmodxx_1,
        (mp_limb_t& to, mp_limb_t a_, mp_limb_t b_, nmod_t nm))
{
    nmodxx_ctx_srcref ctx = nmodxx_ctx_srcref::make(nm);
    nmodxx a = nmodxx::make_nored(a_, ctx);
    nmodxx b = nmodxx::make_nored(b_, ctx);
    to = inv((a + b) * (a * b)).to<mp_limb_t>();
}
DEFINE_FUNC(test_nmodxx_2,
        (mp_limb_t& to, mp_limb_t a, mp_limb_t b, nmod_t nm))
{
    mp_limb_t tmp1, tmp2;
    tmp1 = nmod_add(a, b, nm);
    tmp2 = nmod_mul(a, b, nm);
    tmp1 = nmod_mul(tmp1, tmp2, nm);
    tmp1 = nmod_inv(tmp1, nm);
    to = tmp1;
}


DEFINE_FUNC(test_nmod_polyxx_1,
        (nmod_polyxx& to, const nmod_polyxx& p1, const nmod_polyxx& p2))
{
    to = (p1*p1) + (p2*p2);
}

DEFINE_FUNC(test_nmod_polyxx_2,
        (nmod_polyxx& to, const nmod_polyxx& p1, const nmod_polyxx& p2))
{
    nmod_poly_t tmp1, tmp2;
    nmod_poly_init_preinv(tmp1, p1.modulus(), p1.estimate_ctx()._nmod().ninv);
    nmod_poly_init_preinv(tmp2, p1.modulus(), p1.estimate_ctx()._nmod().ninv);

    nmod_poly_mul(tmp1, p1._poly(), p1._poly());
    nmod_poly_mul(tmp2, p2._poly(), p2._poly());
    nmod_poly_add(to._poly(), tmp1, tmp2);

    nmod_poly_clear(tmp1);
    nmod_poly_clear(tmp2);
}

DEFINE_FUNC(test_nmod_polyxx_3,
        (nmod_polyxx& to, const nmod_polyxx& p1, const nmod_polyxx& p2))
{
    to = ((p1*p1) + (p2*p2)) + ((p1*p2) + (p2*p1));
}
DEFINE_FUNC(test_nmod_polyxx_4,
        (nmod_polyxx& to, const nmod_polyxx& p1, const nmod_polyxx& p2))
{
    nmod_poly_t tmp1, tmp2, tmp3;
    nmod_poly_init_preinv(tmp1, p1.modulus(), p1.estimate_ctx()._nmod().ninv);
    nmod_poly_init_preinv(tmp2, p1.modulus(), p1.estimate_ctx()._nmod().ninv);
    nmod_poly_init_preinv(tmp3, p1.modulus(), p1.estimate_ctx()._nmod().ninv);

    nmod_poly_mul(tmp1, p1._poly(), p1._poly());
    nmod_poly_mul(tmp2, p2._poly(), p2._poly());
    nmod_poly_add(tmp1, tmp1, tmp2);
    nmod_poly_mul(tmp2, p1._poly(), p2._poly());
    nmod_poly_mul(tmp3, p2._poly(), p1._poly());
    nmod_poly_add(tmp2, tmp2, tmp3);
    nmod_poly_add(to._poly(), tmp1, tmp2);

    nmod_poly_clear(tmp1);
    nmod_poly_clear(tmp2);
    nmod_poly_clear(tmp3);
}

DEFINE_FUNC(test_nmod_matxx_1,
        (nmodxx& to, const nmod_matxx& A))
{
    to = trace(transpose(A)*A);
}
DEFINE_FUNC(test_nmod_matxx_2,
        (nmodxx& to, const nmod_matxx& A))
{
    nmod_mat_t tmp1, tmp2;
    nmod_mat_init(tmp1, nmod_mat_nrows(A._mat()), nmod_mat_ncols(A._mat()),
            A.modulus());
    nmod_mat_init(tmp2, nmod_mat_ncols(A._mat()), nmod_mat_ncols(A._mat()),
            A.modulus());
    nmod_mat_transpose(tmp1, A._mat());
    nmod_mat_mul(tmp1, tmp2, A._mat());
    nmod_mat_clear(tmp1);
    to.set_nored(nmod_mat_trace(tmp2));
    nmod_mat_clear(tmp2);
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
test_fmpzxx()
{
    std::string ass1 = disass(program, "test_fmpzxx_symadd_1");
    std::string ass2 = disass(program, "test_fmpzxx_symadd_2");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(fuzzy_equals(count(stripnop(ass1), "\n"),
                count(stripnop(ass2), "\n"), 0.1));

    ass1 = disass(program, "test_fmpzxx_asymadd_1");
    ass2 = disass(program, "test_fmpzxx_asymadd_2");
    std::string ass3 = disass(program, "test_fmpzxx_asymadd_3");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(count(ass1, "mov") == count(ass2, "mov"));
    tassert(count(ass3, "call") == count(ass2, "call"));
    tassert(count(ass3, "mov") == count(ass2, "mov"));

    // This is probably not the best idea.
    // (Actually the same code is generated, up to jump targets, register names
    //  and addresses.)
    tassert(count(ass1, "\n") == count(ass2, "\n"));
    tassert(count(ass3, "\n") == count(ass2, "\n"));

    ass1 = disass(program, "test_fmpzxx_ternary_1");
    ass2 = disass(program, "test_fmpzxx_ternary_2");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(count(ass1, "mov") == count(ass2, "mov"));
    tassert(count(ass1, "\n") == count(ass2, "\n")); // XXX

    ass1 = disass(program, "test_fmpzxx_ternary_3");
    ass2 = disass(program, "test_fmpzxx_ternary_4");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(count(ass1, "je") == count(ass2, "je"));
    tassert(count(ass1, "jmp") == count(ass2, "jmp"));
    tassert(fuzzy_equals(count(ass1, "mov") + count(ass1, "lea"),
                         count(ass2, "mov") + count(ass2, "lea"), 0.1));
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.1));

    ass1 = disass(program, "test_fmpzxx_ternary_6");
    ass2 = disass(program, "test_fmpzxx_ternary_5");
    tassert(count(ass1, "\n") == count(ass2, "\n")); // XXX

    ass1 = disass(program, "test_fmpzxx_ternary_7");
    ass2 = disass(program, "test_fmpzxx_ternary_8");
    tassert(count(ass1, "call") == count(ass2, "call"));
    // NB: ass1 is actually shorter?!?
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.1));

    ass1 = disass(program, "test_fmpzxx_ternary_9");
    ass2 = disass(program, "test_fmpzxx_ternary_10");
    tassert(count(ass1, "call") == count(ass2, "call"));
    // NB: ass1 is actually shorter?!?
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.1));

    ass1 = disass(program, "test_fmpzxx_ref_1");
    ass2 = disass(program, "test_fmpzxx_ref_2");
    tassert(count(ass1, "\n") == count(ass2, "\n"));
    ass1 = disass(program, "test_fmpzxx_ltuple_comp_1");
    ass2 = disass(program, "test_fmpzxx_ltuple_comp_2");
    tassert(count(ass1, "\n") == count(ass2, "\n"));
}

void
test_vector()
{
    std::string ass1 = disass(program, "test_fmpzxx_vec_1");
    std::string ass2 = disass(program, "test_fmpzxx_vec_2");
    tassert(stripaddr(ass1) == stripaddr(ass2));
}

void
test_mat()
{
    std::string ass1 = disass(program, "test_fmpz_matxx_1");
    std::string ass2 = disass(program, "test_fmpz_matxx_2");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(fuzzy_equals(count(stripnop(ass1), "\n"),
                count(stripnop(ass2), "\n"), 0.1));

    ass1 = disass(program, "test_fmpz_matxx_manip_1");
    ass2 = disass(program, "test_fmpz_matxx_manip_2");
    tassert(count(ass1, "call") == count(ass2, "call")); // 1
    tassert(count(ass1, "\n") == count(ass2, "\n"));
}

void
test_poly()
{
    std::string ass1 = disass(program, "test_fmpz_polyxx_divrem_1");
    std::string ass2 = disass(program, "test_fmpz_polyxx_divrem_2");
    tassert(stripaddr(ass1) == stripaddr(ass2));

    ass1 = disass(program, "test_fmpz_polyxx_divrem_1");
    ass2 = disass(program, "test_fmpz_polyxx_divrem_2");
    tassert(count(ass1, "call") == count(ass2, "call"));
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.1));
}

void
test_nmod()
{
    std::string ass1 = disass(program, "test_nmodxx_1");
    std::string ass2 = disass(program, "test_nmodxx_2");
    tassert(stripaddr(ass1) == stripaddr(ass2));

    ass1 = disass(program, "test_nmod_polyxx_1");
    ass2 = disass(program, "test_nmod_polyxx_2");
    tassert(fuzzy_equals(count(ass1, "\n"), count(ass2, "\n"), 0.2));
    tassert(count(ass1, "call") == count(ass2, "call"));

    ass1 = disass(program, "test_nmod_polyxx_3");
    ass2 = disass(program, "test_nmod_polyxx_4");
    tassert(fuzzy_equals(count(stripnop(ass1), "\n"),
                count(stripnop(ass2), "\n"), 0.2));
    tassert(count(ass1, "call") == count(ass2, "call"));

    ass1 = disass(program, "test_nmod_matxx_1");
    ass2 = disass(program, "test_nmod_matxx_2");
    tassert(fuzzy_equals(count(stripnop(ass1), "\n"),
                count(stripnop(ass2), "\n"), 0.05));
    tassert(count(ass1, "call") == count(ass2, "call"));
}

int
main(int argc, char** argv)
{
    std::cout << "codegen....";
    program = argv[0];
    try
    {
        test_tuple();
        test_fmpzxx();
        test_vector();
        test_mat();
        test_poly();
        test_nmod();
    }
    catch(skippable_exception e)
    {
        std::cout << "SKIP (" << e.what() << ")\n";
        return 0;
    }
    std::cout << "PASS" << std::endl;
    return 0;
}

