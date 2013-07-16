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

#ifndef CXX_FMPZXX_H
#define CXX_FMPZXX_H

#include <cstdlib>

#include "cxx/evaluation_tools.h"
#include "cxx/expression.h"
#include "cxx/expression_traits.h"
#include "cxx/flint_classes.h"
#include "cxx/stdmath.h"

#include "fmpz.h"

namespace flint {

template<class Operation, class Data>
class fmpzxx_expression
    : public expression<derived_wrapper<fmpzxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpzxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS(fmpzxx_expression)
    FLINTXX_DEFINE_CTORS(fmpzxx_expression)
    FLINTXX_DEFINE_C_REF(fmpzxx_expression, fmpz, _fmpz)
};

namespace detail {
struct fmpz_data;
}

typedef fmpzxx_expression<operations::immediate, detail::fmpz_data> fmpzxx;
typedef fmpzxx_expression<operations::immediate,
            flint_classes::ref_data<fmpzxx, fmpz> > fmpzxx_ref;
typedef fmpzxx_expression<operations::immediate,
            flint_classes::srcref_data<fmpzxx, fmpzxx_ref, fmpz> > fmpzxx_srcref;

namespace detail {
struct fmpz_data
{
    typedef fmpz_t& data_ref_t;
    typedef const fmpz_t& data_srcref_t;

    fmpz_t inner;

    fmpz_data() {fmpz_init(inner);}
    ~fmpz_data() {fmpz_clear(inner);}
    fmpz_data(const fmpz_data& o) {fmpz_init_set(inner, o.inner);}

    fmpz_data(const char* str)
    {
        fmpz_init(inner);
        fmpz_set_str(inner, str, 10);
    }

    template<class T>
    fmpz_data(const T& t)
    {
        init(t);
    }

    template<class T>
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type init(T t)
    {
        fmpz_init_set_ui(inner, t);
    }
    template<class T>
    typename mp::enable_if<traits::is_signed_integer<T> >::type init(T t)
    {
        fmpz_init(inner);
        fmpz_set_si(inner, t);
    }

    void init(const fmpzxx_srcref& r)
    {
        fmpz_init_set(inner, r._fmpz());
    }
};

} // detail

///////////////////////////////////////////////////////////////////
// HELPERS
///////////////////////////////////////////////////////////////////
namespace traits {
template<class T> struct is_fmpzxx : is_T_expr<T, fmpzxx> { };
} // traits
namespace mp {
template<class Out, class T1, class T2 = void>
struct enable_all_fmpzxx
    : mp::enable_if<mp::and_<traits::is_fmpzxx<T1>, traits::is_fmpzxx<T2> >, Out> { };
template<class Out, class T>
struct enable_all_fmpzxx<Out, T, void>
    : mp::enable_if<traits::is_fmpzxx<T>, Out> { };
} // mp

///////////////////////////////////////////////////////////////////
// RULES
///////////////////////////////////////////////////////////////////
namespace rules {

#define FMPZXX_COND_S FLINTXX_COND_S(fmpzxx)
#define FMPZXX_COND_T FLINTXX_COND_T(fmpzxx)

FLINT_DEFINE_DOIT_COND2(assignment, FMPZXX_COND_T, FMPZXX_COND_S,
        fmpz_set(to._fmpz(), from._fmpz()))

FLINT_DEFINE_DOIT_COND2(assignment,
        FMPZXX_COND_T, traits::is_unsigned_integer,
        fmpz_set_ui(to._fmpz(), from))

FLINT_DEFINE_DOIT_COND2(assignment,
        FMPZXX_COND_T, traits::is_signed_integer,
        fmpz_set_si(to._fmpz(), from))

template<class T, int n>
struct assignment<T, char[n],
    typename mp::enable_if<FMPZXX_COND_S<T> >::type>
{
    static void doit(T& target, const char* source)
    {
        fmpz_set_str(target._fmpz(), const_cast<char*>(source), 10);
    }
};

FLINTXX_DEFINE_CMP(fmpzxx, fmpz_cmp(e1._fmpz(), e2._fmpz()))

template<class T, class U>
struct cmp<T, U,
    typename mp::enable_if<mp::and_<
        FMPZXX_COND_S<T>, traits::is_signed_integer<U> > >::type>
{
    static int get(const T& v, const U& t)
    {
        return fmpz_cmp_si(v._fmpz(), t);
    }
};

template<class T>
struct cmp<fmpzxx, T,
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static int get(const fmpzxx& v, const T& t)
    {
        return fmpz_cmp_ui(v._fmpz(), t);
    }
};

FLINTXX_DEFINE_TO_STR(fmpzxx, fmpz_get_str(0,  base, from._fmpz()))

FLINT_DEFINE_GET_COND(conversion, slong, FMPZXX_COND_S,
        fmpz_get_si(from._fmpz()))
FLINT_DEFINE_GET_COND(conversion, ulong, FMPZXX_COND_S,
        fmpz_get_ui(from._fmpz()))
FLINT_DEFINE_GET_COND(conversion, double, FMPZXX_COND_S,
        fmpz_get_d(from._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_add(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_CBINARY_EXPR_COND2(plus, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_add_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_mul(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_mul_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpzxx,
        FMPZXX_COND_S, traits::is_signed_integer,
        fmpz_mul_si(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_sub(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_sub_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_fdiv_q(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_fdiv_q_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpzxx,
        FMPZXX_COND_S, traits::is_signed_integer,
        fmpz_fdiv_q_si(to._fmpz(), e1._fmpz(), e2))

// TODO this interpretation of mod is not the same as for builtin types!
FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_mod(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_mod_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpzxx, FMPZXX_COND_S,
        fmpz_neg(to._fmpz(), from._fmpz()))
} // rules

FLINTXX_DEFINE_TERNARY(fmpzxx,
        fmpz_addmul(to._fmpz(), e1._fmpz(), e2._fmpz()),
        fmpz_submul(to._fmpz(), e1._fmpz(), e2._fmpz()),
        FLINTXX_UNADORNED_MAKETYPES)

///////////////////////////////////////////////////////////////////////////
// FUNCTIONS
///////////////////////////////////////////////////////////////////////////

// These functions evaluate immediately and do not yield fmpzxxs

template<class T1, class T2>
inline typename mp::enable_all_fmpzxx<bool, T1, T2>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible(t1.evaluate()._fmpz(), t2.evaluate()._fmpz());
}
template<class T1, class T2>
inline typename mp::enable_if<mp::and_<
    traits::is_fmpzxx<T1>, traits::fits_into_slong<T2> >, bool>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible_si(t1.evaluate()._fmpz(), t2);
}

// These functions are evaluated lazily

FLINT_DEFINE_UNOP(fac)
FLINT_DEFINE_BINOP(rfac)
FLINT_DEFINE_BINOP(bin)
namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(rfac_op, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_rfac_ui(to._fmpz(), e1._fmpz(), e2))
FLINT_DEFINE_UNARY_EXPR_COND(fac_op, fmpzxx, traits::is_unsigned_integer,
        fmpz_fac_ui(to._fmpz(), from))

template<class T1, class T2>
struct binary_expression<
    T1,
    typename mp::enable_if<
        mp::and_<
            traits::is_unsigned_integer<T1>,
            traits::is_unsigned_integer<T2> >,
        operations::bin_op>::type,
    T2>
{
    typedef fmpzxx return_t;
    template<class V>
    static void doit(V& to, const T1& t1, const T2& t2)
    {
        fmpz_bin_uiui(to._fmpz(), t1, t2);
    }
};

// standard math functions (c/f stdmath.h)
FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_pow_ui(to._fmpz(), e1._fmpz(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(root_op, fmpzxx,
        FMPZXX_COND_S, traits::fits_into_slong,
        fmpz_root(to._fmpz(), e1._fmpz(), e2))
FLINT_DEFINE_UNARY_EXPR_COND(sqrt_op, fmpzxx, FMPZXX_COND_S,
        fmpz_sqrt(to._fmpz(), from._fmpz()))
} // rules

// TODO many more functions

} // flint

#endif
