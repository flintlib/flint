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

// This file contains the definitions of all rules used by the expression class.
// (Some generally useful implementations can be found in default_rules.h.)
// This file also contains some helper traits, metaprogramming tools and macros.

#ifndef CXX_RULES_H
#define CXX_RULES_H

#include "cxx/mp.h"
#include "cxx/traits.h"

namespace flint {
namespace rules {
// TODO document these
struct no_op
{
    template<class U>
    static void doit(const U&) {}
};

struct UNIMPLEMENTED
{
    static const bool unimplemented_marker = true;
};

template<class T, class Enable = void>
struct print : UNIMPLEMENTED { };

template<class T, class Enable = void>
struct to_string : UNIMPLEMENTED { };
// static std::string get(const T&, int base)

template<class T, class U, class Enable = void>
struct assignment : UNIMPLEMENTED { };

// C-style cmp.
template<class T, class U, class Enable = void>
struct cmp : UNIMPLEMENTED { };

// Rule for equals. Implemented in terms of cmp by default.
template<class T, class U, class Enable = void>
struct equals : UNIMPLEMENTED { };

// Rule for type conversion.
template<class To, class From, class Enable = void>
struct conversion
{
    static To get(const From& from)
    {
        return To(from);
    }
};

// If result_is_temporary is true, then the result coincides with the
// first temporary (provided these have the same type)
// Priorities 2, 1, 0 can be used to resolve conflicts.
template<
    class Op, class Data,
    bool result_is_temporary,
    unsigned priority,
    class Enable = void>
struct evaluation : UNIMPLEMENTED { };
//{
//    typedef X return_t;
//    typedef Y temporaries_t; // a tuple of *pointers*
//    static void doit(const T& input, temporaries_t temps, return_t* output);
//};

// Instantiate temporaries for evaluation. The default implementation does the
// following:
// - find a subexpression t of e which evaluates to type T
// -- return t.create_temporary()
// - if no such subexpression can be found, return T()
// Additionally, the expression class implements a version of create_temporary
// which just returns T(), so if your class is default constructible,
// everything works automatically.
template<class Expr, class T, class Enable = void>
struct instantiate_temporaries;
//{
//    static T get(const Expr& e);
//};

// Convenience helpers, instantiate by evaluation if necessary
// (needs default rules)
template<class T, class Op, class U>
struct binary_expression : UNIMPLEMENTED { };
// typedef X return_t;
// static void doit(return_t& to, const T&, const U&);
template<class T, class Op, class U>
struct commutative_binary_expression : UNIMPLEMENTED { };
// similarly

template<class Op, class T>
struct unary_expression : UNIMPLEMENTED { };
// similarly
} // rules


///////////////////////////////////////////////////////////////////////////////////
// HELPER TRAITS
///////////////////////////////////////////////////////////////////////////////////

namespace traits {
// Compute if the rule T is implemented.
template<class T>
struct is_implemented : mp::not_<_is_convertible<rules::UNIMPLEMENTED, T> > { };
} // traits
} // flint


///////////////////////////////////////////////////////////////////////////////////
// HELPER MACROS
///////////////////////////////////////////////////////////////////////////////////

// These macros should be called in namespace flint::rules

// Specialise a getter called "name". The getter has one argument called "from"
// of type "fromtype", and "eval" which should yield "totype". 
#define FLINT_DEFINE_GET2(name, totype, fromtype1, fromtype2, eval) \
template<> \
struct name<fromtype1, fromtype2> \
{ \
    static totype get(const fromtype1& e1, const fromtype2& e2) \
    { \
        return eval; \
    } \
};

#define FLINT_DEFINE_GET(name, totype, fromtype, eval) \
    FLINT_DEFINE_GET2(name, totype, fromtype, fromtype, eval)

#define FLINT_DEFINE_GET_COND(name, totype, cond, eval) \
template<class T> \
struct name<totype, T, typename mp::enable_if< cond >::type> \
{ \
    static totype get(const T& from) \
    { \
        return eval; \
    } \
};

// Specialise a doit rule called "name"
#define FLINT_DEFINE_DOIT(name, totype, fromtype, eval) \
template<> \
struct name<totype, fromtype> \
{ \
    static void doit(totype& to, const fromtype& from) \
    { \
        eval; \
    } \
};

// Specialise a doit rule called "name" which yields totype. It will
// accept any type "T" which satisfies "cond".
#define FLINT_DEFINE_DOIT_COND(name, totype, cond, eval) \
template<class T> \
struct name<totype, T, typename mp::enable_if< cond >::type> \
{ \
    static void doit(totype& to, const T& from) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_DOIT_COND2(name, cond, eval) \
template<class T, class U> \
struct name<T, U, typename mp::enable_if< cond >::type> \
{ \
    static void doit(T& to, const U& from) \
    { \
        eval; \
    } \
};

// Specialise the unary expression rule type->type.
#define FLINT_DEFINE_UNARY_EXPR(name, type, eval) \
template<> \
struct unary_expression<operations::name, type> \
{ \
    typedef type return_t; \
    template<class V> \
    static void doit(V& to, const type& from) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_UNARY_EXPR_COND(name, ret_type, cond, eval) \
template<class T> \
struct unary_expression<typename mp::enable_if<cond, operations::name>::type, T> \
{ \
    typedef ret_type return_t; \
    template<class V> \
    static void doit(V& to, const T& from) \
    { \
        eval; \
    } \
};

// Specialise the binary expression rule (type, type) -> type
#define FLINT_DEFINE_BINARY_EXPR2(name, rtype, type1, type2, eval) \
template<> \
struct binary_expression<type1, operations::name, type2> \
{ \
    typedef rtype return_t; \
    template<class V> \
    static void doit(V& to, const type1& e1, const type2& e2) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_BINARY_EXPR(name, type, eval) \
    FLINT_DEFINE_BINARY_EXPR2(name, type, type, type, eval)

// Specialise the commutative binary expression rule (type, type) -> type
#define FLINT_DEFINE_CBINARY_EXPR(name, type, eval) \
template<> \
struct commutative_binary_expression<type, operations::name, type> \
{ \
    typedef type return_t; \
    template<class V> \
    static void doit(V& to, const type& e1, const type& e2) \
    { \
        eval; \
    } \
};

// Specialise the commutative binary expression rule (Type, T) -> Type,
// where T must satisfy "cond".
#define FLINT_DEFINE_CBINARY_EXPR_COND(name, Type, cond, eval) \
template<class T> \
struct commutative_binary_expression<Type, \
    typename mp::enable_if<cond, operations::name>::type, T> \
{ \
    typedef Type return_t; \
    template<class V> \
    static void doit(V& to, const Type& e1, const T& e2) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_CBINARY_EXPR_COND2(name, rettype, cond, eval) \
template<class T, class U> \
struct commutative_binary_expression<T, \
    typename mp::enable_if< cond, operations::name>::type, U> \
{ \
    typedef rettype return_t; \
    template<class V> \
    static void doit(V& to, const T& e1, const U& e2) \
    { \
        eval; \
    } \
};

// Specialise the (non-commutative) binary expression rule (Type, T) -> Type,
// where T must satisfy "cond".
#define FLINT_DEFINE_BINARY_EXPR_COND(name, Type, cond, eval) \
template<class T> \
struct binary_expression<Type, \
    typename mp::enable_if<cond, operations::name>::type, T> \
{ \
    typedef Type return_t; \
    template<class V> \
    static void doit(V& to, const Type& e1, const T& e2) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_BINARY_EXPR_COND2(name, rettype, cond, eval) \
template<class T, class U> \
struct binary_expression<T, \
    typename mp::enable_if< cond, operations::name>::type, U> \
{ \
    typedef rettype return_t; \
    template<class V> \
    static void doit(V& to, const T& e1, const U& e2) \
    { \
        eval; \
    } \
};

#endif
