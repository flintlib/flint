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

#include "mp.h"
#include "traits.h"

namespace flint {
namespace rules {
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

// Rule for c-style printing
template<class T, class Enable = void>
struct cprint : UNIMPLEMENTED { };
// static int doit(FILE*, const T&)

template<class T, class Enable = void>
struct print_pretty : UNIMPLEMENTED { };
// static int doit(FILE*, const T&)

template<class T, class Enable = void>
struct read : UNIMPLEMENTED { };
// static int doit(FILE*, T&)

// Rule for swapping
template<class To, class From, class Enable = void>
struct swap : UNIMPLEMENTED { };

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
template<class T, class Op, class U, class Enable = void>
struct binary_expression : UNIMPLEMENTED { };
// typedef X return_t;
// static void doit(return_t& to, const T&, const U&);
template<class T, class Op, class U, class Enable = void>
struct commutative_binary_expression : UNIMPLEMENTED { };
// similarly

template<class Op, class T, class Enable = void>
struct unary_expression : UNIMPLEMENTED { };
// similarly

// Rules for more arguments.
template<class Op, class T, class U, class V, class Enable = void>
struct threeary_expression : UNIMPLEMENTED { };
template<class Op, class T, class U, class V, class W, class Enable = void>
struct fourary_expression : UNIMPLEMENTED { };
template<class Op, class T, class U, class V, class W, class X, class Enable = void>
struct fiveary_expression : UNIMPLEMENTED { };
template<class Op, class T, class U, class V, class W, class X, class Y, class Enable = void>
struct sixary_expression : UNIMPLEMENTED { };
template<class Op, class T, class U, class V, class W, class X, class Y, class Z, class Enable = void>
struct sevenary_expression : UNIMPLEMENTED { };
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
// In general the easiest way to find out what they do is to read the
// definition directly.

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
struct name<totype, T, typename mp::enable_if< cond<T> >::type> \
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
struct name<totype, T, typename mp::enable_if< cond<T> >::type> \
{ \
    static void doit(totype& to, const T& from) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_DOIT_COND2(name, cond1, cond2, eval) \
template<class T, class U> \
struct name<T, U, typename mp::enable_if<mp::and_< cond1 <T>, cond2 <U> > >::type> \
{ \
    static void doit(T& to, const U& from) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_PRINT_COND_(name, cond, eval) \
template<class T> \
struct name<T, typename mp::enable_if< cond <T> >::type> \
{ \
    static int doit(FILE* to, const T& from) \
    { \
        return eval; \
    } \
};
#define FLINT_DEFINE_READ_COND_(name, cond, eval) \
template<class T> \
struct name<T, typename mp::enable_if< cond <T> >::type> \
{ \
    static int doit(FILE* from, T& to) \
    { \
        return eval; \
    } \
};
#define FLINT_DEFINE_PRINT_COND(cond, eval) \
    FLINT_DEFINE_PRINT_COND_(cprint, cond, eval)
#define FLINT_DEFINE_PRINT_PRETTY_COND(cond, eval) \
    FLINT_DEFINE_PRINT_COND_(print_pretty, cond, eval)
#define FLINT_DEFINE_READ_COND(cond, eval) \
    FLINT_DEFINE_READ_COND_(read, cond, eval)
#define FLINT_DEFINE_PRINT_PRETTY_COND_2(cond, extratype, eval) \
template<class T> \
struct print_pretty<T, typename mp::enable_if< cond <T> >::type> \
{ \
    static int doit(FILE* to, const T& from, extratype extra) \
    { \
        return eval; \
    } \
};

// Specialise the unary expression rule type->type.
#define FLINT_DEFINE_UNARY_EXPR_(name, rtype, type, eval) \
template<> \
struct unary_expression<operations::name, type> \
{ \
    typedef rtype return_t; \
    template<class V> \
    static void doit(V& to, const type& from) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_UNARY_EXPR(name, type, eval) \
    FLINT_DEFINE_UNARY_EXPR_(name, type, type, eval)

#define FLINT_DEFINE_UNARY_EXPR_COND(name, ret_type, cond, eval) \
template<class T> \
struct unary_expression<operations::name, T, \
    typename mp::enable_if<cond<T> >::type> \
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
struct commutative_binary_expression<Type, operations::name, T, \
    typename mp::enable_if<cond<T> >::type> \
{ \
    typedef Type return_t; \
    template<class V> \
    static void doit(V& to, const Type& e1, const T& e2) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_CBINARY_EXPR_COND2(name, rettype, cond1, cond2, eval) \
template<class T, class U> \
struct commutative_binary_expression<T, operations::name, U, \
    typename mp::enable_if<mp::and_< cond1 <T>, cond2 <U> > >::type> \
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
struct binary_expression<Type, operations::name, T, \
    typename mp::enable_if<cond<T> >::type> \
{ \
    typedef Type return_t; \
    template<class V> \
    static void doit(V& to, const Type& e1, const T& e2) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_BINARY_EXPR_COND2(name, rettype, cond1, cond2, eval) \
template<class T, class U> \
struct binary_expression<T, operations::name, U, \
    typename mp::enable_if< mp::and_< cond1 <T>, cond2 <U> > >::type> \
{ \
    typedef rettype return_t; \
    template<class V> \
    static void doit(V& to, const T& e1, const U& e2) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_THREEARY_EXPR_COND3(name, rettype, cond1, cond2, cond3, eval) \
template<class T, class U, class V> \
struct threeary_expression<operations::name, T, U, V, \
    typename mp::enable_if< mp::and_< cond1 <T>, cond2 <U>, cond3 <V> > >::type> \
{ \
    typedef rettype return_t; \
    template<class R> \
    static void doit(R& to, const T& e1, const U& e2, const V& e3) \
    { \
        eval; \
    } \
};
#define FLINT_DEFINE_THREEARY_EXPR(name, rettype, T1, T2, T3, eval) \
template<> \
struct threeary_expression<operations::name, T1, T2, T3> \
{ \
    typedef rettype return_t; \
    template<class R> \
    static void doit(R& to, const T1& e1, const T2& e2, const T3& e3) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_FOURARY_EXPR_COND4(name, rettype, cond1, cond2, cond3, cond4, eval) \
template<class T, class U, class V, class W> \
struct fourary_expression<operations::name, T, U, V, W, \
    typename mp::enable_if< mp::and_< cond1 <T>, cond2 <U>, cond3 <V>, cond4 <W> > >::type> \
{ \
    typedef rettype return_t; \
    template<class R> \
    static void doit(R& to, const T& e1, const U& e2, const V& e3, const W& e4) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_FIVEARY_EXPR_COND5(name, rettype, cond1, cond2, cond3, cond4, cond5, eval) \
template<class T, class U, class V, class W, class X> \
struct fiveary_expression<operations::name, T, U, V, W, X, \
    typename mp::enable_if< mp::and_< cond1 <T>, cond2 <U>, cond3 <V>, cond4 <W>, cond5 <X> > >::type> \
{ \
    typedef rettype return_t; \
    template<class R> \
    static void doit(R& to, const T& e1, const U& e2, const V& e3, const W& e4, const X& e5) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_SIXARY_EXPR_COND6(name, rettype, cond1, cond2, cond3, cond4, cond5, cond6, eval) \
template<class T, class U, class V, class W, class X, class Y> \
struct sixary_expression<operations::name, T, U, V, W, X, Y, \
    typename mp::enable_if< mp::and_< cond1 <T>, cond2 <U>, cond3 <V>, cond4 <W>, cond5 <X>, cond6<Y> > >::type> \
{ \
    typedef rettype return_t; \
    template<class R> \
    static void doit(R& to, const T& e1, const U& e2, const V& e3, const W& e4, const X& e5, const Y& e6) \
    { \
        eval; \
    } \
};

#define FLINT_DEFINE_SEVENARY_EXPR_COND7(name, rettype, cond1, cond2, cond3, cond4, cond5, cond6, cond7, eval) \
template<class T, class U, class V, class W, class X, class Y, class Z> \
struct sevenary_expression<operations::name, T, U, V, W, X, Y, Z, \
    typename mp::enable_if< mp::and_< cond1 <T>, cond2 <U>, cond3 <V>, cond4 <W>, cond5 <X>, cond6<Y>, cond7<Z> > >::type> \
{ \
    typedef rettype return_t; \
    template<class R> \
    static void doit(R& to, const T& e1, const U& e2, const V& e3, const W& e4, const X& e5, const Y& e6, const Z& e7) \
    { \
        eval; \
    } \
};

#endif
