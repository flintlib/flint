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

// Helpers to define concrete subclasses of expression.
// Contrary to other parts of this library, they are tailored very
// specifically towards FLINT.

#include "../flint.h"
#include "mp.h"
#include "expression.h"
#include "expression_traits.h"
#include "evaluation_tools.h"
#include "tuple.h"

#ifndef CXX_FLINT_CLASSES_H
#define CXX_FLINT_CLASSES_H

// Flint classes distinguish themselves from "ordinary" expression template
// classes by a public typedef IS_FLINT_CLASS (see also FLINTXX_DEFINE_BASICS
// below). Most functionality in this header disables itself when used on a
// non-flint class.

// For all flint classes, Data of immediates must have typedefs data_ref_t and
// data_srcref_t.

// The immediates of any flint class come in three "flavours": ordinary, ref
// and srcref. Most of the classes below are used to convert these flavours.
// In order for this to work, the expression template class must contain a
// public typedef c_base_t which is the underlying "basic" C type (so not the
// length one array that is usually used), e.g. fmpz_poly_struct. Conversion to
// reference type then yields expression templates which have Data set to
// ref_data<Wrapped, c_base_t> or srcref_data. These implementation work as
// long as all data is stored in c_base_t. If not (e.g. for padic, where there
// is an additional reference to the context), ref_data and srcref_data have to
// be specialised appropriately.

namespace flint {
namespace flint_classes {
template<class Wrapped, class Inner>
struct ref_data
{
    typedef void IS_REF_OR_CREF;
    typedef Wrapped wrapped_t;

    typedef Inner* data_ref_t;
    typedef const Inner* data_srcref_t;

    Inner* inner;

    ref_data(Wrapped& o) : inner(o._data().inner) {}

    static ref_data make(Inner* f) {return ref_data(f);}

private:
    ref_data(Inner* fp) : inner(fp) {}
};

template<class Wrapped, class Ref, class Inner>
struct srcref_data
{
    typedef void IS_REF_OR_CREF;
    typedef Wrapped wrapped_t;

    typedef const Inner* data_ref_t;
    typedef const Inner* data_srcref_t;

    const Inner* inner;

    srcref_data(const Wrapped& o) : inner(o._data().inner) {}
    srcref_data(Ref o) : inner(o._data().inner) {}

    static srcref_data make(const Inner* f) {return srcref_data(f);}

private:
    srcref_data(const Inner* fp) : inner(fp) {}
};

// Helper to determine if T is a flint class.
template<class T, class Enable = void> struct is_flint_class : mp::false_ { };
template<class T>
struct is_flint_class<T, typename T::IS_FLINT_CLASS> : mp::true_ { };

// From a lazy or immediate flint expression, obtain the evaluated
// non-reference type.
// Examples: fmpzxx -> fmpzxx
//           fmpzxx_ref -> fmpzxx
//           fmpzxx + fmpzxx -> fmpzxx
template<class T, class Enable = void>
struct to_nonref {typedef typename T::evaluated_t type;};

template<class T>
struct to_nonref<T, typename T::data_t::IS_REF_OR_CREF>
{
    typedef typename T::data_t::wrapped_t type;
};

template<class T>
struct c_base_t
{
    typedef typename T::c_base_t type;
};

// Given a lazy or non-lazy flint expression, obtain th evaluated reference
// type.
// Examples: fmpzxx -> fmpzxx_ref
//           fmpzxx_ref -> fmpzxx_ref
//           fmpzxx + fmpzxx -> fmpzxx_ref
template<class T>
struct to_ref
{
    typedef typename T::template make_helper<operations::immediate, ref_data<
        typename to_nonref<T>::type, typename c_base_t<T>::type> >::type type;
};
// Similarly for srcref.
template<class T>
struct to_srcref
{
    typedef typename T::template make_helper<operations::immediate, srcref_data<
        typename to_nonref<T>::type, typename to_ref<T>::type,
        typename c_base_t<T>::type> >::type type;
};

// Compute if Ref if the reference type belonging to compare.
// Examples: fmpzxx_ref, fmpzxx + fmpzxx -> true_
//           fmpzxx_srcref, fmpzxx -> false_
template<class Compare, class Ref>
struct is_ref : mp::equal_types<Ref, typename to_ref<Compare>::type> { };

// Similarly for srcref.
template<class Compare, class Ref>
struct is_srcref : mp::equal_types<Ref, typename to_srcref<Compare>::type> { };

// Similarly for non-ref.
template<class T>
struct is_nonref : mp::equal_types<T, typename to_nonref<T>::type > { };

// Flint classes allow implicit conversion only in very special situations.
// This template determines when. Currently, it is used exclusively to allow
// implicit conversion to reference types.
template<class T, class U, class Enable = void>
struct enableimplicit : mp::false_ { };

template<class T, class U>
struct enableimplicit<T, U, typename mp::enable_if<mp::and_<
        is_flint_class<T>, is_flint_class<U>
      > >::type>
    : mp::and_<
        traits::is_immediate_expr<T>,
        traits::is_immediate_expr<U>,
        mp::or_<
            mp::and_<is_ref<U, T>, is_nonref<U> >,
            mp::and_<is_srcref<U, T>, is_nonref<U> >,
            mp::and_<is_srcref<U, T>, is_ref<U, U> >
      > > { };

// Helper template which allows accessing data_(src)ref_t on immediates,
// without causing a compiler error on non-immediates.
// The main use for this are the _fmpz(), _fmpq() etc methods, which only
// work on immediates (but are defined on all instances).
template<class Expr, class Enable = void>
struct maybe_data_ref
{
    typedef void data_ref_t;
    typedef void data_srcref_t;
};
template<class Expr>
struct maybe_data_ref<Expr,
    typename mp::enable_if<
        // NB: cannot use is_immediate, since Expr may be incomplete type!
        mp::equal_types<typename Expr::operation_t, operations::immediate> >::type>
{
    typedef typename Expr::data_t::data_ref_t data_ref_t;
    typedef typename Expr::data_t::data_srcref_t data_srcref_t;
};

// If Base is a non-ref flint class, determine if T is a source operand
// (i.e. non-ref, ref or srcref type belong to Base)
// Examples: fmpzxx, fmpzxx_srcref -> true
//           fmpzxx, fmpzxx -> true
//           fmpzxx, fmpqxx -> false
template<class Base, class T, class Enable = void>
struct is_source : mp::false_ { };
template<class Base, class T>
struct is_source<Base, T, typename mp::enable_if<
    mp::and_<is_flint_class<T>, is_flint_class<Base> > >::type>
    : mp::or_<mp::equal_types<T, Base>, is_ref<Base, T>, is_srcref<Base, T> > { };

// Same with target (i.e. disallow srcref).
template<class Base, class T, class Enable = void>
struct is_target : mp::false_ { };
template<class Base, class T>
struct is_target<Base, T, typename mp::enable_if<
    mp::and_<is_flint_class<T>, is_flint_class<Base> > >::type>
    : mp::or_<mp::equal_types<T, Base>, is_ref<Base, T> > { };

// Predicate version of the above. Useful for FLINT_DEFINE_*_COND.
// See FLINTXX_COND_S and FLINTXX_COND_T
template<class Base>
struct is_source_base
{
    template<class T> struct type : is_source<Base, T> { };
};
template<class Base>
struct is_target_base
{
    template<class T> struct type : is_target<Base, T> { };
};

// Helper for implementing x += y*z etc
template<class T, class Right1, class Right2>
struct ternary_assign_helper
{
    typedef typename mp::make_tuple<Right1, Right2>::type tup_t;
    typedef tools::evaluate_n<tup_t> ev2_t;
    typedef typename ev2_t::temporaries_t temporaries_t;
    typedef mp::back_tuple<temporaries_t> back_t;

    typename back_t::type backing;
    ev2_t ev2;

    static temporaries_t backtemps(typename back_t::type& backing)
    {
        temporaries_t temps;
        back_t::init(temps, backing);
        return temps;
    }

    ternary_assign_helper(const tup_t& tup)
        : backing(mp::htuples::fill<typename back_t::type>(
                    tools::temporaries_filler(
                        tup.first()+tup.second() /* XXX */))),
          ev2(tup, backtemps(backing)) {}
    const T& getleft() {return ev2.template get<0>();}
    const T& getright() {return ev2.template get<1>();}
};

template<class T, class Right1, class Right2>
struct enable_ternary_assign
    : mp::enable_if<mp::and_<
          traits::is_T_expr<typename traits::basetype<Right1>::type, T>,
          traits::is_T_expr<typename traits::basetype<Right2>::type, T> >, T&> { };

// convenience helper
template<class Base, class T> struct is_Base : mp::or_<
     traits::is_T_expr<T, Base>,
     is_source<Base, T> > { };
} // flint_classes

namespace traits {
// Enable evaluation into reference types. See can_evaluate_into in
// expression_traits.h.
// XXX why do we need to disable the case where T, U are equal?
// Is <T, T> not more special?
template<class T, class U>
struct can_evaluate_into<T, U,
    typename mp::enable_if<mp::and_<flint_classes::is_flint_class<T>,
        flint_classes::is_flint_class<U>,
        mp::not_<mp::equal_types<T, U> > > >::type>
    : flint_classes::is_ref<U, T> { };
} // traits


namespace detail {
template<class Expr, class Enable = void>
struct should_enable_extra_ternop : mp::false_ { };
template<class Expr>
struct should_enable_extra_ternop<Expr, typename mp::enable_if<
    flint_classes::is_flint_class<Expr> >::type>
    : mp::equal_types<Expr, typename flint_classes::to_ref<
        typename flint_classes::to_nonref<Expr>::type>::type> { };
} // detail

// We add additional overloads for when the LHS is a reference type. The
// problem is that the standard overloads take LHS via reference, and rvalues
// (such as coming from fmpz_polyxx::get_coeff())
// cannot bind to this. In this case instead objects should be taken by value.
// However, this will make the overload ambiguous. Hence we take by const
// reference and then make an additional copy.
template<class Expr1, class Expr2>
inline typename mp::enable_if<detail::should_enable_extra_ternop<Expr1>, Expr1>::type
operator+=(const Expr1& e1, const Expr2& e2)
{
    Expr1(e1).set(e1 + e2);
    return e1;
}
template<class Expr1, class Expr2>
inline typename mp::enable_if<detail::should_enable_extra_ternop<Expr1>, Expr1>::type
operator-=(const Expr1& e1, const Expr2& e2)
{
    Expr1(e1).set(e1 - e2);
    return e1;
}
template<class Expr1, class Expr2>
inline typename mp::enable_if<detail::should_enable_extra_ternop<Expr1>, Expr1>::type
operator*=(const Expr1& e1, const Expr2& e2)
{
    Expr1(e1).set(e1 * e2);
    return e1;
}
template<class Expr1, class Expr2>
inline typename mp::enable_if<detail::should_enable_extra_ternop<Expr1>, Expr1>::type
operator/=(const Expr1& e1, const Expr2& e2)
{
    Expr1(e1).set(e1 / e2);
    return e1;
}
template<class Expr1, class Expr2>
inline typename mp::enable_if<detail::should_enable_extra_ternop<Expr1>, Expr1>::type
operator%=(const Expr1& e1, const Expr2& e2)
{
    Expr1(e1).set(e1 % e2);
    return e1;
}
template<class Expr1, class Expr2>
inline typename mp::enable_if<detail::should_enable_extra_ternop<Expr1>, Expr1>::type
operator<<=(const Expr1& e1, const Expr2& e2)
{
    Expr1(e1).set(e1 << e2);
    return e1;
}
template<class Expr1, class Expr2>
inline typename mp::enable_if<detail::should_enable_extra_ternop<Expr1>, Expr1>::type
operator>>=(const Expr1& e1, const Expr2& e2)
{
    Expr1(e1).set(e1 >> e2);
    return e1;
}
} // flint

// macros that help defining flint classes

#define FLINTXX_DEFINE_BASICS_NOFLINTCLASS(name)                              \
public:                                                                       \
    typedef typename base_t::evaluated_t evaluated_t;                         \
                                                                              \
    template<class T>                                                         \
    struct doimplicit                                                         \
        : flint_classes::enableimplicit<name, T> { };                         \
                                                                              \
    template<class T>                                                         \
    name& operator=(const T& t)                                               \
    {                                                                         \
        this->set(t);                                                         \
        return *this;                                                         \
    }                                                                         \
                                                                              \
protected:                                                                    \
    explicit name(const Data& d) : base_t(d) {}                               \
                                                                              \
    template<class D, class O, class Da>                                      \
    friend class expression;

// all flint classes should have this
#define FLINTXX_DEFINE_BASICS(name)                                           \
public:                                                                       \
    typedef void IS_FLINT_CLASS;                                              \
    FLINTXX_DEFINE_BASICS_NOFLINTCLASS(name)                                  \

// all flint classes should have this
#define FLINTXX_DEFINE_CTORS(name)                                            \
public:                                                                       \
    name() : base_t() {}                                                      \
    template<class T>                                                         \
    explicit name(const T& t,                                                 \
            typename mp::disable_if<doimplicit<T> >::type* = 0)               \
        : base_t(t) {}                                                        \
    template<class T>                                                         \
    explicit name(T& t,                                                       \
            typename mp::disable_if<doimplicit<T> >::type* = 0)               \
        : base_t(t) {}                                                        \
    template<class T>                                                         \
    name(const T& t,                                                          \
            typename mp::enable_if<doimplicit<T> >::type* = 0)                \
        : base_t(t) {}                                                        \
    template<class T>                                                         \
    name(T& t,                                                                \
            typename mp::enable_if<doimplicit<T> >::type* = 0)                \
        : base_t(t) {}                                                        \
    template<class T, class U>                                                \
    name(const T& t, const U& u) : base_t(t, u) {}                            \
    template<class T, class U>                                                \
    name(T& t, const U& u) : base_t(t, u) {}                            \
    template<class T, class U, class V>                                       \
    name(const T& t, const U& u, const V& v) : base_t(t, u, v) {}             \
    template<class T, class U, class V>                                       \
    name(T& t, const U& u, const V& v) : base_t(t, u, v) {}             \
    template<class T, class U, class V, class W>                              \
    name(const T& t, const U& u, const V& v, const W& w)                      \
        : base_t(t, u, v, w) {}                                               \
    template<class T, class U, class V, class W>                              \
    name(T& t, const U& u, const V& v, const W& w)                      \
        : base_t(t, u, v, w) {}

// Enable the flint reference type scheme. This typedefs c_base_t to ctype,
// and adds the data access wrapper (like _fmpz(), _fmpq()) called accessname.
// It also provides reference constructors from C types.
// All flint classes should have this.
#define FLINTXX_DEFINE_C_REF(name, ctype, accessname)                         \
public:                                                                       \
    typedef ctype c_base_t;                                                   \
    typedef flint_classes::maybe_data_ref<name> wrapped_traits;               \
    typename wrapped_traits::data_ref_t accessname()                          \
    {                                                                         \
        return this->_data().inner;                                           \
    }                                                                         \
    typename wrapped_traits::data_srcref_t accessname() const                 \
    {                                                                         \
        return this->_data().inner;                                           \
    }                                                                         \
                                                                              \
    /* These only make sense with the reference types */                      \
    template<class T>                                                         \
    static name make(T& f)                                                    \
    {                                                                         \
        return name(Data::make(f));                                           \
    }                                                                         \
    template<class T>                                                         \
    static name make(const T& f)                                              \
    {                                                                         \
        return name(Data::make(f));                                           \
    }                                                                         \
    template<class T, class U>                                                \
    static name make(T& f, const U& u)                                        \
    {                                                                         \
        return name(Data::make(f, u));                                        \
    }                                                                         \
    template<class T, class U>                                                \
    static name make(const T& f, const U& u)                                  \
    {                                                                         \
        return name(Data::make(f, u));                                        \
    }                                                                         \
    template<class T, class U, class V>                                       \
    static name make(const T& f, const U& u, const V& v)                      \
    {                                                                         \
        return name(Data::make(f, u, v));                                     \
    }

// Add a statically forwarded constructor called name. (Forwarded to data_t).
#define FLINTXX_DEFINE_FORWARD_STATIC(name)                                   \
    template<class T>                                                         \
    static typename base_t::derived_t name(const T& f)                        \
    {                                                                         \
        return typename base_t::derived_t(Data::name(f));                     \
    }                                                                         \
    template<class T, class U>                                                \
    static typename base_t::derived_t name(const T& f, const U& u)            \
    {                                                                         \
        return typename base_t::derived_t(Data::name(f, u));                  \
    }

// Add a static randomisation function.
// XXX this is not really useful because the arguments are often different.
#define FLINTXX_DEFINE_RANDFUNC(CBase, name) \
static CBase##xx_expression name(frandxx& state, mp_bitcnt_t bits) \
{ \
    CBase##xx_expression res; \
    CBase##_##name(res._data().inner, state._data(), bits); \
    return res; \
}

// Add a forwarded unary operation to the class. Suppose there is a unary
// operation foo() which returns my_typeA, and takes an argument of type
// my_typeB. Now on instances my_typeB, you want to write x.bar() for foo(x).
// Then add FLINTXX_DEFINE_MEMBER_UNOP_RTYPE_(my_typeA, bar, foo) to my_typeB.
//
// XXX due to circular definition problems, this cannot use the usual
// unary_op_helper type approach, and the unop must return the same type
// of expression
#define FLINTXX_DEFINE_MEMBER_UNOP_RTYPE_(rettype, name, funcname) \
FLINT_UNOP_BUILD_RETTYPE(funcname, rettype, typename base_t::derived_t) \
name() const \
{ \
    return flint::funcname(*this); \
}

// Convenience version when name==funcname
#define FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(rettype, name) \
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE_(rettype, name, name)

// Convenience version when rettype==argtype
#define FLINTXX_DEFINE_MEMBER_UNOP_(name, funcname) \
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE_(typename base_t::derived_t, name, funcname)

// Convenience version when rettype==argtype and name==funcname
#define FLINTXX_DEFINE_MEMBER_UNOP(name) FLINTXX_DEFINE_MEMBER_UNOP_(name, name)

// Add a forwarded binary operation. It is not necessary to specify the return
// type.
#define FLINTXX_DEFINE_MEMBER_BINOP_(name, funcname) \
template<class T> \
typename detail::binary_op_helper<typename base_t::derived_t, \
    operations::funcname##_op, T>::enable::type \
name(const T& t) const \
{ \
    return flint::funcname(*this, t); \
}

// Convenience version when funcname==name.
#define FLINTXX_DEFINE_MEMBER_BINOP(name) \
    FLINTXX_DEFINE_MEMBER_BINOP_(name, name)

#define FLINTXX_DEFINE_MEMBER_3OP_(name, funcname) \
template<class T, class U> \
typename detail::nary_op_helper2<operations::funcname##_op, \
    typename base_t::derived_t, T, U>::enable::type \
name(const T& t, const U& u) const \
{ \
    return flint::funcname(*this, t, u); \
}

#define FLINTXX_DEFINE_MEMBER_3OP(name) \
    FLINTXX_DEFINE_MEMBER_3OP_(name, name)

#define FLINTXX_DEFINE_MEMBER_4OP_(name, funcname) \
template<class T, class U, class V> \
typename detail::nary_op_helper2<operations::funcname##_op, \
    typename base_t::derived_t, T, U, V>::enable::type \
name(const T& t, const U& u, const V& v) const \
{ \
    return flint::funcname(*this, t, u, v); \
}

#define FLINTXX_DEFINE_MEMBER_4OP(name) \
    FLINTXX_DEFINE_MEMBER_4OP_(name, name)

#define FLINTXX_DEFINE_MEMBER_5OP_(name, funcname) \
template<class T, class U, class V, class W> \
typename detail::nary_op_helper2<operations::funcname##_op, \
    typename base_t::derived_t, T, U, V, W>::enable::type \
name(const T& t, const U& u, const V& v, const W& w) const \
{ \
    return flint::funcname(*this, t, u, v, w); \
}

#define FLINTXX_DEFINE_MEMBER_5OP(name) \
    FLINTXX_DEFINE_MEMBER_5OP_(name, name)

// Helper macros for FLINT_DEFINE_*_COND?.
#define FLINTXX_COND_S(Base) flint_classes::is_source_base<Base>::template type
#define FLINTXX_COND_T(Base) flint_classes::is_target_base<Base>::template type

// Convenience rules. These all take a Base class as argument, and will
// automatically apply to the related reference types as well.

// Add a to_string() conversion rule, empolying the common flint idiom where
// the string is allocated by the to_string function.
#define FLINTXX_DEFINE_TO_STR(Base, eval) \
template<class T> \
struct to_string<T, \
    typename mp::enable_if< FLINTXX_COND_S(Base)<T> >::type> \
{ \
    static std::string get(const T& from, int base) \
    { \
        char* str = eval; \
        std::string res(str); \
        flint_free(str); \
        return res; \
    } \
};

// Add a swap rule.
#define FLINTXX_DEFINE_SWAP(Base, eval) \
template<class T, class U> \
struct swap<T, U, typename mp::enable_if< mp::and_< \
    FLINTXX_COND_T(Base)<T>, FLINTXX_COND_T(Base)<U> > >::type> \
{ \
    static void doit(T& e1, U& e2) \
    { \
        eval; \
    } \
};

// Define a conversion rule through a default-constructed temporary object.
#define FLINTXX_DEFINE_CONVERSION_TMP(totype, Base, eval) \
template<class T> \
struct conversion<totype, T, \
    typename mp::enable_if< FLINTXX_COND_S(Base)<T> >::type> \
{ \
    static totype get(const T& from) \
    { \
        totype to; \
        eval; \
        return to; \
    } \
};

// Define a cmp rule.
#define FLINTXX_DEFINE_CMP(Base, eval) \
template<class T, class U> \
struct cmp<T, U, \
    typename mp::enable_if<mp::and_<FLINTXX_COND_S(Base)<T>, \
        FLINTXX_COND_S(Base)<U> > >::type> \
{ \
    static int get(const T& e1, const U& e2) \
    { \
        return eval; \
    } \
};

// Define an equals rule.
#define FLINTXX_DEFINE_EQUALS(Base, eval) \
template<class T, class U> \
struct equals<T, U, typename mp::enable_if<mp::and_< \
    FLINTXX_COND_S(Base)<T>, FLINTXX_COND_S(Base)<U> > >::type> \
{ \
    static bool get(const T& e1, const U& e2) \
    { \
        return eval; \
    } \
};

// Define a string assignment rule (c/f many polynomial classes).
#define FLINTXX_DEFINE_ASSIGN_STR(Base, eval) \
template<class T, class U> \
struct assignment<T, U, \
    typename mp::enable_if<mp::and_< \
        FLINTXX_COND_T(Base)<T>, traits::is_string<U> > >::type> \
{ \
    static void doit(T& to, const char* from) \
    { \
        eval; \
    } \
};

#define FLINTXX_UNADORNED_MAKETYPES(Base, left, op, right) \
    Base##_expression< op, tuple< left, tuple< right, empty_tuple> > >

// Optimized evaluation rules using ternary arithmetic (addmul, submul)
// NB: this has to be called in namespace flint, not flint::rules!
#define FLINTXX_DEFINE_TERNARY(Base, addmuleval, submuleval, maketypes)       \
namespace rules {                                                             \
/* a +- b*c */                                                                \
template<class Op, class Left, class Right1, class Right2>                    \
struct evaluation<Op,                                                         \
    tuple<Left, tuple<                                                        \
        maketypes(Base, Right1, operations::times, Right2),                   \
        /* NB: there is no particular reason to have the enable_if here,      \
               many other similar places would do */                          \
        typename mp::enable_if<mp::or_<                                       \
                mp::equal_types<Op, operations::plus>,                        \
                mp::equal_types<Op, operations::minus> >,                     \
            empty_tuple>::type> >,                                            \
    true, 1,                                                                  \
    typename tools::ternary_helper<Base, Left, Right1, Right2>::enable::type> \
{                                                                             \
    /* Helpful for testing. */                                                \
    static const unsigned TERNARY_OP_MARKER = 0;                              \
                                                                              \
    typedef Base return_t;                                                    \
    typedef tools::ternary_helper<Base, Left, Right1, Right2> th;             \
    typedef typename th::temporaries_t temporaries_t;                         \
    typedef tuple<Left, tuple<                                                \
        maketypes(Base, Right1, operations::times, Right2),                   \
        empty_tuple> > data_t;                                                \
    static const bool is_add = mp::equal_types<Op, operations::plus>::val;    \
                                                                              \
    static void doit(const data_t& input, temporaries_t temps, return_t* res) \
    {                                                                         \
        const Base* left = 0;                                                 \
        const Base* right = 0;                                                \
        th::doit(input.first(), input.second()._data().first(),               \
                input.second()._data().second(), temps, res, right, left);    \
        const Base& e1 = *left;                                               \
        const Base& e2 = *right;                                              \
        Base& to = *res;                                                      \
        if(is_add)                                                            \
        {                                                                     \
            addmuleval;                                                       \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            submuleval;                                                       \
        }                                                                     \
    }                                                                         \
};                                                                            \
                                                                              \
/* b*c + a */                                                                 \
template<class Right, class Left1, class Left2>                               \
struct evaluation<operations::plus,                                           \
    tuple<maketypes(Base, Left1, operations::times, Left2),                   \
        tuple<Right, empty_tuple> >,                                          \
    true, 1,                                                                  \
    typename tools::ternary_helper<Base,                                      \
        Right, Left1, Left2, operations::times>::enable::type>                \
{                                                                             \
    /* Helpful for testing. */                                                \
    static const unsigned TERNARY_OP_MARKER = 0;                              \
                                                                              \
    typedef Base return_t;                                                    \
    typedef tools::ternary_helper<Base, Right, Left1, Left2> th;              \
    typedef typename th::temporaries_t temporaries_t;                         \
    typedef tuple<maketypes(Base, Left1, operations::times, Left2),           \
        tuple<Right, empty_tuple> > data_t;                                   \
                                                                              \
    static void doit(const data_t& input, temporaries_t temps, return_t* res) \
    {                                                                         \
        const Base* left = 0;                                                 \
        const Base* right = 0;                                                \
        th::doit(input.second(), input.first()._data().first(),               \
                input.first()._data().second(), temps, res, right, left);     \
        const Base& e1 = *left;                                               \
        const Base& e2 = *right;                                              \
        Base& to = *res;                                                      \
        addmuleval;                                                           \
    }                                                                         \
};                                                                            \
} /* rules */                                                                 \
                                                                              \
/* TODO enable these with references on left hand side(?) */                  \
/* a += b*c */                                                                \
template<class Right1, class Right2>                                          \
inline typename flint_classes::enable_ternary_assign<Base, Right1, Right2>::type \
operator+=(Base& to,                                                          \
        const maketypes(Base, Right1, operations::times, Right2)& other)      \
{                                                                             \
    flint_classes::ternary_assign_helper<Base, Right1, Right2> tah(           \
            other._data());                                                   \
    const Base& e1 = tah.getleft();                                           \
    const Base& e2 = tah.getright();                                          \
    addmuleval;                                                               \
    return to;                                                                \
}                                                                             \
                                                                              \
/* a -= b*c */                                                                \
template<class Right1, class Right2>                                          \
inline typename flint_classes::enable_ternary_assign<Base, Right1, Right2>::type \
operator-=(Base& to,                                                          \
        const maketypes(Base, Right1, operations::times, Right2)& other)      \
{                                                                             \
    flint_classes::ternary_assign_helper<Base, Right1, Right2> tah(           \
            other._data());                                                   \
    const Base& e1 = tah.getleft();                                           \
    const Base& e2 = tah.getright();                                          \
    submuleval;                                                               \
    return to;                                                                \
}

#endif
