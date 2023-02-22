/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
     Demo FLINTXX header to illustrate flintxx extension.
*/

///////////////////////////////////////////////////////////////////////////////
// FAKE C DATA TYPE
// (This would normally reside in foo.h.)
///////////////////////////////////////////////////////////////////////////////

#ifndef FOO_H
#define FOO_H
#include <stdio.h>

extern "C" { // usually only #ifdef __cplusplus etc
typedef slong foo;
typedef slong foo_t[1];

static __inline__ void foo_init(foo_t f)
{
    *f = 0l;
}

static __inline__ void foo_clear(foo_t f)
{
}

static __inline__ void foo_set(foo_t to, const foo_t from)
{
    *to = *from;
}

static __inline__ void foo_set_si(foo_t f, slong e)
{
    *f = e;
}

static __inline__ void foo_add(foo_t to, const foo_t e1, const foo_t e2)
{
    *to = *e1 + *e2;
}

static __inline__ void foo_add_si(foo_t to, const foo_t e1, slong e2)
{
    *to = *e1 + e2;
}

static __inline__ int foo_cmp(const foo_t e1, const foo_t e2)
{
    if(*e1 == *e2)
        return 0;
    return *e1 > *e2 ? 1 : -1;
}

static __inline__ int foo_is_zero(const foo_t f)
{
    return *f == 0;
}

static __inline__ void foo_magic(foo_t to, const foo_t from)
{
    *to = 2 * (*from) + 1;
}
}

#endif

///////////////////////////////////////////////////////////////////////////////
// C++ wrapper
// (This would normally reside in fooxx.h.)
///////////////////////////////////////////////////////////////////////////////

#ifndef FOOXX_H
#define FOOXX_H

#include <iostream>

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"

namespace flint {
// fooxx_expression is an "all-purpose" expression template class. In
// principle, both Operation and Data can be arbitrary types (Data has to be
// copy constructible), but in this generality the objects will be not much
// use. In practice, Operation is an empty type, which is just used as a "tag",
// and Data is a rather primitive type holding essentially just some payload.
// Even more practically speaking, the only instantiations the FLINT developer
// should have have to make explicitly are when Operation is
// operations::immediate.
// The flintxx library will create other instantiations automatically, with
// more complicated Data arguments, and different Operation-s.
template<class Operation, class Data>
class fooxx_expression

// In order for the flintxx library to do its work, your class must derive from
// flint::expression. If your class has just two template parameters Operation
// and Data, then the following line is sufficient.
    : public expression<derived_wrapper<fooxx_expression>, Operation, Data>
{
public:

    // This line is formulaic, and just makes the base class available.
    // The typedef is used by the FLINTXX_DEFINE_* macros below, and is
    // necessary because of namespace injection bugs in gcc<4.5.
    typedef expression<derived_wrapper< ::flint::fooxx_expression>,
              Operation, Data> base_t;

    // The next two lines are formulaic, and most likely required in any
    // concrete class.
    FLINTXX_DEFINE_BASICS(fooxx_expression)
    FLINTXX_DEFINE_CTORS(fooxx_expression)

    // This line enables reference types for your class. The second argument is
    // the underlying C type (note this is foo, not foo_t). The third argument
    // is the name under which to make the underlying C type available.
    // All of fooxx, fooxx_ref and fooxx_srcref will have methods _foo() which
    // can be used to manipulate the underlying C data type.
    FLINTXX_DEFINE_C_REF(fooxx_expression, foo, _foo)

    // Now custom methods can be added. The typical pattern is to call a C
    // function with argument this->evaluate()._foo(). The evaluate() method is
    // inherited from the expression class (this is why it needs to be
    // qualified by "this"). It obtains a reference to self if self is an
    // immediate object, and otherwise evaluates self into a temporary
    // immediate object.
    // If you leave out the evaluate() step, then the method will only work
    // on immediates (which may be desirable).
    bool is_zero() const {return foo_is_zero(this->evaluate()._foo());}
};

// This is formulaic. The class fooxx will be an instantiation of
// fooxx_expression, with Operation operations::immediate and Data
// detail::foo_data. We need to forward-declare this because of cyclic
// dependencies among the immediate types (e.g. fooxx_srcref can be
// constructed from fooxx, and vice versa).
namespace detail {
struct foo_data;
}

// This line just carries out the plan of definition of fooxx explained above.
typedef fooxx_expression<operations::immediate, detail::foo_data> fooxx;

// If you want reference types (i.e. if you had FLINTXX_DEFINE_C_REF above),
// these lines are again formulaic.
typedef fooxx_expression<operations::immediate,
            flint_classes::ref_data<fooxx, foo> > fooxx_ref;
typedef fooxx_expression<operations::immediate,
            flint_classes::srcref_data<fooxx, fooxx_ref, foo> > fooxx_srcref;

namespace detail {

// We now define the actual immediate Data type. This is not just foo_t (the
// underlying C data type), because want it to behave nicely "in a C++ world".
struct foo_data
{
    // In general, your data type can contain members and member types in any
    // way you want. However, to work with the automatic reference type system,
    // the following three lines are necessary.
    foo_t inner;
    typedef foo_t& data_ref_t;
    typedef const foo_t& data_srcref_t;

    // Default constructor. If this is not provided, fooxx will not be default
    // constructible (this is OK but requires some additional care, see e.g.
    // padicxx).
    foo_data() {foo_init(inner);}

    // Destructor. You most likely want this.
    ~foo_data() {foo_clear(inner);}

    // Copy constructor. You must provide this.
    foo_data(const foo_data& o)
    {
        foo_init(inner);
        foo_set(inner, o.inner);
    }

    // Instantiation from srcref. This is basically the same as the copy,
    // constructor, but unfortunately has to be repeated. This also takes care
    // of instantiation from ref, since ref->srcref is an implicit conversion
    // path.
    foo_data(fooxx_srcref r)
    {
        foo_init(inner);
        foo_set(inner, r._foo());
    }

    // Now you can add more constructors, or in fact any methods you like.
    // This one allows constructing fooxx directly from long, int,
    // unsigned short etc.
    template<class T>
    foo_data(T t,
            typename mp::enable_if<traits::fits_into_slong<T> >::type* = 0)
    {
        foo_init(inner);
        foo_set_si(inner, t);
    }
};
} // detail

// By now our data type is instantiable, but nothing can be done with it.
// The flintxx library would be able to create expression templates involving
// fooxx, but will not do so because it has no way of evaluating them. We
// need to provides evaluation (and other) *rules* to the library. These
// (have to) live in namespace flint::rules.
//
// All possible rules are defined in flintxx/rules.h.

namespace rules {

// These two lines are convenient, are not formulaic except that they are used
// in all code below.
#define FOOXX_COND_S FLINTXX_COND_S(fooxx)
#define FOOXX_COND_T FLINTXX_COND_T(fooxx)

// Define a conditional assignment rule. The general pattern is
//
//     FLINT_DEFINE_DOIT_COND2(name, cond1, cond2, eval).
//
// This will define a "doit" rule for "name", which takes one input and
// one output argument. The result looks something like
//
//     template<class T, class U>
//     struct assignment<T, U, enable if cond1<T> and cond2<U> are satisfied>
//     {
//         static void doit(T& to, const U& from)
//         eval;
//     };
//
// In our case, we are defining an assignment rule, i.e. an explanation on
// how to execute operator=. If the right hand side is an expression template,
// flintxx will automatically evaluate it first. Thus we need only treat the
// case where the LHS is fmpzxx or fmpzxx_ref, and the RHS is fmpzxx, fmpzxx_ref
// or fmpzxx_srcref. This is precisely what the conditions FOOXX_COND_T
// and FOOXX_COND_S (conditions "fooxx target" and "fooxx source") mean.
FLINT_DEFINE_DOIT_COND2(assignment, FOOXX_COND_T, FOOXX_COND_S,
        foo_set(to._foo(), from._foo()))

// This line defines assignment of integral PODs to fooxx. Since the underlying
// C library only defines fooxx_set_si, we can only safely allow this if the
// right hand side can always be losslessly converted into a signed long,
// so we use the condition traits::fits_into_slong. Traits are defined all
// throughout flintxx, but the most general purpose ones (like fits_into_slong,
// is_unsigned_integer etc) can be found in flintxx/traits.h
FLINT_DEFINE_DOIT_COND2(assignment, FOOXX_COND_T, traits::fits_into_slong, 
        foo_set_si(to._foo(), from, 1))

// We now define evaluation rules. In full generality, the rule evaluation<...>
// can be used to define how to evaluate any kind of expression. But this is
// difficult to use. Moreover, much evaluation logic is shared among most
// data types. For example, to evaluate an expression like a + b + c,
// one typically first has to evaluate (say) a + b into a temporary t, and then
// evaluate t + c. The only step that is specific to fooxx here is how to
// add two immediates.
// For this reason, flintxx has special convenience forms of the evaluation
// rule, called binary and unary expressions. Defining a binary expression
// f(x, y) tells flintxx how to evaluate operation "f" on types "x" and "y",
// typically immediates. Then flintxx will figure out how to evaluate the
// arguments into temporaries first etc.
// There is a common special case, when f(x, y) is always the same as f(y, x),
// even though x and y may be of different types. Letting flintxx know of this
// avoids defining the rule both ways round.
//
// Here we define a commutative binary expression rule, for operation "plus",
// to be executed on to objects of types T and U, both satisfying FOOXX_COND_S.
// The result is to be of type foooxx (the second argument).
// In this case the types are fully symmetric, so we could have used
// FLINT_DEFINE_BINARY_EXPR_COND2 without adverse effects.
//
// The eval statement should have the effect of to = e1 + e2.
FLINT_DEFINE_CBINARY_EXPR_COND2(plus, fooxx, FOOXX_COND_S, FOOXX_COND_S,
        foo_add(to._foo(), e1._foo(), e2._foo()))

// Addation of fooxx and PODs. This time CBINARY instead of BINARY is vital.
FLINT_DEFINE_CBINARY_EXPR_COND2(plus, fooxx, FOOXX_COND_S,
        traits::fits_into_slong,
        foo_add_si(to._foo(), e1._foo(), e2))

// Next we define relational operators. A convenient way of doing so is using
// a "cmp" function, which is handily provided by the underlying C library.
// This has a somewhat peculiar signature, so cannot be defined using one of
// the standard macros. However, it comes up with many actual FLINT data types,
// so we have a special FLINTXX macro just for defining cmp.
FLINTXX_DEFINE_CMP(fooxx, foo_cmp(e1._foo(), e2._foo()))

// Now we define a rule how to print fooxx. There is no macro for this, because
// normally instead we define conversion to string, and flintxx takes care of
// printing. However, the C library for fooxx provides neither printing nor
// conversion to string, so we have to do our own implementation.
template<class T>
struct print<T, typename mp::enable_if<FOOXX_COND_S<T> >::type>
{
    static void doit(const T& i, std::ostream& o)
    {
        o << *i._foo();
    }
};
} // rules

// By now fooxx is a pretty passable wrapper type. In fact the only thing left
// to do is to expose foo_magic. This is a special function which can be
// executed on instances of foo, and yields another instance of foo. It is
// essentially just another unary expression, just with an unusual name, so
// this is how we treat it.

// This line introduces a new type of unary operation, called "magic_op",
// together with a function flint::magic(T), which creates expression templates
// with this new operation. In principle, any expression template data type is
// now allowed to define rules how to performa magic on itself.
FLINT_DEFINE_UNOP(magic)

// Finally, we need to explain how to perform magic on flintxx. This is again
// a rule.
namespace rules {

// The pattern should be familiar by now.
FLINT_DEFINE_UNARY_EXPR_COND(magic_op, fooxx, FOOXX_COND_S,
        foo_magic(to._foo(), from._foo()))
} // rules
} // flint

#endif

///////////////////////////////////////////////////////////////////////////////
// Example program
///////////////////////////////////////////////////////////////////////////////

using namespace flint;

int
main()
{
    fooxx a, b(4);
    fooxx_ref ar(a);
    fooxx_srcref br(b);

    ar = 1 + br + 1; // a=6
    std::cout << magic(a + (-1)) << '\n'; // 2*(6-1)+1 = 11

    return 0;
}
