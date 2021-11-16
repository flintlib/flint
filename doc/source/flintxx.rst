.. _flintxx:

**flintxx** -- C++ Wrapper
===============================================================================

Introduction
-------------------------------------------------------------------------------

Flint provides a C++ wrapper which makes extensive use of metaprogramming.
It is currently not maintained for the whole of Flint due to lack of a C++
developer willing to step up and maintain it. Therefore it is provided as-is,
though we do try to do sufficient maintenance to keep the
`currently wrapped functions <https://github.com/wbhart/flint2/blob/trunk/doc/source/flintxx_functions.txt>`_
compiling.

In this section we describe how to use ``flintxx`` the Flint C++ wrapper.

We begin with a simple example:

.. code-block:: C

    #include "fmpzxx.h"

    using namespace flint;
    fmpzxx x, y;
    x = 7u;
    y = x*x;
    std::cout << x << "^2 = " << y << std::endl; 

As can be seen, if a FLINT C interface is called ``foo``and resides in
``foo.h``, then the C++ version is called ``fooxx`` and resides in
``fooxx.h``. All flintxx classes live inside ``namespace flint``.

Functions which operate on wrapper classes are usually implemented as
overloaded stand-alone functions, with the type prefix dropped. For example,
a call to ``flint::gcd(f1, f2)`` yields an expression template evaluating via
``fmpz_gcd``, provided ``f1`` and ``f2`` evaluate to instances of
``fmpzxx``.

Sometimes we felt that dropping the type prefix would yield incomprehensible
names, as for example in ``fmpq_next_minimal``, or ``fmpq_reconstruct``. In
these cases the type prefix is swapped for the flintxx equivalent, so the
flintxx version would be called ``fmpqxx_reconstruct``.

In this situation, usually the same functionality is also exposed as a
(possibly static) member function, and this is the preferred way of
accessing the functionality. Thus one should write
``fmpqxx::reconstruct(a, m)`` or ``fmpqxx(0, 1u).next_minimal()``.

Expression templates
-------------------------------------------------------------------------------

The implementation of flintxx tries very hard not to incur any overhead over
using the native C interface. For this reason, we use ``expression templates``
for lazily evaluating expressions. This allows us to avoid creating
excessively many temporaries, for example.

This means that even if ``x`` and ``y`` are of type ``fmpzxx``, then ``x + y``
will not be of type ``fmpzxx``. Instead it will be an object which for most
purposes behaves just like ``fmpzxx``, but really only expresses the fact
that it represents the sum of ``x`` and ``y``.

This distinction almost never matters, since expression templates are evaluated
automatically in most cases. Thus ``cout << x + y`` or ``x + y == 7`` will
work just as one might expect.

There are ways to request explicit evaluation of an expression template, most
notably ``(x + y).evaluate()`` and ``fmpzxx(x + y)``.

One caveat of the expression template approach is that compiler error messages
can be long and hard to understand.

In flintxx we strive to type-check template parameters as early as possible in
order to keep error messages short and close to the actual user error.
Excessively long error messages are often indicative of a bug in flintxx.

Tuples
-------------------------------------------------------------------------------

Many FLINT functions naturally return two or more arguments. A typical example
is ``divrem``. The underlying C function is
``void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)``.

Mapping this directly to C++ would yield something like
``void divrem(fmpz_polyxx& Q, fmpz_polyxx& R, const fmpz_polyxx& A, const fmpz_polyxx& B)}``.

While functional, this is not particularly nice; the syntax
``divrem(Q, R, A, B)``, where the first two arguments are modified, is just
very reminiscent of C. We prefer an expression closer to the python analogue
``(Q, R) = divrem(A, B)``.

For this purpose, flintxx uses lazy tuples and the following is a valid
flintxx expression: ``ltupleref(Q, R) = divrem(A, B)``.

For the purpose of this documentation, ltuple types are denoted
as follows: ``Ltuple<Type1, Type2, ..., Typer>``.

Thus, ``divrem`` would return an object of type
``Ltuple<fmpz_polyxx, fmpz_polyxx>``.

The user should never try to construct such types by hand; instead use the
function ``ltupleref`` (and perhaps occasionally ``ltuple``; both documented
later).

One thing to observe is that ltuples are typed fairly weakly. Thus assignments
and equality comparisons can be performed as long as both sides have the same
length, and the operation can be performed on all components (whether or not
the component types match).

Another interesting feature of ltuples is the type
``flint::detail::IGNORED_TYPE``. In an ltuple assignment, where the left hand
side contains a reference to this type, the relevant entry is just discarded.

Since the ``ltuple.h`` header automatically creates a static instance ``_`` of
this type, in the following listing, the lines marked (1) and (2) have the same
effect (but the second is potentially more efficient).

.. code-block:: C

    #include "fmpz_polyxx.h"

    using namespace flint;

    fmpz_polyxx f, g;

    fmpz_polyxx R;
    ltupleref(_, R) = divrem(f, g); // (1)
    R = f % g;                      // (2)

Note finally that using ``ltuple`` intermediates often results in more
copies than necessary. For example the expression
``ltupleref(num, _) = divrem(a, b)`` assigns the quotient to ``num``,
creating just a temporary ``fmpzxx`` to hold the remainder. In contrast,
``num = divrem(a, b).get<0>()`` creates two temporary instances of
``fmpzxx``.

Reference types
-------------------------------------------------------------------------------

One subtlety in wrapping a C library is that references do not work as easily
as one might expect. For example, consider the class ``fmpqxx``, wrapping
``fmpq_t``, i.e. rational numbers. As such, an instance of ``fmpqxx`` has a
numerator and denominator. In C, these are accessible via macros
``fmpq_numref`` and ``fmpq_denref``, which yield ``fmpz*``, which can be used
essentially interchangeably with ``fmpz_t``. In particular, any library
function which operates on ``fmpz_t`` can operate on the numerator or
denominator of an ``fmpq_t``. In C++, we would like to have a member functions
``den` and ``num`` which return an object of type ``fmpzxx&`` (i.e.
a reference to ``fmpzxx``).

However, this is not possible, since ``fmpqxx`` is not implemented as a pair
of ``fmpzxx``, and instead simply contains an ``fmpq_t``.

For this reason, for every C interface ``foo``, flintxx provides two
additional types, called ``fooxx_ref`` and ``fooxx_srcref``, acting as
replacements for ``fooxx&`` and ``const foox&``, respectively, in
situations where no underlying C++ object exists.

Instances of ``fooxx_ref`` or ``fooxx_srcref`` behave exactly like instances
of ``fooxx``. In fact, the user should never notice a difference. Any flintxx
operation or expression which works on objects of type ``foo`` also works on
objects of type ``fooxx_ref`` and ``fooxx_srcref``.

Moreover, instances of ``foo`` can be converted implicitly to ``fooxx_ref``
or ``fooxx_srcref``, and ``fooxx_ref`` can be converted implicitly to
``fooxx_srcref``.

It is also possible to explicitly convert reference types ``fooxx_*ref`` to
``fooxx`` (since this entails copying, we provide no implicit conversion).

In summary, the class ``fooxx_ref`` behaves like a reference to an object of
type ``fooxx``. As such it can be used both as a right hand side and as a
left hand side, just like ``fooxx``.

The class ``fooxx_srcref`` behaves like a reference to a constant object of
type ``fooxx``, and so cannot be used as a left hand side. These objects are
created by flintxx automatically,`for example upon calling
``fmpqxx::num()``.

Unified coefficient access
-------------------------------------------------------------------------------

Consider again the ``x.num()`` method of ``fmpqxx``. In various situations,
this can have different return types. Namely, if ``x`` is a writable
expression, then ``x.num()`` returns an ``fmpzxx_ref``. In particular the
return value behaves just like ``fmpzxx``, no evaluation is necessary to
obtain it, there are no copies, and it is possible to change the
return value (and thus change ``x``).

If on the other hand ``x`` is a readonly immediate, then the return value of
``x.num()`` has type ``fmpzxx_srcref``. This again behaves just like
``fmpzxx`` and no evaluations or copies are necessary, but this time it is
not possible to change the return value (and so it is not possible to change
``x``, either).

Finally, if ``x`` is a lazy expression, then the return value is actually a
lazy expression template. Thus to obtain the "actual" value of ``x.num()``,
evaluations are necessary, and potentially so are copies.

Thus in any case the return value behaves just like ``fmpqxx``, but apart
from that the behaviour of ``x.num()`` varies quite drastically in the
different situations. We call this "unified coefficient access" (the
coefficients of a ``fmpqxx`` being ``num(), den()``), and the same
behaviour occurs in many other flintxx types, e.g. in
``fmpz_polyxx.coeff()``, etc.

Type conversion
-------------------------------------------------------------------------------

As a rule, flintxx does not perform automatic type conversions (except when
related to the promotion to reference types, c/f earlier discussion). In
expression templates, operands can be automatically promoted if the underlying
C interface provides this facility. Beyond that, types have to be converted
explicitly.

There are two ways of doing this. The preferred one is using static
constructor functions. Typical examples are
``fmpz_polyxx::from_ground(fmpzarg)`` and
``nmod_polyxx::reduce(mplimbarg, nmodctxarg)``. The former takes an (expression
template evaluating to) ``fmpzxx`` and returns an ``fmpz_polyxx`` representing
the constant polynomial with value the ``fmpzxx``. The latter takes an argument
of type ``mp_limb_t`` and one of type ``nmodxx_ctx_srcref`` (essentially a
word-sized modulus) and returns an ``nmod_polyxx`` representing the constant
polynomial obtained by reducing ``mplimbarg``.

The general format for this is ``totype::constructorname(arg1, arg2, ...)``.
We prefer this because it makes explicit the type that is being converted to,
and the way the arguments are to be interpreted.

This format only works if the target type is part of flintxx. In other cases,
we use a ``.to<totype>()`` syntax, as in ``fmpzexpr.to<slong>()``.

Input and output
-------------------------------------------------------------------------------

In C++ it is customary to provide input and output via iostreams, and
overloading the operators ``<<`` and ``>>``. When wrapping a C library which
works on the ``FILE`` interface, this is rather hard to accomplish.

For this reason, flintxx only provides streaming output (i.e. ``<<``), and
only when there is a ``to_string`` method. Unfortunately this applies to only
a small subset of the FLINT types.

For output in other cases, and input in all cases, we provide C-like functions.
Namely, the functions ``print``, ``print_pretty``, ``read`` and ``read_pretty``
can be used similarly to the C ``flint_printf`` and ``scanf``.

For example, ``print(x)`` where ``x`` is an ``fmpz`` has the same effect as
``std::cout << x``.

Inheritance and flintxx
-------------------------------------------------------------------------------

The flintxx classes are not designed for inheritance. If you want to modify
behaviour, you should wrap flintxx types into your own classes (extension by
aggregation, not inheritance).

Notation and conventions in flintxx documentation
-------------------------------------------------------------------------------

As explained above, the flintxx classes and functions perform quite a number of
operations which should be invisible to the user. Some template types implement
methods which only make sense for some template arguments, etc.

For example, every expression template built from ``fmpq_polyxx`` (polynomials
with rational coefficients) has a method ``set_coeff``. However, this method
only makes sense for objects of type ``fmpq_polyxx`` or ``fmpq_polyxx_ref``
(calling it on other types will result in a compilation error), and its
existence in objects of other types should be considered an implementation
detail.

In what follows, we document a "virtual" set of classes and functions, which
explain how the user should expect its objects to behave, and which we
guarantee to maintain. Other interfaces should be considered implementation
details and subject to change.

Consider the interface ``fmpzxx``, and more concretely an instance ``a``.
As in the above discussion, we see that from ``a`` we can build a lot of
different objects: expression templates like ``a+a``, constant objects like
``const fmpzxx& b = a;``, reference objects like ``fmpzxx_ref c(a)``, etc.
These by nature behave somewhat differently. For our purposes, we classify
types into "targets" (things which can be assigned to), "sources" (things
which contain actual computed data, or references thereto, as opposed to lazy
expression templates) and "expressions" (sources or expression templates).

Note that every target is a source, and every source is an expression.

We denote any type which can act as a target for ``fmpzxx`` as ``Fmpz_target``
(note the initial capital letter!), any ``fmpzxx`` source as ``Fmpz_source``
and any ``fmpzxx`` expression as ``Fmpz_expr``. Such made up type names
(always with initial capital letter) are referred to as "virtual types" in the
documentation. These are used for all flint classes (e.g. ``Fmpq_expr`` or
``Fmpz_polyxx_src``).

When using virtual types, we will suppress reference notation. No flintxx types
are ever copied automatically, unless the documentation explicitly says so.
This is a general philosophy of flintxx: the library does as many things
automatically as it can, without introducing additional calls to underlying
Flint C functions. So for example, it is not possible to implicitly convert
``int`` to ``fmpzxx`` (since doing so requires a C call). Of course explicit
conversions (or assignments) work completely fine.

It is also often the case that flintxx functions are conditionally enabled
templates. A notation such as ``void foo(T:is_signed_integer)`` denotes a
template function which is enabled whenever the template parameter ``T``
satisfies the type trait ``is_signed_integer``. These type traits should be
self-explanatory.

In what follows, we will never document copy constructors, or implicit
conversion constructors pertaining to reference types. We will also not
document assignment operators for expressions of the same type. Thus if
``x`` is an ``fmpzxx`` and ``y`` is an ``fmpqxx``, then ``x = x`` and
``y = x`` are both valid, but only the second assignment operator is
documented explicitly.

Most flintxx functions and methods wrap underlying C functions in a way which
is evident from the signature of the flintxx function/method. If this is the
case, no further documentation is provided. For example, the function
``double dlog(Fmpz_expr x)`` simply wraps ``double fmpz_dlog(const fmpz_t)``. 

As is evident from the return type, ``dlog`` immediately evaluates its
argument, and then computes the logarithm. In contrast, a function like
``Fmpz_expr gcd(Fmpz_expr, Fmpz_expr)`` returns a lazily evaluated expression
template and wraps ``void fmpz_gcd(fmpz_t, const fmpz_t, const fmpz_t)``.

In case a Flint C function has more than one return value in the form of
arguments passed in by reference, the C++ wrapper returns an ``ltuple``. In
this case, the order of the ``ltuple`` arguments is the same as the order of
the function arguments; so for example ``ltupleref(Q, R) = divrem(A, B)`` has
the same effect as ``fmpz_poly_divrem(q, r, a, b)``, provided ``Q, R, A, B``
are ``fmpz_polyxx`` and ``q, r, a, b`` are the underlying ``fmpz_poly_t``.

If such a convention is followed, the documentation below may not further
explain anything. In all other cases, further explanation is provided (this
applies in particular if the C function has return type different from
``void``).

Global functions or member functions?
-------------------------------------------------------------------------------

Often it is not clear if functionality is exposed as a global function,
such as ``gcd(a, b)``, or as a member function, such as ``a.gcd(b)``. In
flintxx, we strive to make both available when feasible.

In the documentation, the global versions are documented in detail (explaining
the allowed types etc), whereas the member function versions are summarised
more briefly under e.g. ``Fmpz_expr::unary operation() const``,
``Fmpz_expr::binary operation(??) const`` etc.




