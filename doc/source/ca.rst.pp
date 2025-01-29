.. _ca:

**ca.h** -- exact real and complex numbers
===============================================================================

A :type:`ca_t` represents a real or complex number in a form suitable
for exact field arithmetic or comparison.
Exceptionally, a :type:`ca_t` may represent a special nonnumerical value,
such as an infinity.

Introduction: numbers
-------------------------------------------------------------------------------

A *Calcium number* is a real or complex number
represented as an element of a formal field `K = \mathbb{Q}(a_1, \ldots, a_n)`
where the symbols `a_k` denote fixed algebraic or transcendental numbers
called *extension numbers*.
For example, `e^{-2 \pi} - 3 i` may be represented as
`(1 - 3 a_2^2 a_1) / a_2^2` in the field `\mathbb{Q}(a_1,a_2)` with
`a_1 = i, a_2 = e^{\pi}`.
Extension numbers and fields are documented
in the following separate modules:

* :ref:`ca-ext`

* :ref:`ca-field`

The user does not need to construct extension numbers or formal
extension fields explicitly: each :type:`ca_t` contains an internal
pointer to its formal field, and operations on Calcium numbers generate
and cache fields automatically as needed to express the results.

This representation is not canonical (in general). A given complex number
can be represented in different ways depending on
the choice of formal field *K*. Even within a fixed field *K*,
a number can have different representations if there are algebraic
relations between the extension numbers.
Two numbers *x* and *y* can be tested for inequality using numerical
evaluation; to test for equality, it may be necessary to
eliminate dependencies between extension numbers.
One of the central goals of Calcium will be to
implement heuristics for such elimination.

Together with each formal field *K*, Calcium stores a
*reduction ideal* `I = \{g_1,\ldots,g_m\}`
with `g_i \in \mathbb{Z}[a_1,\ldots,a_n]`, defining a set of algebraic relations
`g_i(a_1,\ldots,a_n) = 0`.
Relations can be absolute, say `g_i = a_j^2 + 1`, or relative,
say `g_i = 2 a_j - 4 a_k - a_l a_m`.
The reduction ideal effectively partitions `K` into
equivalence classes of
complex numbers
(e.g. `i^2 = -1` or `2 \log(\pi i) = 4 \log(\sqrt{\pi}) + \pi i`),
enabling simplifications and equality proving.

Extension numbers are always sorted `a_1 \succ a_2 \succ \ldots \succ a_n`
where `\succ` denotes a structural ordering (see :func:`ca_cmp_repr`).
If the reduction ideal is triangular and the multivariate polynomial
arithmetic uses lexicographic ordering, reduction by *I*
eliminates numbers `a_i` with higher complexity in the sense of `\succ`.

The reduction ideal is an imperfect computational crutch: it is not guaranteed
to capture *all* algebraic relations, and reduction is not guaranteed
to produce uniquely defined representatives.
However, in the specific case of an absolute number field `K = \mathbb{Q}(a)`
where *a* is a :type:`qqbar_t` extension, the reduction ideal (consisting
of a single minimal polynomial) is canonical and field
elements of *K* can be chosen canonically.

Introduction: special values
-------------------------------------------------------------------------------

In order to provide a closed arithmetic system and express limiting
cases of operators and special functions, a :type:`ca_t` can hold
any of the following special values besides ordinary numbers:

* *Unsigned infinity*, a formal object `{\tilde \infty}` representing
  the value of `1 / 0`. More generally, this is the
  value of meromorphic functions at poles.

* *Signed infinity*, a formal object `a \cdot \infty` where the sign `a`
  is a Calcium number with `|a| = 1`.
  The most common values are `+\infty, -\infty, +i \infty, -i \infty`.
  Signed infinities are used to denote directional limits and logarithmic
  singularities (for example, `\log(0) = -\infty`).

* *Undefined*, a formal object representing the value of indeterminate
  forms such as `0 / 0` and essential singularities such as 
  `\exp(\tilde \infty)`, where a number or infinity would not make sense
  as an answer.

* *Unknown*, a meta-value used to signal that the actual desired value
  could not be computed, either because Calcium does not (yet) have a data
  structure or algorithm for that case, or because doing so would be
  unreasonably expensive. This occurs, for example, if Calcium performs
  a division and is unable to decide whether the result is a number,
  unsigned infinity or undefined (because testing for zero fails).
  Wrappers may want to check output variables for
  *Unknown* and throw an exception (e.g. *NotImplementedError* in Python).

The distinction between *Calcium numbers* (which must represent
elements of `\mathbb{C}`) and the different kinds of nonnumerical values
(infinities, Undefined or Unknown) is essential. Nonnumerical values may
not be used as field extension numbers `a_k`, and the denominator of
a formal field element must always represent a nonzero complex number.
Accordingly, for any given Calcium value *x* that is not *Unknown*,
it is exactly known whether *x* represents A) a number, B) unsigned infinity,
C) a signed infinity, or D) Undefined.

Number objects
-------------------------------------------------------------------------------

For all types, a *type_t* is defined as an array of length one of type
*type_struct*, permitting a *type_t* to be passed by reference.

.. type:: ca_struct

.. type:: ca_t

    A :type:`ca_t` contains an index to a field *K*, and
    data representing an element *x* of *K*.
    The data is either an inline rational number (:type:`fmpq_t`), an
    inline Antic number field element (:type:`nf_elem_t`)
    when *K* is an absolute algebraic number field `\mathbb{Q}(a)`,
    or a pointer to a heap-allocated :type:`fmpz_mpoly_q_t` representing
    an element of a generic field `\mathbb{Q}(a_1,\ldots,a_n)`.
    Special values are encoded using magic bits in the field index.

.. type:: ca_ptr

   Alias for ``ca_struct *``, used for vectors of numbers.

.. type:: ca_srcptr

   Alias for ``const ca_struct *``, used for vectors of numbers
   when passed as constant input to functions.


Context objects
-------------------------------------------------------------------------------

.. type:: ca_ctx_struct

.. type:: ca_ctx_t

    A :type:`ca_ctx_t` context object holds a cache of fields *K* and
    constituent extension numbers `a_k`.
    The field index in an individual :type:`ca_t` instance represents
    a shallow reference to the object defining the field *K* within the
    context object, so creating many elements of the same field is cheap.

    Since context objects are mutable (and may be mutated even when
    performing read-only operations on :type:`ca_t` instances), they must not
    be accessed simultaneously by different threads: in multithreaded
    environments, the user must use a separate context object for each thread.

.. function:: void ca_ctx_init(ca_ctx_t ctx)

    Initializes the context object *ctx* for use.
    Any evaluation options stored in the context object
    are set to default values.

.. function:: void ca_ctx_clear(ca_ctx_t ctx)

    Clears the context object *ctx*, freeing any memory allocated internally.
    This function should only be called after all :type:`ca_t` instances
    referring to this context have been cleared.

.. function:: void ca_ctx_print(ca_ctx_t ctx)

    Prints a description of the context *ctx* to standard output.
    This will give a complete listing of the cached fields in *ctx*.

Memory management for numbers
-------------------------------------------------------------------------------

.. function:: void ca_init(ca_t x, ca_ctx_t ctx)

    Initializes the variable *x* for use, associating it with the
    context object *ctx*. The value of *x* is set to the rational number 0.

.. function:: void ca_clear(ca_t x, ca_ctx_t ctx)

    Clears the variable *x*.

.. function:: void ca_swap(ca_t x, ca_t y, ca_ctx_t ctx)

    Efficiently swaps the variables *x* and *y*.


Symbolic expressions
-------------------------------------------------------------------------------

.. function:: void ca_get_fexpr(fexpr_t res, const ca_t x, ulong flags, ca_ctx_t ctx)

    Sets *res* to a symbolic expression representing *x*.

.. function:: int ca_set_fexpr(ca_t res, const fexpr_t expr, ca_ctx_t ctx)

    Sets *res* to the value represented by the symbolic expression *expr*.
    Returns 1 on success and 0 on failure.
    This function essentially just traverses the expression tree using
    ``ca`` arithmetic; it does not provide advanced symbolic evaluation.
    It is guaranteed to at least be able to parse the output of
    :func:`ca_get_fexpr`.


.. _ca-printing:

Printing
-------------------------------------------------------------------------------

The style of printed output is controlled by
``ctx->options[CA_OPT_PRINT_FLAGS]``
(see :ref:`context-options`) which can be set to any
combination of the following flags:

.. macro:: CA_PRINT_N

    Print a decimal approximation of the number.
    The approximation is guaranteed to be correctly rounded to within
    one unit in the last place.

    If combined with ``CA_PRINT_REPR``, numbers appearing
    within the symbolic representation will also be printed with
    decimal approximations.

    Warning: printing a decimal approximation requires a computation,
    which can be expensive. It can also mutate
    cached data (numerical enclosures of extension numbers), affecting
    subsequent computations.

.. macro:: CA_PRINT_DIGITS

    Multiplied by a positive integer, specifies the number of
    decimal digits to show with ``CA_PRINT_N``. If not given,
    the default precision is six digits.

.. macro:: CA_PRINT_REPR

    Print the symbolic representation of the number (including
    its recursive elements). If used together with ``CA_PRINT_N``,
    field elements will print as ``decimal {symbolic}`` while
    extension numbers will print as ``decimal [symbolic]``.

    All extension numbers appearing in the field defining ``x``
    and in the inner constructions of those extension numbers
    will be given local labels ``a``, ``b``, etc. for this printing.

.. macro:: CA_PRINT_FIELD

    For each field element, explicitly print its formal field
    along with its reduction ideal if present, e.g. ``QQ`` or
    ``QQ(a,b,c) / <a-b, c^2+1>``.

.. macro:: CA_PRINT_DEFAULT

    The default print style. Equivalent to ``CA_PRINT_N | CA_PRINT_REPR``.

.. macro:: CA_PRINT_DEBUG

    Verbose print style for debugging. Equivalent to ``CA_PRINT_N | CA_PRINT_REPR | CA_PRINT_FIELD``.

As a special case, small integers are always printed
as simple literals.

As illustration, here are the numbers
`-7`, `2/3`, `(\sqrt{3}+5)/2` and `\sqrt{2} (\log(\pi) + \pi i)`
printed in various styles::

    # CA_PRINT_DEFAULT
    -7
    0.666667 {2/3}
    3.36603 {(a+5)/2 where a = 1.73205 [a^2-3=0]}
    1.61889 + 4.44288*I {a*c+b*c*d where a = 1.14473 [Log(3.14159 {b})], b = 3.14159 [Pi], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}

    # CA_PRINT_N
    -7
    0.666667
    3.36603
    1.61889 + 4.44288*I

    # CA_PRINT_N | (CA_PRINT_DIGITS * 20)
    -7
    0.66666666666666666667
    3.3660254037844386468
    1.6188925298220266685 + 4.4428829381583662470*I

    # CA_PRINT_REPR
    -7
    2/3
    (a+5)/2 where a = [a^2-3=0]
    a*c+b*c*d where a = Log(b), b = Pi, c = [c^2-2=0], d = [d^2+1=0]

    # CA_PRINT_DEBUG
    -7
    0.666667 {2/3  in  QQ}
    3.36603 {(a+5)/2  in  QQ(a)/<a^2-3> where a = 1.73205 [a^2-3=0]}
    1.61889 + 4.44288*I {a*c+b*c*d  in  QQ(a,b,c,d)/<c^2-2, d^2+1> where a = 1.14473 [Log(3.14159 {b  in  QQ(b)})], b = 3.14159 [Pi], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}

.. function:: void ca_print(const ca_t x, ca_ctx_t ctx)

    Prints *x* to standard output.

.. function:: void ca_fprint(FILE * fp, const ca_t x, ca_ctx_t ctx)

    Prints *x* to the file *fp*.

.. function:: char * ca_get_str(const ca_t x, ca_ctx_t ctx)

    Prints *x* to a string which is returned.
    The user should free this string by calling ``flint_free``.

.. function:: void ca_printn(const ca_t x, slong n, ca_ctx_t ctx)

    Prints an *n*-digit numerical representation of *x* to standard output.

Special values
-------------------------------------------------------------------------------

.. function:: void ca_zero(ca_t res, ca_ctx_t ctx)
              void ca_one(ca_t res, ca_ctx_t ctx)
              void ca_neg_one(ca_t res, ca_ctx_t ctx)

    Sets *res* to the integer 0, 1 or -1. This creates a canonical representation
    of this number as an element of the trivial field `\mathbb{Q}`.

.. function:: void ca_i(ca_t res, ca_ctx_t ctx)
              void ca_neg_i(ca_t res, ca_ctx_t ctx)

    Sets *res* to the imaginary unit `i = \sqrt{-1}`, or its negation `-i`.
    This creates a canonical representation of `i` as the generator of the
    algebraic number field `\mathbb{Q}(i)`.

.. function:: void ca_pi(ca_t res, ca_ctx_t ctx)

    Sets *res* to the constant `\pi`. This creates an element
    of the transcendental number field `\mathbb{Q}(\pi)`.

.. function:: void ca_pi_i(ca_t res, ca_ctx_t ctx)

    Sets *res* to the constant `\pi i`. This creates an element of the
    composite field `\mathbb{Q}(i,\pi)` rather than representing `\pi i`
    (or even `2 \pi i`, which for some purposes would be more elegant)
    as an atomic quantity.

.. function:: void ca_euler(ca_t res, ca_ctx_t ctx)

    Sets *res* to Euler's constant `\gamma`. This creates an element
    of the (transcendental?) number field `\mathbb{Q}(\gamma)`.

.. function:: void ca_unknown(ca_t res, ca_ctx_t ctx)

    Sets *res* to the meta-value *Unknown*.

.. function:: void ca_undefined(ca_t res, ca_ctx_t ctx)

    Sets *res* to *Undefined*.

.. function:: void ca_uinf(ca_t res, ca_ctx_t ctx)

    Sets *res* to unsigned infinity `{\tilde \infty}`.

.. function:: void ca_pos_inf(ca_t res, ca_ctx_t ctx)
              void ca_neg_inf(ca_t res, ca_ctx_t ctx)
              void ca_pos_i_inf(ca_t res, ca_ctx_t ctx)
              void ca_neg_i_inf(ca_t res, ca_ctx_t ctx)

    Sets *res* to the signed infinity `+\infty`, `-\infty`, `+i \infty` or `-i \infty`.

Assignment and conversion
-------------------------------------------------------------------------------

.. function:: void ca_set(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to a copy of *x*.

.. function:: void ca_set_si(ca_t res, slong v, ca_ctx_t ctx)
              void ca_set_ui(ca_t res, ulong v, ca_ctx_t ctx)
              void ca_set_fmpz(ca_t res, const fmpz_t v, ca_ctx_t ctx)
              void ca_set_fmpq(ca_t res, const fmpq_t v, ca_ctx_t ctx)

    Sets *res* to the integer or rational number *v*. This creates a canonical
    representation of this number as an element of the trivial field
    `\mathbb{Q}`.

.. function:: void ca_set_d(ca_t res, double x, ca_ctx_t ctx)
              void ca_set_d_d(ca_t res, double x, double y, ca_ctx_t ctx)

    Sets *res* to the value of *x*, or the complex value `x + yi`.
    NaN is interpreted as *Unknown* (not *Undefined*).

.. function:: void ca_transfer(ca_t res, ca_ctx_t res_ctx, const ca_t src, ca_ctx_t src_ctx)

    Sets *res* to *src* where the corresponding context objects *res_ctx* and
    *src_ctx* may be different.

    This operation preserves the mathematical value represented by *src*,
    but may result in a different internal representation depending on the
    settings of the context objects.

Conversion of algebraic numbers
-------------------------------------------------------------------------------

.. function:: void ca_set_qqbar(ca_t res, const qqbar_t x, ca_ctx_t ctx)

    Sets *res* to the algebraic number *x*.

    If *x* is rational, *res* is set to the canonical representation as
    an element in the trivial field `\mathbb{Q}`.

    If *x* is irrational, this function always sets *res* to an element of
    a univariate number field `\mathbb{Q}(a)`. It will not, for example,
    identify `\sqrt{2} + \sqrt{3}`
    as an element of `\mathbb{Q}(\sqrt{2}, \sqrt{3})`. However, it may
    attempt to find a simpler number field than that generated by *x*
    itself. For example:

    * If *x* is quadratic, it will be expressed as an element of
      `\mathbb{Q}(\sqrt{N})` where *N* has no small repeated factors
      (obtained by performing a smooth factorization of the discriminant).

    * TODO: if possible, coerce *x* to a low-degree cyclotomic field.

.. function:: int ca_get_fmpz(fmpz_t res, const ca_t x, ca_ctx_t ctx)
              int ca_get_fmpq(fmpq_t res, const ca_t x, ca_ctx_t ctx)
              int ca_get_qqbar(qqbar_t res, const ca_t x, ca_ctx_t ctx)

    Attempts to evaluate *x* to an explicit integer, rational or
    algebraic number. If successful, sets *res* to this number and
    returns 1. If unsuccessful, returns 0.

    The conversion certainly fails if *x* does not represent an integer,
    rational or algebraic number (respectively), but can also fail if *x*
    is too expensive to compute under
    the current evaluation limits.
    In particular, the evaluation will be aborted if an intermediate
    algebraic number (or more precisely, the resultant polynomial prior
    to factorization) exceeds ``CA_OPT_QQBAR_DEG_LIMIT``
    or the coefficients exceed some multiple of ``CA_OPT_PREC_LIMIT``.
    Note that evaluation may hit those limits even if the minimal polynomial
    for *x* itself is small. The conversion can also fail if no algorithm
    has been implemented for the functions appearing in the construction
    of *x*.

.. function:: int ca_can_evaluate_qqbar(const ca_t x, ca_ctx_t ctx)

    Checks if :func:`ca_get_qqbar` has a chance to succeed. In effect,
    this checks if all extension numbers are manifestly algebraic
    numbers (without doing any evaluation).

Random generation
-------------------------------------------------------------------------------

.. function:: void ca_randtest_rational(ca_t res, flint_rand_t state, slong bits, ca_ctx_t ctx)

    Sets *res* to a random rational number with numerator and denominator
    up to *bits* bits in size.

.. function:: void ca_randtest(ca_t res, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx)

    Sets *res* to a random number generated by evaluating a random expression.
    The algorithm randomly selects between generating a "simple" number
    (a random rational number or quadratic field element with coefficients
    up to *bits* in size, or a random builtin constant),
    or if *depth* is nonzero, applying a random arithmetic operation or
    function to operands produced through recursive calls with
    *depth* - 1. The output is guaranteed to be a number, not a special value.

.. function:: void ca_randtest_special(ca_t res, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx)

    Randomly generates either a special value or a number.

.. function:: void ca_randtest_same_nf(ca_t res, flint_rand_t state, const ca_t x, slong bits, slong den_bits, ca_ctx_t ctx)

    Sets *res* to a random element in the same number field as *x*,
    with numerator coefficients up to *bits* in size and denominator
    up to *den_bits* in size. This function requires that *x* is an
    element of an absolute number field.

Representation properties
-------------------------------------------------------------------------------

The following functions deal with the representation of a :type:`ca_t` and
hence can always be decided quickly and unambiguously. The return value
for predicates is 0 for false and 1 for true.

.. function:: int ca_equal_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Returns whether *x* and *y* have identical representation. For field
    elements, this checks if *x* and *y* belong to the same formal field
    (with generators having identical representation) and are represented by
    the same rational function within that field.

    For special values, this tests
    equality of the special values, with *Unknown* handled as if it were
    a value rather than a meta-value: that is, *Unknown* = *Unknown* gives 1,
    and *Unknown* = *y* gives 0 for any other kind of value *y*.
    If neither *x* nor *y* is *Unknown*, then representation equality
    implies that *x* and *y* describe to the same mathematical value, but if
    either operand is *Unknown*, the result is meaningless for
    mathematical comparison.

.. function:: int ca_cmp_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Compares the representations of *x* and *y* in a canonical sort order,
    returning -1, 0 or 1. This only performs a lexicographic comparison
    of the representations of *x* and *y*; the return value does not say
    anything meaningful about the numbers represented by *x* and *y*.

.. function:: ulong ca_hash_repr(const ca_t x, ca_ctx_t ctx)

    Hashes the representation of *x*.

.. function:: int ca_is_unknown(const ca_t x, ca_ctx_t ctx)

    Returns whether *x* is Unknown.

.. function:: int ca_is_special(const ca_t x, ca_ctx_t ctx)

    Returns whether *x* is a special value or metavalue
    (not a field element).

.. function:: int ca_is_qq_elem(const ca_t x, ca_ctx_t ctx)

    Returns whether *x* is represented as an element of the
    rational field `\mathbb{Q}`.

.. function:: int ca_is_qq_elem_zero(const ca_t x, ca_ctx_t ctx)
              int ca_is_qq_elem_one(const ca_t x, ca_ctx_t ctx)
              int ca_is_qq_elem_integer(const ca_t x, ca_ctx_t ctx)

    Returns whether *x* is represented as the element 0, 1 or
    any integer in the rational field `\mathbb{Q}`.

.. function:: int ca_is_nf_elem(const ca_t x, ca_ctx_t ctx)

    Returns whether *x* is represented as an element of a univariate
    algebraic number field `\mathbb{Q}(a)`.

.. function:: int ca_is_cyclotomic_nf_elem(slong * p, ulong * q, const ca_t x, ca_ctx_t ctx)

    Returns whether *x* is represented as an element of a univariate
    cyclotomic field, i.e. `\mathbb{Q}(a)` where *a* is a root of unity.
    If *p* and *q* are not *NULL* and *x* is represented as an
    element of a cyclotomic field, this also sets *p* and *q* to the
    minimal integers with `0 \le p < q` such that the generating
    root of unity is `a = e^{2 \pi i p / q}`.
    Note that the answer 0 does not prove that *x* is not
    a cyclotomic number,
    and the order *q* is also not necessarily the generator of the
    *smallest* cyclotomic field containing *x*.
    For the purposes of this function, only nontrivial
    cyclotomic fields count; the return value is 0 if *x*
    is represented as a rational number.

.. function:: int ca_is_generic_elem(const ca_t x, ca_ctx_t ctx)

    Returns whether *x* is represented as a generic field element;
    i.e. it is not a special value, not represented as
    an element of the rational field, and not represented as
    an element of a univariate algebraic number field.


Value predicates
-------------------------------------------------------------------------------

The following predicates check a mathematical property which might
not be effectively decidable. The result is a :type:`truth_t` to allow
representing an unknown outcome.

.. function:: truth_t ca_check_is_number(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a number. The result is ``T_TRUE`` is *x* is
    a field element (and hence a complex number), ``T_FALSE`` if *x* is
    an infinity or *Undefined*, and ``T_UNKNOWN`` if *x* is *Unknown*.

.. function:: truth_t ca_check_is_zero(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_one(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_one(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_i(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_i(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is equal to the number `0`, `1`, `-1`, `i`, or `-i`.

.. function:: truth_t ca_check_is_algebraic(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_rational(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_integer(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is respectively an algebraic number, a rational number,
    or an integer.

.. function:: truth_t ca_check_is_real(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a real number. Warning: this returns ``T_FALSE`` if *x* is an
    infinity with real sign.

.. function:: truth_t ca_check_is_negative_real(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a negative real number. Warning: this returns ``T_FALSE``
    if *x* is negative infinity.

.. function:: truth_t ca_check_is_imaginary(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is an imaginary number. Warning: this returns ``T_FALSE`` if
    *x* is an infinity with imaginary sign.

.. function:: truth_t ca_check_is_undefined(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is the special value *Undefined*.

.. function:: truth_t ca_check_is_infinity(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is any infinity (unsigned or signed).

.. function:: truth_t ca_check_is_uinf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is unsigned infinity `{\tilde \infty}`.

.. function:: truth_t ca_check_is_signed_inf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is any signed infinity.

.. function:: truth_t ca_check_is_pos_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_pos_i_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_i_inf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is equal to the signed infinity `+\infty`, `-\infty`,
    `+i \infty`, `-i \infty`, respectively.

Comparisons
-------------------------------------------------------------------------------

.. function:: truth_t ca_check_equal(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Tests `x = y` as a mathematical equality.
    The result is ``T_UNKNOWN`` if either operand is *Unknown*.
    The result may also be ``T_UNKNOWN`` if *x* and *y* are numerically
    indistinguishable and cannot be proved equal or unequal by
    an exact computation.

.. function:: truth_t ca_check_lt(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_le(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_gt(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_ge(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Compares *x* and *y*, implementing the respective operations
    `x < y`, `x \le y`, `x > y`, `x \ge y`.
    Only real numbers and `-\infty` and `+\infty` are considered comparable.
    The result is ``T_FALSE`` (not ``T_UNKNOWN``) if either operand is not
    comparable (being a nonreal complex number, unsigned infinity, or
    undefined).

Field structure operations
-------------------------------------------------------------------------------

.. function:: void ca_merge_fields(ca_t resx, ca_t resy, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *resx* and *resy* to copies of *x* and *y* coerced to a common field.
    Both *x* and *y* must be field elements (not special values).

    In the present implementation, this simply merges the lists of generators,
    avoiding duplication. In the future, it will be able to eliminate
    generators satisfying algebraic relations.

.. function:: void ca_condense_field(ca_t res, ca_ctx_t ctx)

    Attempts to demote the value of *res* to a trivial subfield of its
    current field by removing unused generators. In particular, this demotes
    any obviously rational value to the trivial field `\mathbb{Q}`.

    This function is applied automatically in most operations
    (arithmetic operations, etc.).

.. function:: ca_ext_ptr ca_is_gen_as_ext(const ca_t x, ca_ctx_t ctx)

    If *x* is a generator of its formal field, `x = a_k \in \mathbb{Q}(a_1,\ldots,a_n)`,
    returns a pointer to the extension number defining `a_k`. If *x* is
    not a generator, returns *NULL*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void ca_neg(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the negation of *x*.
    For numbers, this operation amounts to a direct negation within the
    formal field.
    For a signed infinity `c \infty`, negation gives `(-c) \infty`; all other
    special values are unchanged.

.. function:: void ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_add_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_add_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_add(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the sum of *x* and *y*. For special values, the following
    rules apply (`c \infty` denotes a signed infinity, `|c| = 1`):

    * `c \infty + d \infty = c \infty` if `c = d`

    * `c \infty + d \infty = \text{Undefined}` if `c \ne d`

    * `\tilde \infty + c \infty = \tilde \infty + \tilde \infty = \text{Undefined}`

    * `c \infty + z = c \infty` if `z \in \mathbb{C}`

    * `\tilde \infty + z = \tilde \infty` if `z \in \mathbb{C}`

    * `z + \text{Undefined} = \text{Undefined}` for any value *z* (including *Unknown*)

    In any other case involving special values, or if the specific case cannot
    be distinguished, the result is *Unknown*.

.. function:: void ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_sub_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_sub_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_fmpq_sub(ca_t res, const fmpq_t x, const ca_t y, ca_ctx_t ctx)
              void ca_fmpz_sub(ca_t res, const fmpz_t x, const ca_t y, ca_ctx_t ctx)
              void ca_ui_sub(ca_t res, ulong x, const ca_t y, ca_ctx_t ctx)
              void ca_si_sub(ca_t res, slong x, const ca_t y, ca_ctx_t ctx)
              void ca_sub(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the difference of *x* and *y*.
    This is equivalent to computing `x + (-y)`.

.. function:: void ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_mul_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_mul_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_mul(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the product of *x* and *y*. For special values, the following
    rules apply (`c \infty` denotes a signed infinity, `|c| = 1`):

    * `c \infty \cdot d \infty = c d \infty`

    * `c \infty \cdot \tilde \infty = \tilde \infty`

    * `\tilde \infty \cdot \tilde \infty = \tilde \infty`

    * `c \infty \cdot z = \operatorname{sgn}(z) c \infty` if `z \in \mathbb{C} \setminus \{0\}`

    * `c \infty \cdot 0 = \text{Undefined}`

    * `\tilde \infty \cdot 0 = \text{Undefined}`

    * `z \cdot  \text{Undefined} = \text{Undefined}` for any value *z* (including *Unknown*)

    In any other case involving special values, or if the specific case cannot
    be distinguished, the result is *Unknown*.

.. function:: void ca_inv(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the multiplicative inverse of *x*. In a univariate algebraic
    number field, this always produces a rational denominator, but the
    denominator might not be rationalized in a multivariate
    field.
    For special values
    and zero, the following rules apply:

    * `1 / (c \infty) = 1 / \tilde \infty = 0`

    * `1 / 0 = \tilde \infty`

    * `1 / \text{Undefined} = \text{Undefined}`

    * `1 / \text{Unknown} = \text{Unknown}`

    If it cannot be determined whether *x* is zero or nonzero,
    the result is *Unknown*.

.. function:: void ca_fmpq_div(ca_t res, const fmpq_t x, const ca_t y, ca_ctx_t ctx)
              void ca_fmpz_div(ca_t res, const fmpz_t x, const ca_t y, ca_ctx_t ctx)
              void ca_ui_div(ca_t res, ulong x, const ca_t y, ca_ctx_t ctx)
              void ca_si_div(ca_t res, slong x, const ca_t y, ca_ctx_t ctx)
              void ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_div_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_div_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_div_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_div(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the quotient of *x* and *y*. This is equivalent
    to computing `x \cdot (1 / y)`. For special values and division
    by zero, this implies the following rules
    (`c \infty` denotes a signed infinity, `|c| = 1`):

    * `(c \infty) / (d \infty) = (c \infty) / \tilde \infty = \tilde \infty / (c \infty) = \tilde \infty / \tilde \infty = \text{Undefined}`

    * `c \infty / z = (c / \operatorname{sgn}(z)) \infty` if `z \in \mathbb{C} \setminus \{0\}`

    * `c \infty / 0 = \tilde \infty / 0 = \tilde \infty`

    * `z / (c \infty) = z / \tilde \infty = 0` if `z \in \mathbb{C}`

    * `z / 0 = \tilde \infty` if `z \in \mathbb{C} \setminus \{0\}`

    * `0 / 0 = \text{Undefined}`

    * `z / \text{Undefined} = \text{Undefined}` for any value *z* (including *Unknown*)

    * `\text{Undefined} / z = \text{Undefined}` for any value *z* (including *Unknown*)

    In any other case involving special values, or if the specific case cannot
    be distinguished, the result is *Unknown*.

.. function:: void ca_dot(ca_t res, const ca_t initial, int subtract, ca_srcptr x, slong xstep, ca_srcptr y, slong ystep, slong len, ca_ctx_t ctx)

    Computes the dot product of the vectors *x* and *y*, setting
    *res* to `s + (-1)^{subtract} \sum_{i=0}^{len-1} x_i y_i`.

    The initial term *s* is optional and can be
    omitted by passing *NULL* (equivalently, `s = 0`).
    The parameter *subtract* must be 0 or 1.
    The length *len* is allowed to be negative, which is equivalent
    to a length of zero.
    The parameters *xstep* or *ystep* specify a step length for
    traversing subsequences of the vectors *x* and *y*; either can be
    negative to step in the reverse direction starting from
    the initial pointer.
    Aliasing is allowed between *res* and *s* but not between
    *res* and the entries of *x* and *y*.

.. function:: void ca_fmpz_poly_evaluate(ca_t res, const fmpz_poly_t poly, const ca_t x, ca_ctx_t ctx)
              void ca_fmpq_poly_evaluate(ca_t res, const fmpq_poly_t poly, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the polynomial *poly* evaluated at *x*.

.. function:: void ca_fmpz_mpoly_evaluate_horner(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)
              void ca_fmpz_mpoly_evaluate_iter(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)
              void ca_fmpz_mpoly_evaluate(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)

    Sets *res* to the multivariate polynomial *f* evaluated at the vector of arguments *x*.

.. function:: void ca_fmpz_mpoly_q_evaluate(ca_t res, const fmpz_mpoly_q_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)

    Sets *res* to the multivariate rational function *f* evaluated at the vector of arguments *x*.

.. function:: void ca_fmpz_mpoly_q_evaluate_no_division_by_zero(ca_t res, const fmpz_mpoly_q_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)
              void ca_inv_no_division_by_zero(ca_t res, const ca_t x, ca_ctx_t ctx)

    These functions behave like the normal arithmetic functions,
    but assume (and do not check) that division by zero cannot occur.
    Division by zero will result in undefined behavior.


Powers and roots
-------------------------------------------------------------------------------

.. function:: void ca_sqr(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the square of *x*.

.. function:: void ca_pow_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_pow_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_pow_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_pow_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_pow(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to *x* raised to the power *y*.
    Handling of special values is not yet implemented.

.. function:: void ca_pow_si_arithmetic(ca_t res, const ca_t x, slong n, ca_ctx_t ctx)

    Sets *res* to *x* raised to the power *n*. Whereas :func:`ca_pow`,
    :func:`ca_pow_si` etc. may create `x^n` as an extension number
    if *n* is large, this function always perform the exponentiation
    using field arithmetic.

.. function:: void ca_sqrt_inert(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_sqrt_nofactor(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_sqrt_factor(ca_t res, const ca_t x, ulong flags, ca_ctx_t ctx)
              void ca_sqrt(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the principal square root of *x*.

    For special values, the following definitions apply:

    * `\sqrt{c \infty} = \sqrt{c} \infty`

    * `\sqrt{\tilde \infty} = \tilde \infty`.

    * Both *Undefined* and *Unknown* map to themselves.

    The *inert* version outputs the generator in the formal field
    `\mathbb{Q}(\sqrt{x})` without simplifying.

    The *factor* version writes `x = A^2 B` in `K` where `K` is
    the field of *x*, and outputs `A \sqrt{B}` or
    `-A \sqrt{B}` (whichever gives the correct sign) as an element of
    `K(\sqrt{B})` or some subfield thereof.
    This factorization is only a heuristic and is not guaranteed
    to make `B` minimal.
    Factorization options can be passed through to *flags*: see
    :func:`ca_factor` for details.

    The *nofactor* version will not perform a general factorization, but
    may still perform other simplifications. It may in particular attempt to
    simplify `\sqrt{x}` to a single element in `\overline{\mathbb{Q}}`.

.. function:: void ca_sqrt_ui(ca_t res, ulong n, ca_ctx_t ctx)

    Sets *res* to the principal square root of *n*.


Complex parts
-------------------------------------------------------------------------------

.. function:: void ca_abs(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the absolute value of *x*.

    For special values, the following definitions apply:

    * `|c \infty| = |\tilde \infty| = +\infty`.

    * Both *Undefined* and *Unknown* map to themselves.

    This function will attempt to simplify its argument through an exact
    computation. It may in particular attempt to simplify `|x|` to
    a single element in `\overline{\mathbb{Q}}`.

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(|x|)`.

.. function:: void ca_sgn(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the sign of *x*, defined by

    .. math::

        \operatorname{sgn}(x) = \begin{cases} 0 & x = 0 \\ \frac{x}{|x|} & x \ne 0 \end{cases}

    for numbers. For special values, the following definitions apply:

    * `\operatorname{sgn}(c \infty) = c`.

    * `\operatorname{sgn}(\tilde \infty) = \operatorname{Undefined}`.

    * Both *Undefined* and *Unknown* map to themselves.

    This function will attempt to simplify its argument through an exact
    computation. It may in particular attempt to simplify `\operatorname{sgn}(x)` to
    a single element in `\overline{\mathbb{Q}}`.

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(\operatorname{sgn}(x))`.

.. function:: void ca_csgn(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the extension of the real sign function taking the
    value 1 for *z* strictly in the right half plane, -1 for *z* strictly
    in the left half plane, and the sign of the imaginary part when *z* is on
    the imaginary axis. Equivalently, `\operatorname{csgn}(z) = z / \sqrt{z^2}`
    except that the value is 0 when *z* is exactly zero.
    This function gives *Undefined* for unsigned infinity
    and `\operatorname{csgn}(\operatorname{sgn}(c \infty)) = \operatorname{csgn}(c)` for
    signed infinities.

.. function:: void ca_arg(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the complex argument (phase) of *x*, normalized
    to the range `(-\pi, +\pi]`. The argument of 0 is defined as 0.
    For special values, the following definitions apply:

    * `\operatorname{arg}(c \infty) = \operatorname{arg}(c)`.

    * `\operatorname{arg}(\tilde \infty) = \operatorname{Undefined}`.

    * Both *Undefined* and *Unknown* map to themselves.

.. function:: void ca_re(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the real part of *x*. The result is *Undefined* if *x*
    is any infinity (including a real infinity).

.. function:: void ca_im(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the imaginary part of *x*. The result is *Undefined* if *x*
    is any infinity (including an imaginary infinity).

.. function:: void ca_conj_deep(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_conj_shallow(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_conj(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the complex conjugate of *x*.
    The *shallow* version creates a new extension element
    `\overline{x}` unless *x* can be trivially conjugated in-place
    in the existing field.
    The *deep* version recursively conjugates the extension numbers
    in the field of *x*.

.. function:: void ca_floor(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the floor function of *x*. The result is *Undefined* if *x*
    is any infinity (including a real infinity).
    For complex numbers, this is presently defined to take the floor of the
    real part.

.. function:: void ca_ceil(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the ceiling function of *x*. The result is *Undefined* if *x*
    is any infinity (including a real infinity).
    For complex numbers, this is presently defined to take the ceiling of the
    real part.


Exponentials and logarithms
-------------------------------------------------------------------------------

.. function:: void ca_exp(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the exponential function of *x*.

    For special values, the following definitions apply:

    * `e^{+\infty} = +\infty`

    * `e^{c \infty} = \tilde \infty` if `0 < \operatorname{Re}(c) < 1`.

    * `e^{c \infty} = 0` if `\operatorname{Re}(c) < 0`.

    * `e^{c \infty} = \text{Undefined}` if `\operatorname{Re}(c) = 0`.

    * `e^{\tilde \infty} = \text{Undefined}`.

    * Both *Undefined* and *Unknown* map to themselves.

    The following symbolic simplifications are performed automatically:

    * `e^0 = 1`

    * `e^{\log(z)} = z`

    * `e^{(p/q) \log(z)} = z^{p/q}` (for rational `p/q`)

    * `e^{(p/q) \pi i}` = algebraic root of unity (for small rational `p/q`)

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(e^x)`.

.. function:: void ca_log(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the natural logarithm of *x*.

    For special values and at the origin, the following definitions apply:

    * For any infinity, `\log(c\infty) = \log(\tilde \infty) = +\infty`.

    * `\log(0) = -\infty`. The result is *Unknown* if deciding `x = 0` fails.

    * Both *Undefined* and *Unknown* map to themselves.

    The following symbolic simplifications are performed automatically:

    * `\log(1) = 0`

    * `\log\left(e^z\right) = z + 2 \pi i k`

    * `\log\left(\sqrt{z}\right) = \tfrac{1}{2} \log(z) + 2 \pi i k`

    * `\log\left(z^a\right) = a \log(z) + 2 \pi i k`

    * `\log(x) = \log(-x) + \pi i` for negative real *x*

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(\log(x))`.


Trigonometric functions
-------------------------------------------------------------------------------

.. function:: void ca_sin_cos_exponential(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
              void ca_sin_cos_direct(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
              void ca_sin_cos_tangent(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
              void ca_sin_cos(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)

    Sets *res1* to the sine of *x* and *res2* to the cosine of *x*.
    Either *res1* or *res2* can be *NULL* to compute only the other function.
    Various representations are implemented:

    * The *exponential* version expresses the sine and cosine in terms
      of complex exponentials. Simple algebraic values will simplify to
      rational numbers or elements of cyclotomic fields.

    * The *direct* method expresses the sine and cosine in terms of
      the original functions (perhaps after applying some symmetry
      transformations, which may interchange sin and cos). Extremely
      simple algebraic values will automatically simplify to elements
      of real algebraic number fields.

    * The *tangent* version expresses the sine and cosine in terms
      of `\tan(x/2)`, perhaps after applying some symmetry
      transformations. Extremely simple algebraic values will
      automatically simplify to elements of real algebraic number
      fields.

    By default, the standard function uses the *exponential*
    representation as this typically works best for field arithmetic
    and simplifications, although it has the disadvantage of
    introducing complex numbers where real numbers would be sufficient.
    The behavior of the standard function can be changed using the
    :macro:`CA_OPT_TRIG_FORM` context setting.

    For special values, the following definitions apply:

    * `\sin(\pm i \infty) = \pm i \infty`

    * `\cos(\pm i \infty) = +\infty`

    * All other infinities give `\operatorname{Undefined}`

.. function:: void ca_sin(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_cos(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the sine or cosine of *x*.
    These functions are shortcuts for :func:`ca_sin_cos`.

.. function:: void ca_tan_sine_cosine(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_tan_exponential(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_tan_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_tan(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the tangent of *x*.
    The *sine_cosine* version evaluates the tangent as a quotient of
    a sine and cosine, the *direct* version evaluates it directly
    as a tangent (possibly after transforming the variable),
    and the *exponential* version evaluates it in terms of complex
    exponentials.
    Simple algebraic values will automatically simplify to elements
    of trigonometric or cyclotomic number fields.

    By default, the standard function uses the *exponential*
    representation as this typically works best for field arithmetic
    and simplifications, although it has the disadvantage of
    introducing complex numbers where real numbers would be sufficient.
    The behavior of the standard function can be changed using the
    :macro:`CA_OPT_TRIG_FORM` context setting.

    For special values, the following definitions apply:

    * At poles, `\tan((n+\tfrac{1}{2}) \pi) = \tilde \infty`

    * `\tan(e^{i \theta} \infty) = +i, \quad 0 < \theta < \pi`

    * `\tan(e^{i \theta} \infty) = -i, \quad -\pi < \theta < 0`

    * `\tan(\pm \infty) = \tan(\tilde \infty) = \operatorname{Undefined}`


.. function:: void ca_cot(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the cotangent *x*. This is equivalent to
    computing the reciprocal of the tangent.

.. function:: void ca_atan_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_atan_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_atan(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the inverse tangent of *x*.

    The *direct* version expresses the result as an inverse tangent
    (possibly after transforming the variable). The *logarithm*
    version expresses it in terms of complex logarithms.
    Simple algebraic inputs will automatically simplify to
    rational multiples of `\pi`.

    By default, the standard function uses the *logarithm*
    representation as this typically works best for field arithmetic
    and simplifications, although it has the disadvantage of
    introducing complex numbers where real numbers would be sufficient.
    The behavior of the standard function can be changed using the
    :macro:`CA_OPT_TRIG_FORM` context setting (exponential mode
    results in logarithmic forms).

    For special values, the following definitions apply:

    * `\operatorname{atan}(\pm i) = \pm i \infty`

    * `\operatorname{atan}(c \infty) = \operatorname{csgn}(c) \pi / 2`

    * `\operatorname{atan}(\tilde \infty) = \operatorname{Undefined}`

.. function:: void ca_asin_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_acos_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_asin_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_acos_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_asin(ca_t res, const ca_t x, ca_ctx_t ctx)
              void ca_acos(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the inverse sine (respectively, cosine) of *x*.

    The *direct* version expresses the result as an inverse sine or
    cosine (possibly after transforming the variable). The *logarithm*
    version expresses it in terms of complex logarithms.
    Simple algebraic inputs will automatically simplify to
    rational multiples of `\pi`.

    By default, the standard function uses the *logarithm*
    representation as this typically works best for field arithmetic
    and simplifications, although it has the disadvantage of
    introducing complex numbers where real numbers would be sufficient.
    The behavior of the standard function can be changed using the
    :macro:`CA_OPT_TRIG_FORM` context setting (exponential mode
    results in logarithmic forms).

    The inverse cosine is presently implemented as
    `\operatorname{acos}(x) = \pi/2 - \operatorname{asin}(x)`.


Special functions
-------------------------------------------------------------------------------

.. function:: void ca_gamma(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the gamma function of *x*.

.. function:: void ca_erf(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the error function of *x*.

.. function:: void ca_erfc(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the complementary error function of *x*.

.. function:: void ca_erfi(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the imaginary error function of *x*.


Numerical evaluation
-------------------------------------------------------------------------------

.. function:: void ca_get_acb_raw(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)

    Sets *res* to an enclosure of the numerical value of *x*.
    A working precision of *prec* bits is used internally for the evaluation,
    without adaptive refinement.
    If *x* is any special value, *res* is set to *acb_indeterminate*.

.. function:: void ca_get_acb(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)
              void ca_get_acb_accurate_parts(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)

    Sets *res* to an enclosure of the numerical value of *x*.
    The working precision is increased adaptively to try to ensure *prec*
    accurate bits in the output. The *accurate_parts* version tries to ensure
    *prec* accurate bits for both the real and imaginary part separately.

    The refinement is stopped if the working precision exceeds
    ``CA_OPT_PREC_LIMIT`` (or twice the initial precision, if this is larger).
    The user may call *acb_rel_accuracy_bits* to check is the calculation
    was successful.

    The output is not rounded down to *prec* bits (to avoid unnecessary
    double rounding); the user may call *acb_set_round* when rounding
    is desired.

.. function:: char * ca_get_decimal_str(const ca_t x, slong digits, ulong flags, ca_ctx_t ctx)

    Returns a decimal approximation of *x* with precision up to
    *digits*. The output is guaranteed to be correct within 1 ulp
    in the returned digits, but the number of returned digits may
    be smaller than *digits* if the numerical evaluation does
    not succeed.

    If *flags* is set to 1, attempts to achieve full accuracy for
    both the real and imaginary parts separately.

    If *x* is not finite or a finite enclosure cannot be produced,
    returns the string "?".

    The user should free the returned string with ``flint_free``.


Rewriting and simplification
-------------------------------------------------------------------------------

.. function:: void ca_rewrite_complex_normal_form(ca_t res, const ca_t x, int deep, ca_ctx_t ctx)

    Sets *res* to *x* rewritten using standardizing transformations
    over the complex numbers:

    * Elementary functions are rewritten in terms of (complex) exponentials, roots and logarithms
    * Complex parts are rewritten using logarithms, square roots, and (deep) complex conjugates
    * Algebraic numbers are rewritten in terms of cyclotomic fields where applicable

    If *deep* is set, the rewriting is applied recursively to the tower
    of extension numbers; otherwise, the rewriting is only applied
    to the top-level extension numbers.

    The result is not a normal form in the strong sense (the same
    number can have many possible representations even after applying
    this transformation), but in practice this is a powerful heuristic
    for simplification.


Factorization
-------------------------------------------------------------------------------

.. type:: ca_factor_struct

.. type:: ca_factor_t

    Represents a real or complex number in factored form
    `b_1^{e_1} b_2^{e_2} \cdots b_n^{e_n}` where `b_i` and `e_i`
    are :type:`ca_t` numbers (the exponents need not be integers).

.. function:: void ca_factor_init(ca_factor_t fac, ca_ctx_t ctx)

    Initializes *fac* and sets it to the empty factorization
    (equivalent to the number 1).

.. function:: void ca_factor_clear(ca_factor_t fac, ca_ctx_t ctx)

    Clears the factorization structure *fac*.

.. function:: void ca_factor_one(ca_factor_t fac, ca_ctx_t ctx)

    Sets *fac* to the empty factorization (equivalent to the number 1).

.. function:: void ca_factor_print(const ca_factor_t fac, ca_ctx_t ctx)

    Prints a description of *fac* to standard output.

.. function:: void ca_factor_insert(ca_factor_t fac, const ca_t base, const ca_t exp, ca_ctx_t ctx)

    Inserts `b^e` into *fac* where *b* is given by *base* and *e*
    is given by *exp*. If a base element structurally identical to *base*
    already exists in *fac*, the corresponding exponent is incremented by *exp*;
    otherwise, this factor is appended.

.. function:: void ca_factor_get_ca(ca_t res, const ca_factor_t fac, ca_ctx_t ctx)

    Expands *fac* back to a single :type:`ca_t` by evaluating the powers
    and multiplying out the result.

.. function:: void ca_factor(ca_factor_t res, const ca_t x, ulong flags, ca_ctx_t ctx)

    Sets *res* to a factorization of *x*
    of the form `x = b_1^{e_1} b_2^{e_2} \cdots b_n^{e_n}`.
    Requires that *x* is not a special value.
    The type of factorization is controlled by *flags*, which can be set
    to a combination of constants in the following section.

Factorization options
................................................................................

The following flags select the structural polynomial factorization to
perform over formal fields `\mathbb{Q}(a_1,\ldots,a_n)`.
Each flag in the list strictly encompasses the factorization power of
the preceding flag, so it is unnecessary to pass more than one flag.

.. macro:: CA_FACTOR_POLY_NONE

    No polynomial factorization at all.

.. macro:: CA_FACTOR_POLY_CONTENT

    Only extract the rational content.

.. macro:: CA_FACTOR_POLY_SQF

    Perform a squarefree factorization in addition to extracting
    the rational content.

.. macro:: CA_FACTOR_POLY_FULL

    Perform a full multivariate polynomial factorization.

The following flags select the factorization to perform over `\mathbb{Z}`.
Integer factorization is applied if *x* is an element of `\mathbb{Q}`, and to
the extracted rational content of polynomials.
Each flag in the list strictly encompasses the factorization power of
the preceding flag, so it is unnecessary to pass more than one flag.

.. macro:: CA_FACTOR_ZZ_NONE

    No integer factorization at all.

.. macro:: CA_FACTOR_ZZ_SMOOTH

    Perform a smooth factorization to extract small prime factors
    (heuristically up to ``CA_OPT_SMOOTH_LIMIT`` bits) in addition to
    identifying perfect powers.

.. macro:: CA_FACTOR_ZZ_FULL

    Perform a complete integer factorization into prime numbers.
    This is prohibitively slow for general integers exceeding 70-80 digits.

.. _context-options:

Context options
-------------------------------------------------------------------------------

The *options* member of a :type:`ca_ctx_t` object is an array of *slong*
values controlling simplification behavior and various other settings.
The values of the array at the following indices can be changed by the user
(example: ``ctx->options[CA_OPT_PREC_LIMIT] = 65536``).

It is recommended to set options controlling evaluation only at the time
when a context object is created. Changing such options later should
normally be harmless, but since the update will not apply
retroactively to objects that have already been computed and cached,
one might not see the expected behavior.
Superficial options (printing) can be changed at any time.

.. macro:: CA_OPT_VERBOSE

    Whether to print debug information. Default value: 0.

.. macro:: CA_OPT_PRINT_FLAGS

    Printing style. See :ref:`ca-printing` for details.
    Default value: ``CA_PRINT_DEFAULT``.

.. macro:: CA_OPT_MPOLY_ORD

    Monomial ordering to use for multivariate polynomials. Possible
    values are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.
    Default value: ``ORD_LEX``.
    This option must be set before doing any computations.

.. macro:: CA_OPT_PREC_LIMIT

    Maximum precision to use internally for numerical evaluation with Arb,
    and in some cases for the magntiude of exact coefficients.
    This parameter affects the possibility to prove inequalities
    and find simplifications between related extension numbers.
    This is not a strict limit; some calculations may use higher precision
    when there is a good reason to do so.
    Default value: 4096.

.. macro:: CA_OPT_QQBAR_DEG_LIMIT

    Maximum degree of :type:`qqbar_t` elements allowed internally during
    simplification of algebraic numbers. This limit may be exceeded
    when the user provides explicit :type:`qqbar_t` input of higher degree.
    Default value: 120.

.. macro:: CA_OPT_LOW_PREC

    Numerical precision to use for fast checks (typically, before attempting
    more expensive operations). Default value: 64.

.. macro:: CA_OPT_SMOOTH_LIMIT

    Size in bits for factors in smooth integer factorization. Default value: 32.

.. macro:: CA_OPT_LLL_PREC

    Precision to use to find integer relations using LLL. Default value: 128.

.. macro:: CA_OPT_POW_LIMIT

    Largest exponent to expand powers automatically. This only applies
    in multivariate and transcendental fields: in number fields,
    ``CA_OPT_PREC_LIMIT`` applies instead. Default value: 20.

.. macro:: CA_OPT_USE_GROEBNER

    Boolean flag for whether to use Gröbner basis computation.
    This flag and the following limits affect the ability to
    prove multivariate identities.
    Default value: 1.

.. macro:: CA_OPT_GROEBNER_LENGTH_LIMIT

    Maximum length of ideal basis allowed in Buchberger's algorithm.
    Default value: 100.

.. macro:: CA_OPT_GROEBNER_POLY_LENGTH_LIMIT

    Maximum length of polynomials allowed in Buchberger's algorithm.
    Default value: 1000.

.. macro:: CA_OPT_GROEBNER_POLY_BITS_LIMIT

    Maximum coefficient size in bits of polynomials allowed in
    Buchberger's algorithm.
    Default value: 10000.

.. macro:: CA_OPT_VIETA_LIMIT

    Maximum degree *n* of algebraic numbers for which to add Vieta's
    formulas to the reduction ideal.
    This must be set relatively low
    since the number of terms in Vieta's formulas is `O(2^n)`
    and the resulting Gröbner basis computations can be expensive.
    Default value: 6.

.. macro:: CA_OPT_TRIG_FORM

    Default representation of trigonometric functions.
    The following values are possible:

    .. macro:: CA_TRIG_DIRECT

        Use the direct functions (with some exceptions).

    .. macro:: CA_TRIG_EXPONENTIAL

        Use complex exponentials.

    .. macro:: CA_TRIG_SINE_COSINE

        Use sines and cosines.

    .. macro:: CA_TRIG_TANGENT

        Use tangents.

    Default value: ``CA_TRIG_EXPONENTIAL``.

    The *exponential* representation is currently used by default
    as typically works best for field arithmetic
    and simplifications, although it has the disadvantage of
    introducing complex numbers where real numbers would be sufficient.
    This may change in the future.



Internal representation
-------------------------------------------------------------------------------

.. macro:: CA_FMPQ(x)

.. macro:: CA_FMPQ_NUMREF(x)

.. macro:: CA_FMPQ_DENREF(x)

    Assuming that *x* holds an element of the trivial field `\mathbb{Q}`,
    this macro returns a pointer which can be used as an :type:`fmpq_t`,
    or respectively to the numerator or denominator as an :type:`fmpz_t`.

.. macro:: CA_MPOLY_Q(x)

    Assuming that *x* holds a generic field element as data,
    this macro returns a pointer which can be used as
    an :type:`fmpz_mpoly_q_t`.

.. macro:: CA_NF_ELEM(x)

    Assuming that *x* holds an Antic number field element as data,
    this macro returns a pointer which can be used as
    an :type:`nf_elem_t`.


.. function:: void _ca_make_field_element(ca_t x, ca_field_srcptr new_index, ca_ctx_t ctx)

    Changes the internal representation of *x* to that of an element
    of the field with index *new_index* in the context object *ctx*.
    This may destroy the value of *x*.

.. function:: void _ca_make_fmpq(ca_t x, ca_ctx_t ctx)

    Changes the internal representation of *x* to that of an element of
    the trivial field `\mathbb{Q}`. This may destroy the value of *x*.


.. raw:: latex

    \newpage

