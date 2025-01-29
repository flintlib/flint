.. _fexpr-builtin:

**fexpr_builtin.h** -- builtin symbols
===============================================================================

This module defines symbol names with a predefined meaning for
use in symbolic expressions. These symbols will eventually all
support LaTeX rendering as well as symbolic and numerical evaluation
(where applicable).

By convention, all builtin symbol names are at least two characters
long and start with an uppercase letter. Single-letter symbol names
and symbol names beginning with a lowercase letter are reserved for
variables.

For any builtin symbol name ``Symbol``, the header file
``fexpr_builtin.h`` defines a C constant ``FEXPR_Symbol`` as an
index to a builtin symbol table.
The symbol will be documented as ``Symbol`` below.

C helper functions
------------------------------------------------------------------------

.. function:: slong fexpr_builtin_lookup(const char * s)

    Returns the internal index used to encode the builtin symbol
    with name *s* in expressions. If *s* is not the name of a builtin
    symbol, returns -1.

.. function:: const char * fexpr_builtin_name(slong n)

    Returns a read-only pointer for a string giving the name of the
    builtin symbol with index *n*.

.. function:: slong fexpr_builtin_length(void)

    Returns the number of builtin symbols.

Variables and iteration
------------------------------------------------------------------------

Expressions involving the following symbols have a special role
in binding variables.

.. macro:: For

    Generator expression. This is a syntactical construct which does
    not represent a mathematical object on its own.
    In general, ``For(x, ...)`` defines the symbol ``x`` as a locally bound
    variable in the scope of the parent expression.
    The following arguments ``...`` specify an evaluation range,
    set or point. Their interpretation depends on the parent
    operator. The following cases are possible.

    Case 1: ``For(x, S)`` specifies iteration or comprehension for ``x``
    ranging over the values of the set ``S``.
    This interpretation is used in operators that aggregate values
    over a set. The ``For`` expression may be followed by a filter
    predicate ``P(x)`` restricting the range to a subset of ``S``.
    Examples:

        ``Set(f(x), For(x, S))`` denotes `\{f(x) : x \in S\}`.

        ``Set(f(x), For(x, S), P(x))`` denotes `\{f(x) : x \in S \operatorname{and} P(x)\}`.

        ``Sum(f(x), For(x, S))`` denotes `\sum_{x \in S} f(x)`.

        ``Sum(f(x), For(x, S), P(x))`` denotes `\sum_{x \in S, \, P(x)} f(x)`.

    Case 2: ``For(x, a, b)`` specifies that ``x`` ranges between
    the endpoints ``a`` and ``b`` in the context of ``Sum``,
    ``Product``, ``Integral``, and similar operators.
    Examples:

        ``Sum(f(n), For(n, a, b))`` denotes `\sum_{n=a}^b f(n)`.
        The iteration is empty if `b < a`.

        ``Integral(f(x), For(x, a, b))`` denotes `\int_a^b f(x) dx`,
        where the integral follows a straight-line path from *a*
        to *b*. Swapping *a* and *b* negates the value.

    Case 3: ``For(x, a)`` specifies that ``x`` approaches the
    point ``a`` in the context of ``Limit``-type operator, or
    differentiation with respect to ``x`` at the point ``a`` in the
    context of a ``Derivative``-type operator. Examples:

        ``Derivative(f(x), For(x, a))`` denotes `f'(a)`.

        ``Limit(f(x), For(x, a))`` denotes `\lim_{x \to a} f(x)`.

    Case 4: ``For(x, a, n)`` specifies differentiation with respect
    to ``x`` at the point ``a`` to order ``n`` in the context of
    a ``Derivative``-type operator. Examples:

        ``Derivative(f(x), For(x, a, n))`` denotes `f^{(n)}(a)`.

.. macro:: Where

    ``Where(f(x), Def(x, a))`` defines the symbol ``x`` as an alias for
    the expression ``a`` and evaluates the expression ``f(x)`` with
    this bound value of ``x``. This is equivalent to ``f(a)``.
    This may be rendered as `f(x) \; \operatorname{where} x = a`.

    ``Where(f(x), Def(f(t), a))`` defines the symbol ``f`` as
    a function mapping the dummy variable ``t`` to ``a``.

    ``Where(Add(a, b), Def(Tuple(a, b), T))`` is a destructuring
    assignment.

.. macro:: Def

    Definition expression. This is a syntactical construct which does
    not represent a mathematical object on its own.
    The ``Def`` expression is used only within a ``Where``-expression;
    see that documentation of that symbol for more examples.

    ``Def(x, a)`` defines the symbol ``x`` as an alias for the
    expression ``a``.

    ``Def(f(x, y, z), a)`` defines the symbol ``f`` as a function
    of three variables. The dummy variables ``x``, ``y`` and ``z``
    may appear within the expression ``a``.

.. macro:: Fun

    ``Fun(x, expr)`` defines an anonymous univariate function mapping
    the symbol ``x`` to the expression ``expr``.
    The symbol ``x`` becomes locally bound within this ``Fun``
    expression.

.. macro:: Step

.. macro:: Repeat

Booleans and logic
------------------------------------------------------------------------

.. macro:: Equal

    ``Equal(a, b)``, signifying `a = b`, is ``True`` if ``a`` and
    ``b`` represent the same object, and ``False`` otherwise.
    This operator can be called with any number of arguments,
    in which case it evaluates whether all arguments are
    equal.

.. macro:: NotEqual

    ``NotEqual(a, b)``, signifying `a \ne b`, is equivalent to
    ``Not(Equal(a, b))``.

.. macro:: Same
 
    ``Same(a, b)`` gives ``a`` (or equivalently ``b``) if ``a`` and
    ``b`` represent the same object, and ``Undefined`` otherwise.
    This can be used to assert or emphasize that two expressions
    represent the same value within a formula.
    This operator can be called with any number of arguments,
    in which case it asserts that all arguments are equal.

.. macro:: True

    ``True`` is a logical constant.

.. macro:: False

    ``False`` is a logical constant.

.. macro:: Not

    ``Not(x)`` is the logical negation of ``x``.

.. macro:: And

    ``And(x, y)`` is the logical AND of ``x`` and ``y``. This function
    can be called with any number of arguments.

.. macro:: Or

    ``Or(x, y)`` is the logical OR of ``x`` and ``y``. This function
    can be called with any number of arguments.

.. macro:: Equivalent

    ``Equivalent(x, y)`` denotes the logical equivalence `x \Leftrightarrow y`.
    Semantically, this is the same as ``Equal``
    called with logical arguments.

.. macro:: Implies

    ``Implies(x, y)`` denotes the logical implication `x \implies y`.

.. macro:: Exists

    Existence quantifier.

    ``Exists(f(x), For(x, S))`` denotes `f(x) \;\text{ for some } x \in S`.

    ``Exists(f(x), For(x, S), P(x))`` denotes `f(x) \;\text{ for some } x \in S \text{ with } P(x)`.

.. macro:: All

    Universal quantifier.

    ``All(f(x), For(x, S))`` denotes `f(x) \;\text{ for all } x \in S`.

    ``All(f(x), For(x, S), P(x))`` denotes `f(x) \;\text{ for all } x \in S \text{ with } P(x)`.


.. macro:: Cases

    ``Cases(Case(f(x), P(x)), Case(g(x), Otherwise))`` denotes:

    .. math::

        \begin{cases} f(x), & P(x)\\g(x), & \text{otherwise}\\ \end{cases}

    ``Cases(Case(f(x), P(x)), Case(g(x), Q(x)), Case(h(x), Otherwise))`` denotes:

    .. math::

        \begin{cases} f(x), & P(x)\\g(x), & Q(x)\\h(x), & \text{otherwise}\\ \end{cases}

    If both `P(x)` and `Q(x)` are true simultaneously, no ordering is implied;
    it is assumed that `f(x)` and `g(x)` give the same value for any such `x`.
    More generally, this operator can be called with any number of case
    distinctions.

    If the *Otherwise* case is omitted, the result is undefined if neither
    predicate is true.

.. macro:: Case

    See ``Cases``.

.. macro:: Otherwise

    See ``Cases``.

Tuples, lists and sets
------------------------------------------------------------------------

.. macro:: Tuple

.. macro:: List

.. macro:: Set

.. macro:: Item

.. macro:: Element

.. macro:: NotElement

.. macro:: EqualAndElement

.. macro:: Length

.. macro:: Cardinality

.. macro:: Concatenation

.. macro:: Union

.. macro:: Intersection

.. macro:: SetMinus

.. macro:: Subset

.. macro:: SubsetEqual

.. macro:: CartesianProduct

.. macro:: CartesianPower

.. macro:: Subsets

    ``Subsets(S)`` is the power set `\mathscr{P}(S)` comprising
    all subsets of the set ``S``.

.. macro:: Sets

    ``Sets`` is the class `\operatorname{Sets}` of all sets.

.. macro:: Tuples

    ``Tuples`` is the class of all tuples.

    ``Tuples(S)`` is the set of all tuples with elements in the
    set ``S``.

    ``Tuples(S, n)`` is the set of all length-``n`` tuples with elements in the
    set ``S``.


Numbers and arithmetic
------------------------------------------------------------------------

Undefined
........................................................................

.. macro:: Undefined

    ``Undefined`` is the special value `\mathfrak{u}` (undefined).

Particular numbers
........................................................................

.. macro:: Pi

    ``Pi`` is the constant `\pi`.

.. macro:: NumberI

    ``NumberI`` is the imaginary unit `i`.
    The verbose name leaves ``i`` and ``I`` to be used as a
    variable names.

.. macro:: NumberE

    ``NumberE`` is the base of the natural logarithm `e`.
    The verbose name leaves ``e`` and ``E`` to be used as a variable
    names.

.. macro:: GoldenRatio

    ``GoldenRatio`` is the golden ratio `\varphi`.

.. macro:: Euler

    ``Euler`` is Euler's constant `\gamma`.

.. macro:: CatalanConstant

    ``CatalanConstant`` is Catalan's constant `G`.

.. macro:: KhinchinConstant

    ``KhinchinConstant`` is Khinchin's constant `K`.

.. macro:: GlaisherConstant

    ``GlaisherConstant`` is Glaisher's constant `A`.

.. macro:: RootOfUnity

    ``RootOfUnity(n)`` is the principal complex *n*-th root of unity `\zeta_n = e^{2 \pi i / n}`.

    ``RootOfUnity(n, k)`` is the complex *n*-th root of unity `\zeta_n^k`.

Number constructors
........................................................................

Remark: the rational number with numerator *p* and denominator *q*
can be constructed as ``Div(p, q)``.

.. macro:: Decimal

    ``Decimal(str)`` gives the rational number specified by the
    string *str* in ordinary decimal floating-point notation
    (for example ``-3.25e-725``).

.. macro:: AlgebraicNumberSerialized

.. macro:: PolynomialRootIndexed

.. macro:: PolynomialRootNearest

.. macro:: Enclosure

.. macro:: Approximation

.. macro:: Guess

.. macro:: Unknown


Arithmetic operations
........................................................................

.. macro:: Pos

.. macro:: Neg

.. macro:: Add

.. macro:: Sub

.. macro:: Mul

.. macro:: Div

.. macro:: Pow

.. macro:: Sqrt

.. macro:: Root


Inequalities
........................................................................

.. macro:: Less

.. macro:: LessEqual

.. macro:: Greater

.. macro:: GreaterEqual

.. macro:: EqualNearestDecimal


Sets of numbers
........................................................................

.. macro:: NN

    ``NN`` is the set of natural numbers (including 0), `\mathbb{N}`.

.. macro:: ZZ

    ``ZZ`` is the set of integers, `\mathbb{Z}`.

.. macro:: QQ

    ``QQ`` is the set of rational numbers, `\mathbb{Q}`.

.. macro:: RR

    ``RR`` is the set of real numbers, `\mathbb{R}`.

.. macro:: CC

    ``CC`` is the set of complex numbers, `\mathbb{C}`.

.. macro:: Primes

    ``Primes`` is the set of positive prime numbers, `\mathbb{P}`

.. macro:: IntegersGreaterEqual

    ``IntegersGreaterEqual(x)``, given an extended real number *x*,
    gives the set `\mathbb{Z}_{\ge x}`
    of integers greater than or equal to *x*.

.. macro:: IntegersLessEqual

    ``IntegersLessEqual(x)``, given an extended real number *x*,
    gives the set `\mathbb{Z}_{\le x}`
    of integers less than or equal to *x*.

.. macro:: Range

    ``Range(a, b)``, given integers *a* and *b*, gives the set
    `\{a, a+1, \ldots, b\}` of integers between *a* and *b*.
    This is the empty set if *a* is greater than *b*.

.. macro:: AlgebraicNumbers

    The set of complex algebraic numbers `\overline{\mathbb{Q}}`.

.. macro:: RealAlgebraicNumbers

    The set of real algebraic numbers `\overline{\mathbb{Q}}_{\mathbb{R}}`.

.. macro:: Interval

    ``Interval(a, b)``, given extended real numbers *a* and *b*, gives
    the closed interval `[a, b]`.

.. macro:: OpenInterval

    ``OpenInterval(a, b)``, given extended real numbers *a* and *b*, gives
    the open interval `(a, b)`.

.. macro:: ClosedOpenInterval

    ``ClosedOpenInterval(a, b)``, given extended real numbers *a* and *b*, gives
    the closed-open interval `[a, b)`.

.. macro:: OpenClosedInterval

    ``OpenClosedInterval(a, b)``, given extended real numbers *a* and *b*, gives
    the closed-open interval `(a, b]`.

.. macro:: RealBall

    ``RealBall(m, r)``, given a real number *m* and an extended real number *r*, gives the
    the closed real ball `[m \pm r]` with center *m* and radius *r*.

.. macro:: OpenRealBall

    ``OpenRealBall(m, r)``, given a real number *m* and an extended real number *r*, gives the
    the open real ball `(m \pm r)` with center *m* and radius *r*.

.. macro:: OpenComplexDisk

    ``OpenComplexDisk(m, r)``, given a complex number *m* and an extended real number *r*,
    gives the open complex disk `D(m, r)` with center *m* and radius *r*.

.. macro:: ClosedComplexDisk

    ``ClosedComplexDisk(m, r)``, given a complex number *m* and a real number *r*,
    gives the closed complex disk `\overline{D}(m, r)` with center *m* and radius *r*.

.. macro:: UpperHalfPlane

    ``UpperHalfPlane`` is the set `\mathbb{H}` of complex numbers
    with positive imaginary part.

.. macro:: UnitCircle

.. macro:: BernsteinEllipse

.. macro:: Lattice


Infinities and extended numbers
........................................................................

.. macro:: Infinity

    ``Infinity`` is the positive signed infinity `\infty`.

.. macro:: UnsignedInfinity

    ``UnsignedInfinity`` is the unsigned infinity `\tilde \infty`.

.. macro:: RealSignedInfinities

    ``RealSignedInfinities`` is the set of real signed infinities
    `\{+\infty, -\infty\}`.

.. macro:: ComplexSignedInfinities

    ``ComplexSignedInfinities`` is the set of complex signed
    infinities `\{e^{i \theta} \cdot \infty : \theta \in \mathbb{R}\}`.

.. macro:: RealInfinities

    ``RealInfinities`` is the set of real infinities (signed
    and unsigned)
    `\{+\infty, -\infty\} \cup \{\tilde \infty\}`.

.. macro:: ComplexInfinities

    ``ComplexInfinities`` is the set of complex infinities (signed
    and unsigned)
    `\{e^{i \theta} \cdot \infty : \theta \in \mathbb{R}\} \cup \{\tilde \infty\}`.

.. macro:: ExtendedRealNumbers

    ``ExtendedRealNumbers`` is the set of extended real numbers
    `\mathbb{R} \cup \{+\infty, -\infty\}`.

.. macro:: ProjectiveRealNumbers

    ``ProjectiveRealNumbers`` is the set of projectively extended real numbers
    `\mathbb{R} \cup \{\tilde \infty\}`.

.. macro:: SignExtendedComplexNumbers

    ``SignExtendedComplexNumbers`` is the set of complex numbers
    extended with signed infinities
    `\mathbb{C} \cup \{e^{i \theta} \cdot \infty : \theta \in \mathbb{R}\}`.

.. macro:: ProjectiveComplexNumbers

    ``ProjectiveComplexNumbers`` is the set of projectively
    extended complex numbers (also known as the Riemann sphere)
    `\mathbb{C} \cup \{\tilde \infty\}`.

.. macro:: RealSingularityClosure

    ``RealSingularityClosure`` is the Calcium singularity closure for real
    functions, encompassing real numbers, signed infinities,
    unsigned infinity, and *undefined* (u). This set is defined as
    `\mathbb{R}_{\text{Sing}} = \mathbb{R} \cup \{+\infty, -\infty\} \cup \{\tilde \infty\} \cup \{ \mathfrak{u} \}`.

.. macro:: ComplexSingularityClosure

    ``ComplexSingularityClosure`` is the Calcium singularity closure for complex
    functions, encompassing complex numbers, signed infinities,
    unsigned infinity, and *undefined* (u). This set is defined as
    `\mathbb{C}_{\text{Sing}} = \mathbb{C} \cup \{e^{i \theta} \cdot \infty : \theta \in \mathbb{R}\} \cup \{\tilde \infty\} \cup \{ \mathfrak{u} \}`.


Operators and calculus
------------------------------------------------------------------------

Sums and products
........................................................................

.. macro:: Sum

.. macro:: Product

.. macro:: PrimeSum

.. macro:: PrimeProduct

.. macro:: DivisorSum

.. macro:: DivisorProduct

Solutions and zeros
........................................................................

.. macro:: Zeros

.. macro:: UniqueZero

.. macro:: Solutions

.. macro:: UniqueSolution

Extreme values
........................................................................

.. macro:: Supremum

.. macro:: Infimum

.. macro:: Minimum

.. macro:: Maximum

.. macro:: ArgMin

.. macro:: ArgMax

.. macro:: ArgMinUnique

.. macro:: ArgMaxUnique

Limits
........................................................................

.. macro:: Limit

.. macro:: SequenceLimit

.. macro:: RealLimit

.. macro:: LeftLimit

.. macro:: RightLimit

.. macro:: ComplexLimit

.. macro:: MeromorphicLimit

.. macro:: SequenceLimitInferior

.. macro:: SequenceLimitSuperior

.. macro:: AsymptoticTo

Derivatives
........................................................................

.. macro:: Derivative

.. macro:: RealDerivative

.. macro:: ComplexDerivative

.. macro:: ComplexBranchDerivative

.. macro:: MeromorphicDerivative

Integrals
........................................................................

.. macro:: Integral

Complex analysis
........................................................................

.. macro:: Path

.. macro:: CurvePath

.. macro:: Poles

.. macro:: IsHolomorphicOn

.. macro:: IsMeromorphicOn

.. macro:: Residue

.. macro:: ComplexZeroMultiplicity

.. macro:: AnalyticContinuation

Matrices and linear algebra
------------------------------------------------------------------------

.. macro:: Matrix

.. macro:: Row

.. macro:: Column

.. macro:: RowMatrix

.. macro:: ColumnMatrix

.. macro:: DiagonalMatrix

.. macro:: Matrix2x2

.. macro:: ZeroMatrix

.. macro:: IdentityMatrix

.. macro:: Det

.. macro:: Spectrum

.. macro:: SingularValues

.. macro:: Matrices

.. macro:: SL2Z

.. macro:: PSL2Z

.. macro:: SpecialLinearGroup

.. macro:: GeneralLinearGroup

.. macro:: HilbertMatrix


Polynomials, series and rings
------------------------------------------------------------------------

.. macro:: Pol

.. macro:: Ser

.. macro:: Polynomial

.. macro:: Coefficient

.. macro:: PolynomialDegree

.. macro:: Polynomials

.. macro:: PolynomialFractions

.. macro:: FormalPowerSeries

.. macro:: FormalLaurentSeries

.. macro:: FormalPuiseuxSeries

.. macro:: Zero

.. macro:: One

.. macro:: Characteristic

.. macro:: Rings

.. macro:: CommutativeRings

.. macro:: Fields

.. macro:: QuotientRing

.. macro:: FiniteField

.. macro:: EqualQSeriesEllipsis

.. macro:: IndefiniteIntegralEqual

.. macro:: QSeriesCoefficient

.. macro:: Call

.. macro:: CallIndeterminate

Special functions
------------------------------------------------------------------------

Number parts and step functions
........................................................................

.. macro:: Abs

.. macro:: Sign

.. macro:: Re

.. macro:: Im

.. macro:: Arg

.. macro:: Conjugate

.. macro:: Csgn

.. macro:: RealAbs

.. macro:: Max

.. macro:: Min

.. macro:: Floor

.. macro:: Ceil

.. macro:: KroneckerDelta

Primes and divisibility
........................................................................

.. macro:: IsOdd

.. macro:: IsEven

.. macro:: CongruentMod

.. macro:: Divides

.. macro:: Mod

.. macro:: GCD

.. macro:: LCM

.. macro:: XGCD

.. macro:: IsPrime

.. macro:: Prime

.. macro:: PrimePi

.. macro:: DivisorSigma

.. macro:: MoebiusMu

.. macro:: EulerPhi

.. macro:: DiscreteLog

.. macro:: LegendreSymbol

.. macro:: JacobiSymbol

.. macro:: KroneckerSymbol

.. macro:: SquaresR

.. macro:: LiouvilleLambda


Elementary functions
........................................................................

.. macro:: Exp

.. macro:: Log

.. macro:: Sin

.. macro:: Cos

.. macro:: Tan

.. macro:: Cot

.. macro:: Sec

.. macro:: Csc

.. macro:: Sinh

.. macro:: Cosh

.. macro:: Tanh

.. macro:: Coth

.. macro:: Sech

.. macro:: Csch

.. macro:: Asin

.. macro:: Acos

.. macro:: Atan

.. macro:: Acot

.. macro:: Asec

.. macro:: Acsc

.. macro:: Asinh

.. macro:: Acosh

.. macro:: Atanh

.. macro:: Acoth

.. macro:: Asech

.. macro:: Acsch

.. macro:: Atan2

.. macro:: Sinc

.. macro:: LambertW


Combinatorial functions
........................................................................

.. macro:: SloaneA

.. macro:: SymmetricPolynomial

.. macro:: Cyclotomic

.. macro:: Fibonacci

.. macro:: BernoulliB

.. macro:: BernoulliPolynomial

.. macro:: StirlingCycle

.. macro:: StirlingS1

.. macro:: StirlingS2

.. macro:: EulerE

.. macro:: EulerPolynomial

.. macro:: BellNumber

.. macro:: PartitionsP

.. macro:: LandauG


Gamma function and factorials
........................................................................

.. macro:: Factorial

.. macro:: Binomial

.. macro:: Gamma

.. macro:: LogGamma

.. macro:: DoubleFactorial

.. macro:: RisingFactorial

.. macro:: FallingFactorial

.. macro:: HarmonicNumber

.. macro:: DigammaFunction

.. macro:: DigammaFunctionZero

.. macro:: BetaFunction

.. macro:: BarnesG

.. macro:: LogBarnesG

.. macro:: StirlingSeriesRemainder

.. macro:: LogBarnesGRemainder

Orthogonal polynomials
........................................................................

.. macro:: ChebyshevT

.. macro:: ChebyshevU

.. macro:: LegendreP

.. macro:: JacobiP

.. macro:: HermiteH

.. macro:: LaguerreL

.. macro:: GegenbauerC

.. macro:: SphericalHarmonicY

.. macro:: LegendrePolynomialZero

.. macro:: GaussLegendreWeight

Exponential integrals
........................................................................

.. macro:: Erf

.. macro:: Erfc

.. macro:: Erfi

.. macro:: UpperGamma

.. macro:: LowerGamma

.. macro:: IncompleteBeta

.. macro:: IncompleteBetaRegularized

.. macro:: LogIntegral

.. macro:: ExpIntegralE

.. macro:: ExpIntegralEi

.. macro:: SinIntegral

.. macro:: SinhIntegral

.. macro:: CosIntegral

.. macro:: CoshIntegral

.. macro:: FresnelC

.. macro:: FresnelS


Bessel and Airy functions
........................................................................

.. macro:: AiryAi

.. macro:: AiryBi

.. macro:: AiryAiZero

.. macro:: AiryBiZero

.. macro:: BesselJ

.. macro:: BesselI

.. macro:: BesselY

.. macro:: BesselK

.. macro:: HankelH1

.. macro:: HankelH2

.. macro:: BesselJZero

.. macro:: BesselYZero

.. macro:: CoulombF

.. macro:: CoulombG

.. macro:: CoulombH

.. macro:: CoulombC

.. macro:: CoulombSigma


Hypergeometric functions
........................................................................

.. macro:: Hypergeometric0F1

.. macro:: Hypergeometric1F1

.. macro:: Hypergeometric1F2

.. macro:: Hypergeometric2F1

.. macro:: Hypergeometric2F2

.. macro:: Hypergeometric2F0

.. macro:: Hypergeometric3F2

.. macro:: HypergeometricU

.. macro:: HypergeometricUStar

.. macro:: HypergeometricUStarRemainder

.. macro:: Hypergeometric0F1Regularized

.. macro:: Hypergeometric1F1Regularized

.. macro:: Hypergeometric1F2Regularized

.. macro:: Hypergeometric2F1Regularized

.. macro:: Hypergeometric2F2Regularized

.. macro:: Hypergeometric3F2Regularized


Zeta and L-functions
........................................................................

.. macro:: RiemannZeta

.. macro:: RiemannZetaZero

.. macro:: RiemannHypothesis

.. macro:: RiemannXi

.. macro:: HurwitzZeta

.. macro:: LerchPhi

.. macro:: PolyLog

.. macro:: MultiZetaValue

.. macro:: DirichletL

.. macro:: DirichletLZero

.. macro:: DirichletLambda

.. macro:: DirichletCharacter

.. macro:: DirichletGroup

.. macro:: PrimitiveDirichletCharacters

.. macro:: GeneralizedRiemannHypothesis

.. macro:: ConreyGenerator

.. macro:: GeneralizedBernoulliB

.. macro:: StieltjesGamma

.. macro:: KeiperLiLambda

.. macro:: GaussSum

Elliptic integrals
........................................................................

.. macro:: AGM

.. macro:: AGMSequence

.. macro:: EllipticK

.. macro:: EllipticE

.. macro:: EllipticPi

.. macro:: IncompleteEllipticF

.. macro:: IncompleteEllipticE

.. macro:: IncompleteEllipticPi

.. macro:: CarlsonRF

.. macro:: CarlsonRG

.. macro:: CarlsonRJ

.. macro:: CarlsonRD

.. macro:: CarlsonRC

.. macro:: CarlsonHypergeometricR

.. macro:: CarlsonHypergeometricT


Elliptic, theta and modular functions
........................................................................

.. macro:: JacobiTheta

.. macro:: JacobiThetaQ

.. macro:: DedekindEta

.. macro:: ModularJ

.. macro:: ModularLambda

.. macro:: EisensteinG

.. macro:: EisensteinE

.. macro:: DedekindSum

.. macro:: WeierstrassP

.. macro:: WeierstrassZeta

.. macro:: WeierstrassSigma

.. macro:: EllipticRootE

.. macro:: HilbertClassPolynomial

.. macro:: EulerQSeries

.. macro:: DedekindEtaEpsilon

.. macro:: ModularGroupAction

.. macro:: ModularGroupFundamentalDomain

.. macro:: ModularLambdaFundamentalDomain

.. macro:: PrimitiveReducedPositiveIntegralBinaryQuadraticForms

.. macro:: JacobiThetaEpsilon

.. macro:: JacobiThetaPermutation

Nonsemantic markup
........................................................................

.. macro:: Ellipsis

    ``Ellipsis`` renders as `\ldots` in LaTeX. It can be used to
    indicate missing function arguments for display purposes,
    but it has no predefined builtin semantics.

.. macro:: Parentheses

    ``Parentheses(x)`` semantically represents ``x``, but renders
    with parentheses (`\left(x\right)`) when converted to LaTeX.

.. macro:: Brackets

    ``Brackets(x)`` semantically represents ``x``, but renders
    with brackets (`\left[x\right]`) when converted to LaTeX.

.. macro:: Braces

    ``Braces(x)`` semantically represents ``x``, but renders
    with braces (`\left\{x\right\}`) when converted to LaTeX.

.. macro:: AngleBrackets

    ``AngleBrackets(x)`` semantically represents ``x``, but renders
    with angle brackets (`\left\langle x\right\rangle`) when
    converted to LaTeX.

.. macro:: Logic

    ``Logic(x)`` semantically represents ``x``, but forces logical
    expressions within *x* to be rendered using symbols instead
    of text.

.. macro:: ShowExpandedNormalForm

    ``ShowExpandedNormalForm(x)`` semantically represents ``x``, but
    displays the expanded normal form of the expression instead of
    rendering the expression verbatim.
    Warning: this triggers a nontrivial (potentially very expensive)
    computation.

.. macro:: Subscript


.. raw:: latex

    \newpage
