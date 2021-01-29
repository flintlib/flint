latex_test_cases = [
    ("""f(0)""", "f(0)"),
    ("""f("Hello, world!")""", "f(\\text{``Hello, world!''})"),
    ("f", r"f"),
    ("f_", r"f"),
    ("f_(0)", r"f_{0}"),
    ("f()", r"f()"),
    ("Add(Add(Add(f(a, b), c_(n)), f_(x, y)), f_())", r"f(a, b) + c_{n} + f_{x, y} + f_{}"),
    ("f(alpha, beta, chi, delta, ell, epsilon, eta)", r"f(\alpha, \beta, \chi, \delta, \ell, \varepsilon, \eta)"),
    ("f(gamma, iota, kappa, lamda, mu, nu, omega, phi)", r"f(\gamma, \iota, \kappa, \lambda, \mu, \nu, \omega, \phi)"),
    ("f(pi, rho, sigma, tau, theta, varphi, vartheta, xi, zeta)", r"f(\pi, \rho, \sigma, \tau, \theta, \varphi, \vartheta, \xi, \zeta)"),
    ("f(Delta, GreekGamma, GreekPi, Lamda, Omega, Phi, Psi, Sigma, Theta, Xi)", r"f(\Delta, \Gamma, \Pi, \Lambda, \Omega, \Phi, \Psi, \Sigma, \Theta, \Xi)"),
    ("f(alpha_, beta_, chi_, delta_, ell_, epsilon_, eta_)", r"f(\alpha, \beta, \chi, \delta, \ell, \varepsilon, \eta)"),
    ("f(gamma_, iota_, kappa_, lamda_, mu_, nu_, omega_, phi_)", r"f(\gamma, \iota, \kappa, \lambda, \mu, \nu, \omega, \phi)"),
    ("f(pi_, rho_, sigma_, tau_, theta_, varphi_, vartheta_, xi_, zeta)", r"f(\pi, \rho, \sigma, \tau, \theta, \varphi, \vartheta, \xi, \zeta)"),
    ("f(Delta_, GreekGamma_, GreekPi_, Lamda_, Omega_, Phi_, Psi_, Sigma_, Theta_, Xi_)", r"f(\Delta, \Gamma, \Pi, \Lambda, \Omega, \Phi, \Psi, \Sigma, \Theta, \Xi)"),
    ("f(alpha(x), beta(x), chi(x), delta(x), ell(x), epsilon(x), eta(x))", r"f\!\left(\alpha(x), \beta(x), \chi(x), \delta(x), \ell(x), \varepsilon(x), \eta(x)\right)"),
    ("f(gamma(x), iota(x), kappa(x), lamda(x), mu(x), nu(x), omega(x), phi(x))", r"f\!\left(\gamma(x), \iota(x), \kappa(x), \lambda(x), \mu(x), \nu(x), \omega(x), \phi(x)\right)"),
    ("f(pi(x), rho(x), sigma(x), tau(x), theta(x), varphi(x), vartheta(x), xi(x), zeta(x))", r"f\!\left(\pi(x), \rho(x), \sigma(x), \tau(x), \theta(x), \varphi(x), \vartheta(x), \xi(x), \zeta(x)\right)"),
    ("f(Delta(x), GreekGamma(x), GreekPi(x), Lamda(x), Omega(x), Phi(x), Psi(x), Sigma(x), Theta(x), Xi(x))", r"f\!\left(\Delta(x), \Gamma(x), \Pi(x), \Lambda(x), \Omega(x), \Phi(x), \Psi(x), \Sigma(x), \Theta(x), \Xi(x)\right)"),
    ("f(alpha_(n), beta_(n), chi_(n), delta_(n), ell_(n), epsilon_(n), eta_(n))", r"f\!\left(\alpha_{n}, \beta_{n}, \chi_{n}, \delta_{n}, \ell_{n}, \varepsilon_{n}, \eta_{n}\right)"),
    ("f(gamma_(n), iota_(n), kappa_(n), lamda_(n), mu_(n), nu_(n), omega_(n), phi_(n))", r"f\!\left(\gamma_{n}, \iota_{n}, \kappa_{n}, \lambda_{n}, \mu_{n}, \nu_{n}, \omega_{n}, \phi_{n}\right)"),
    ("f(pi_(n), rho_(n), sigma_(n), tau_(n), theta_(n), varphi_(n), vartheta_(n), xi_(n), zeta_(n))", r"f\!\left(\pi_{n}, \rho_{n}, \sigma_{n}, \tau_{n}, \theta_{n}, \varphi_{n}, \vartheta_{n}, \xi_{n}, \zeta_{n}\right)"),
    ("f(Delta_(n), GreekGamma_(n), GreekPi_(n), Lamda_(n), Omega_(n), Phi_(n), Psi_(n), Sigma_(n), Theta_(n), Xi_(n))", r"f\!\left(\Delta_{n}, \Gamma_{n}, \Pi_{n}, \Lambda_{n}, \Omega_{n}, \Phi_{n}, \Psi_{n}, \Sigma_{n}, \Theta_{n}, \Xi_{n}\right)"),
    ("f(a)(b)", r"f(a)(b)"),
    ("Mul(f, g)(x)", r"\left(f g\right)(x)"),
    ("Add(f, g)(x)", r"\left(f + g\right)(x)"),
    ("c_(m, n, p)(x, y, z)", r"c_{m, n, p}(x, y, z)"),
    ("f_(Div(-3, 2))(Div(-3, 2))", r"f_{-3 / 2}\!\left(-\frac{3}{2}\right)"),
    ("Mul(2, Pi, NumberI)", r"2 \pi i"),
    ("Mul(-2, Pi, NumberI)", r"-2 \pi i"),
    ("Add(2, x)", r"2 + x"),
    ("Sub(2, x)", r"2 - x"),
    ("Add(2, Neg(x))", r"2 + \left(-x\right)"),
    ("Sub(2, Neg(x))", r"2 - \left(-x\right)"),
    ("Add(2, Sub(x, y))", r"2 + \left(x - y\right)"),
    ("Sub(2, Add(x, y))", r"2 - \left(x + y\right)"),
    ("Sub(2, Sub(x, y))", r"2 - \left(x - y\right)"),
    ("Add(-3, -4, -5)", r"-3-4-5"),
    ("Add(-3, Mul(-4, x), -5)", r"-3-4 x-5"),
    ("Add(-3, Mul(-4, x), Pos(5))", r"-3-4 x+5"),
    ("Add(-3, Div(Mul(-4, x), 7), -5)", r"-3-\frac{4 x}{7}-5"),
    ("Add(-3, Div(Mul(4, x), 7), -5)", r"-3 + \frac{4 x}{7}-5"),
    ("Mul(0, 1)", r"0 \cdot 1"),
    ("Mul(3, Pow(2, n))", r"3 \cdot {2}^{n}"),
    ("Mul(3, Pow(-1, n))", r"3 \cdot {\left(-1\right)}^{n}"),
    ("Mul(-3, Pow(-1, n))", r"-3 \cdot {\left(-1\right)}^{n}"),
    ("Mul(-1, -2, -3)", r"-1 \cdot \left(-2\right) \cdot \left(-3\right)"),
    ("Div(-1, 3)", r"-\frac{1}{3}"),
    ("Div(Mul(-5, Pi), 3)", r"-\frac{5 \pi}{3}"),
    ("Div(Neg(Mul(5, Pi)), 3)", r"-\frac{5 \pi}{3}"),
    ("Div(Add(Add(Mul(-5, Pow(x, 2)), Mul(4, x)), 1), Add(Mul(3, x), y))", r"\frac{-5 {x}^{2} + 4 x + 1}{3 x + y}"),
    ("Pow(2, n)", r"{2}^{n}"),
    ("Pow(-1, n)", r"{\left(-1\right)}^{n}"),
    ("Pow(10, Pow(10, -10))", r"{10}^{{10}^{-10}}"),
    ("Pow(2, Div(-1, 3))", r"{2}^{-1 / 3}"),
    ("Pow(2, Div(-1, Mul(3, n)))", r"{2}^{-1 / \left(3 n\right)}"),
    ("Equal(Add(Pow(Sin(x), 2), Pow(Cos(x), 2)), 1)", r"\sin^{2}\!\left(x\right) + \cos^{2}\!\left(x\right) = 1"),
    ("Set(Set(), Set(1), Set(1, 2, 3))", r"\left\{\left\{\right\}, \left\{1\right\}, \left\{1, 2, 3\right\}\right\}"),
    ("Tuple(Tuple(), Tuple(1), Tuple(1, 2, 3))", r"\left(\left(\right), \left(1\right), \left(1, 2, 3\right)\right)"),
    ("List(List(), List(1), List(1, 2, 3))", r"\left[\left[\right], \left[1\right], \left[1, 2, 3\right]\right]"),
    ("Set(f(x), For(x, CC))", r"\left\{ f(x) : x \in \mathbb{C} \right\}"),
    ("Set(f(x), For(x, CC), NotEqual(x, 0))", r"\left\{ f(x) : x \in \mathbb{C}\,\mathbin{\operatorname{and}}\, x \ne 0 \right\}"),
    ("List(Floor(x), Ceil(x), Abs(x), RealAbs(x), Conjugate(z), Sqrt(x))", r"\left[\left\lfloor x \right\rfloor, \left\lceil x\right\rceil, \left|x\right|, \left|x\right|, \overline{z}, \sqrt{x}\right]"),
    ("List(Floor(Div(1, 2)), Ceil(Div(1, 2)), Abs(Div(1, 2)), RealAbs(Div(1, 2)), Conjugate(Div(1, 2)), Sqrt(Div(1, 2)))", r"\left[\left\lfloor \frac{1}{2} \right\rfloor, \left\lceil \frac{1}{2}\right\rceil, \left|\frac{1}{2}\right|, \left|\frac{1}{2}\right|, \overline{\frac{1}{2}}, \sqrt{\frac{1}{2}}\right]"),
    ("And(Equal(Length(Tuple(1, 2, 3)), 3), Equal(Cardinality(Set()), 0))", r"\# \left(1, 2, 3\right) = 3 \;\mathbin{\operatorname{and}}\; \# \left\{\right\} = 0"),
    ("Tuple(Parentheses(x), Brackets(x), Braces(x), AngleBrackets(x))", r"\left(\left(x\right), \left[x\right], \left\{x\right\}, \left\langle x\right\rangle\right)"),
    ("Tuple(Parentheses(Div(1, 2)), Brackets(Div(1, 2)), Braces(Div(1, 2)), AngleBrackets(Div(1, 2)))", r"\left(\left(\frac{1}{2}\right), \left[\frac{1}{2}\right], \left\{\frac{1}{2}\right\}, \left\langle \frac{1}{2}\right\rangle\right)"),
    ("Concatenation(A, B)", r"A  \,^\frown  B"),
    ("Equal(Concatenation(Tuple(a, b), Tuple(c, d, e), Tuple()), Tuple(a, b, c, d, e))", r"\left(a, b\right)  \,^\frown  \left(c, d, e\right)  \,^\frown  \left(\right) = \left(a, b, c, d, e\right)"),
    ("Equal(f(x), Cases(Case(y, P(x)), Case(Neg(y), Q(x))))", r"f(x) = \begin{cases} y, & P(x)\\-y, & Q(x)\\ \end{cases}"),
    ("Equal(f(x), Cases(Case(y, P(x)), Case(Neg(y), Q(x)), Case(0, Otherwise)))", r"f(x) = \begin{cases} y, & P(x)\\-y, & Q(x)\\0, & \text{otherwise}\\ \end{cases}"),
    ("And(Equal(True, Not(False)), Equal(False, Not(True)))", r"\operatorname{True} = \operatorname{not} \operatorname{False} \;\mathbin{\operatorname{and}}\; \operatorname{False} = \operatorname{not} \operatorname{True}"),
    ("All(Greater(x, 0), For(x, S))", r"x > 0 \;\text{ for all } x \in S"),
    ("All(Greater(x, 0), For(x, S), P(x))", r"x > 0 \;\text{ for all } x \in S \text{ with } P(x)"),
    ("Exists(Greater(x, 0), For(x, S))", r"x > 0 \;\text{ for some } x \in S"),
    ("Exists(Greater(x, 0), For(x, S), P(x))", r"x > 0 \;\text{ for some } x \in S \text{ with } P(x)"),
    ("Logic(All(Greater(x, 0), For(x, S)))", r"\forall x \in S : \, x > 0"),
    ("Logic(All(Greater(x, 0), For(x, S), P(x)))", r"\forall x \in S, \,P(x) : \, x > 0"),
    ("Logic(Exists(Greater(x, 0), For(x, S)))", r"\exists x \in S : \, x > 0"),
    ("Logic(Exists(Greater(x, 0), For(x, S), P(x)))", r"\exists x \in S, \,P(x) : \, x > 0"),
    ("Or(Q, And(P, Q, Not(P), Or(Q, P), Not(Or(Q, P))))", r"Q \;\mathbin{\operatorname{or}}\; \left(P \;\mathbin{\operatorname{and}}\; Q \;\mathbin{\operatorname{and}}\; \operatorname{not} P \;\mathbin{\operatorname{and}}\; \left(Q \;\mathbin{\operatorname{or}}\; P\right) \;\mathbin{\operatorname{and}}\; \operatorname{not} \,\left(Q \;\mathbin{\operatorname{or}}\; P\right)\right)"),
    ("Logic(Or(Q, And(P, Q, Not(P), Or(Q, P), Not(Or(Q, P)))))", r"Q \,\lor\, \left(P \,\land\, Q \,\land\, \neg P \,\land\, \left(Q \,\lor\, P\right) \,\land\, \neg \left(Q \,\lor\, P\right)\right)"),
    ("Equivalent(A, B)", r"A \iff B"),
    ("Equivalent(Not(Equal(x, y)), NotEqual(x, y))", r"\left(\operatorname{not} \,\left(x = y\right)\right) \iff \left(x \ne y\right)"),
    ("Implies(P, And(R, S))", r"P \;\implies\; \left(R \;\mathbin{\operatorname{and}}\; S\right)"),
    ("Implies(Element(x, QQ), Element(x, RR))", r"x \in \mathbb{Q} \;\implies\; x \in \mathbb{R}"),
    ("And(Less(x, y), Less(x, y, z), LessEqual(x, y), LessEqual(x, y, z))", r"x < y \;\mathbin{\operatorname{and}}\; x < y < z \;\mathbin{\operatorname{and}}\; x \le y \;\mathbin{\operatorname{and}}\; x \le y \le z"),
    ("And(Greater(x, y), Greater(x, y, z), GreaterEqual(x, y), GreaterEqual(x, y, z))", r"x > y \;\mathbin{\operatorname{and}}\; x > y > z \;\mathbin{\operatorname{and}}\; x \ge y \;\mathbin{\operatorname{and}}\; x \ge y \ge z"),
    ("Subset(Primes, NN, ZZ, QQ, RR, CC)", r"\mathbb{P} \subset \mathbb{N} \subset \mathbb{Z} \subset \mathbb{Q} \subset \mathbb{R} \subset \mathbb{C}"),
    ("Subset(QQ, AlgebraicNumbers, CC)", r"\mathbb{Q} \subset \overline{\mathbb{Q}} \subset \mathbb{C}"),
    ("SubsetEqual(S, QQ)", r"S \subseteq \mathbb{Q}"),
    ("NotElement(123456789012345678901234567890, SetMinus(QQ, ZZ))", r"123456789012345678901234567890 \notin \mathbb{Q} \setminus \mathbb{Z}"),
    ("KroneckerDelta(x, Div(1, 2))", r"\delta_{(x,1 / 2)}"),
    ("Set(Interval(a, b), OpenInterval(a, b), ClosedOpenInterval(a, b), OpenClosedInterval(a, b))", r"\left\{\left[a, b\right], \left(a, b\right), \left[a, b\right), \left(a, b\right]\right\}"),
    ("Set(Interval(a, Div(1, 2)), OpenInterval(a, Div(1, 2)), ClosedOpenInterval(a, Div(1, 2)), OpenClosedInterval(a, Div(1, 2)))", r"\left\{\left[a, 1 / 2\right], \left(a, 1 / 2\right), \left[a, 1 / 2\right), \left(a, 1 / 2\right]\right\}"),
    ("Set(RealBall(m, r), OpenRealBall(m, r))", r"\left\{\left[m \pm r\right], \left(m \pm r\right)\right\}"),
    ("Set(ClosedComplexDisk(m, r), OpenComplexDisk(m, r))", r"\left\{\overline{D}(m, r), D(m, r)\right\}"),
    ("Set(Undefined, UnsignedInfinity, Pos(Infinity), Neg(Infinity))", r"\left\{\mathfrak{u}, \hat{\infty}, +\infty, -\infty\right\}"),
    ("Equal(RealSignedInfinities, Set(Pos(Infinity), Neg(Infinity)))", r"\{\pm \infty\} = \left\{+\infty, -\infty\right\}"),
    ("Equal(ComplexSignedInfinities, Set(Mul(Exp(Mul(NumberI, theta)), Infinity), For(theta, OpenClosedInterval(Neg(Pi), Pi))))", r"\{[e^{i \theta}] \infty\} = \left\{ e^{i \theta} \infty : \theta \in \left(-\pi, \pi\right] \right\}"),
    ("Equal(RealInfinities, Union(RealSignedInfinities, Set(UnsignedInfinity)))", r"\{\hat{\infty}, \pm \infty\} = \{\pm \infty\} \cup \left\{\hat{\infty}\right\}"),
    ("Equal(ComplexInfinities, Union(ComplexSignedInfinities, Set(UnsignedInfinity)))", r"\{\hat{\infty}, [e^{i \theta}] \infty\} = \{[e^{i \theta}] \infty\} \cup \left\{\hat{\infty}\right\}"),
    ("Equal(ExtendedRealNumbers, Union(RR, RealSignedInfinities))", r"\overline{\mathbb{R}} = \mathbb{R} \cup \{\pm \infty\}"),
    ("Equal(SignExtendedComplexNumbers, Union(CC, ComplexSignedInfinities))", r"\overline{\mathbb{C}}_{[e^{i \theta}] \infty} = \mathbb{C} \cup \{[e^{i \theta}] \infty\}"),
    ("Equal(ProjectiveRealNumbers, Union(RR, Set(UnsignedInfinity)))", r"\hat{\mathbb{R}}_{\infty} = \mathbb{R} \cup \left\{\hat{\infty}\right\}"),
    ("Equal(ProjectiveComplexNumbers, Union(CC, Set(UnsignedInfinity)))", r"\hat{\mathbb{C}}_{\infty} = \mathbb{C} \cup \left\{\hat{\infty}\right\}"),
    ("Equal(RealSingularityClosure, Union(RR, Set(UnsignedInfinity), RealSignedInfinities, Set(Undefined)))", r"\overline{\mathbb{R}}_{\text{Sing}} = \mathbb{R} \cup \left\{\hat{\infty}\right\} \cup \{\pm \infty\} \cup \left\{\mathfrak{u}\right\}"),
    ("Equal(ComplexSingularityClosure, Union(CC, Set(UnsignedInfinity), ComplexSignedInfinities, Set(Undefined)))", r"\overline{\mathbb{C}}_{\text{Sing}} = \mathbb{C} \cup \left\{\hat{\infty}\right\} \cup \{[e^{i \theta}] \infty\} \cup \left\{\mathfrak{u}\right\}"),
    ("ArgMin(Add(f(x), g(x)), For(x, RR), NotEqual(x, 0))", r"\mathop{\operatorname{arg\,min}\,}\limits_{x \in \mathbb{R},\,x \ne 0} \left[f(x) + g(x)\right]"),
    ("List(ArgMin(f(x), For(x, S)), ArgMax(f(x), For(x, S)), ArgMin(f(x), For(x, S), P(x)), ArgMax(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{arg\,min}\,}\limits_{x \in S} f(x), \mathop{\operatorname{arg\,max}\,}\limits_{x \in S} f(x), \mathop{\operatorname{arg\,min}\,}\limits_{x \in S,\,P(x)} f(x), \mathop{\operatorname{arg\,max}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
    ("List(Minimum(f(x), For(x, S)), Maximum(f(x), For(x, S)), Minimum(f(x), For(x, S), P(x)), Maximum(f(x), For(x, S), P(x)))", r"\left[\mathop{\min\,}\limits_{x \in S} f(x), \mathop{\max\,}\limits_{x \in S} f(x), \mathop{\min\,}\limits_{x \in S,\,P(x)} f(x), \mathop{\max\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
    ("List(ArgMinUnique(f(x), For(x, S)), ArgMaxUnique(f(x), For(x, S)), ArgMinUnique(f(x), For(x, S), P(x)), ArgMaxUnique(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{arg\,min*}\,}\limits_{x \in S} f(x), \mathop{\operatorname{arg\,max*}\,}\limits_{x \in S} f(x), \mathop{\operatorname{arg\,min*}\,}\limits_{x \in S,\,P(x)} f(x), \mathop{\operatorname{arg\,max*}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
    ("List(Infimum(f(x), For(x, S)), Supremum(f(x), For(x, S)), Infimum(f(x), For(x, S), P(x)), Supremum(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{inf}\,}\limits_{x \in S} f(x), \mathop{\operatorname{sup}\,}\limits_{x \in S} f(x), \mathop{\operatorname{inf}\,}\limits_{x \in S,\,P(x)} f(x), \mathop{\operatorname{sup}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
    ("List(Solutions(Q(x), For(x, S)), Zeros(f(x), For(x, S)), Solutions(Q(x), For(x, S), P(x)), Zeros(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{solutions}\,}\limits_{x \in S} Q(x), \mathop{\operatorname{zeros}\,}\limits_{x \in S} f(x), \mathop{\operatorname{solutions}\,}\limits_{x \in S,\,P(x)} Q(x), \mathop{\operatorname{zeros}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
    ("List(UniqueSolution(Q(x), For(x, S)), UniqueZero(f(x), For(x, S)), UniqueSolution(Q(x), For(x, S), P(x)), UniqueZero(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{solution*}\,}\limits_{x \in S} Q(x), \mathop{\operatorname{zero*}\,}\limits_{x \in S} f(x), \mathop{\operatorname{solution*}\,}\limits_{x \in S,\,P(x)} Q(x), \mathop{\operatorname{zero*}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
    ("Sum(f(n) + g(n), For(n, a, b))", r"\sum_{n=a}^{b} \left(f(n) + g(n)\right)"),
    ("Sum(f(n), For(n, ZZ))", r"\sum_{n  \in \mathbb{Z}} f(n)"),
    ("Sum(f(n), For(n, ZZ), NotEqual(n, 0))", r"\sum_{\textstyle{n  \in \mathbb{Z} \atop n \ne 0}} f(n)"),
    ("Sum(f(n), For(n, a, b), NotEqual(n, 0))", r"\sum_{\textstyle{n=a \atop n \ne 0}}^{b} f(n)"),
    ("Sum(f(n), For(n, a, b))", r"\sum_{n=a}^{b} f(n)"),
    ("Product(f(n) + g(n), For(n, a, b))", r"\prod_{n=a}^{b} \left(f(n) + g(n)\right)"),
    ("Product(f(n), For(n, NN))", r"\prod_{n  \in \mathbb{N}} f(n)"),
    ("Product(f(n), For(n, NN), NotEqual(g(n), 0))", r"\prod_{\textstyle{n  \in \mathbb{N} \atop g(n) \ne 0}} f(n)"),
    ("Product(f(n), For(n, a, b), NotEqual(n, 0))", r"\prod_{\textstyle{n=a \atop n \ne 0}}^{b} f(n)"),
    ("Product(f(n), For(n, a, b))", r"\prod_{n=a}^{b} f(n)"),
    ("Equal(Set(f(n), For(n, ZZ)), Union(Set(f(n), For(n, ZZ), IsEven(n)), Set(f(n), For(n, ZZ), IsOdd(n))))", r"\left\{ f(n) : n \in \mathbb{Z} \right\} = \left\{ f(n) : n \in \mathbb{Z}\,\mathbin{\operatorname{and}}\, n \text{ even} \right\} \cup \left\{ f(n) : n \in \mathbb{Z}\,\mathbin{\operatorname{and}}\, n \text{ odd} \right\}"),
    ("Equal(Primes, Set(p, For(p, NN), IsPrime(p)))", r"\mathbb{P} = \left\{ p : p \in \mathbb{N}\,\mathbin{\operatorname{and}}\, p \text{ prime} \right\}"),
    ("Equal(Sum(f(n), Element(n, ZZ)), Add(Sum(f(n), Element(n, ZZ), IsOdd(n)), Sum(f(n), Element(n, ZZ), IsEven(n))))", r"\sum_{n  \in \mathbb{Z}} f(n) = \sum_{\textstyle{n  \in \mathbb{Z} \atop n \text{ odd}}} f(n) + \sum_{\textstyle{n  \in \mathbb{Z} \atop n \text{ even}}} f(n)"),
    ("Integral(f(x), For(x, -Infinity, Infinity))", r"\int_{-\infty}^{\infty} f(x) \, dx"),
    ("Integral(f(x), For(x, RR))", r"\int_{x \in \mathbb{R}} f(x) \, dx"),
    ("Integral(f(x) + g(x) / h(x), For(x, a, b))", r"\int_{a}^{b} \left(f(x) + \frac{g(x)}{h(x)}\right) \, dx"),
    ("Set(Derivative(f(x_), For(x_, x)), Derivative(f(x_), For(x_, Div(x, y))), Derivative(Gamma(x_), For(x_, 1)))", r"\left\{f'\!\left(x\right), f'\!\left(\frac{x}{y}\right), \Gamma'\!\left(1\right)\right\}"),
    ("Set(Derivative(f(x_), For(x_, x, 0)), Derivative(f(x_), For(x_, x, 1)), Derivative(f(x_), For(x_, x, 2)), Derivative(f(x_), For(x_, x, 3)), Derivative(f(x_), For(x_, x, 4)), Derivative(f(x_), For(x_, x, n)), Derivative(f(x_), For(x_, x, Add(Mul(2, n), 3))))", r"\left\{{f}^{(0)}\!\left(x\right), f'\!\left(x\right), f''\!\left(x\right), f'''\!\left(x\right), {f}^{(4)}\!\left(x\right), {f}^{(n)}\!\left(x\right), {f}^{(2 n + 3)}\!\left(x\right)\right\}"),
    ("Set(Derivative(f(Add(x, 1)), For(x, x)), Derivative(f(Add(x, 1)), For(x, x, 0)), Derivative(f(Add(x, 1)), For(x, x, 1)), Derivative(f(Add(x, 1)), For(x, x, n)))", r"\left\{\frac{d}{d x}\, f\!\left(x + 1\right), \frac{d^{0}}{{d x}^{0}}\, f\!\left(x + 1\right), \frac{d}{d x}\, f\!\left(x + 1\right), \frac{d^{n}}{{d x}^{n}}\, f\!\left(x + 1\right)\right\}"),
    ("Set(Derivative(Add(f(x), g(x)), For(x, Add(y, 3))), Derivative(Add(f(x), g(x)), For(x, Add(y, 3), 5)))", r"\left\{\left[\frac{d}{d x}\, \left[f(x) + g(x)\right] \right]_{x = y + 3}, \left[\frac{d^{5}}{{d x}^{5}}\, \left[f(x) + g(x)\right] \right]_{x = y + 3}\right\}"),
    ("Set(RealDerivative(f(x), For(x, 1)), ComplexDerivative(f(x), For(x, 1)), ComplexBranchDerivative(f(x), For(x, 1)), MeromorphicDerivative(f(x), For(x, 1)))", r"\left\{f'\!\left(1\right), f'\!\left(1\right), f'\!\left(1\right), f'\!\left(1\right)\right\}"),
    ("Set(Limit(f(x), For(x, a)), Limit(f(x), For(x, a), P(x)))", r"\left\{\lim_{x \to a} f(x), \lim_{x \to a,\,P(x)} f(x)\right\}"),
    ("Set(Limit(f(x), For(x, a)), RealLimit(f(x), For(x, a)), ComplexLimit(f(x), For(x, a)), MeromorphicLimit(f(x), For(x, a)))", r"\left\{\lim_{x \to a} f(x), \lim_{x \to a} f(x), \lim_{x \to a} f(x), \lim_{x \to a} f(x)\right\}"),
    ("Set(LeftLimit(f(x), For(x, 0)), RightLimit(f(x), For(x, 0)))", r"\left\{\lim_{x \to {0}^{-}} f(x), \lim_{x \to {0}^{+}} f(x)\right\}"),
    ("Set(SequenceLimit(f(n), For(n, Infinity)), SequenceLimitInferior(f(n), For(n, Infinity)), SequenceLimitSuperior(f(n), For(n, Infinity)))", r"\left\{\lim_{n \to \infty} f(n), \liminf_{n \to \infty} f(n), \limsup_{n \to \infty} f(n)\right\}"),
    ("Sub(Limit(Add(f(x), g(x)), For(x, a)), Limit(Sub(f(x), g(x)), For(x, a)))", r"\lim_{x \to a} \left[f(x) + g(x)\right] - \lim_{x \to a} \left[f(x) - g(x)\right]"),
    ("Divides(GCD(a, b), LCM(a, b))", r"\gcd(a, b) \mid \operatorname{lcm}(a, b)"),
    ("Set(Exp(x), Exp(Div(3, 2)), Exp(Add(Neg(Pow(x, 2)), x)), Exp(Abs(Im(z))), Exp(Div(3, Add(2, x))), Exp(Sin(x)))", r"\left\{e^{x}, e^{3 / 2}, e^{-{x}^{2} + x}, e^{\left|\operatorname{Im}(z)\right|}, \exp\!\left(\frac{3}{2 + x}\right), \exp\!\left(\sin(x)\right)\right\}"),
    ("Add(Sin(x), Cos(x), Tan(x), Cot(x), Sec(x), Csc(x))", r"\sin(x) + \cos(x) + \tan(x) + \cot(x) + \sec(x) + \csc(x)"),
    ("Add(Sinh(x), Cosh(x), Tanh(x), Coth(x), Sech(x), Csch(x))", r"\sinh(x) + \cosh(x) + \tanh(x) + \coth(x) + \operatorname{sech}(x) + \operatorname{csch}(x)"),
    ("Add(Asin(x), Acos(x), Atan(x), Acot(x), Asec(x), Acsc(x))", r"\operatorname{asin}(x) + \operatorname{acos}(x) + \operatorname{atan}(x) + \operatorname{acot}(x) + \operatorname{asec}(x) + \operatorname{acsc}(x)"),
    ("Add(Asinh(x), Acosh(x), Atanh(x), Acoth(x), Asech(x), Acsch(x))", r"\operatorname{asinh}(x) + \operatorname{acosh}(x) + \operatorname{atanh}(x) + \operatorname{acoth}(x) + \operatorname{asech}(x) + \operatorname{acsch}(x)"),
    ("Exp(Neg(Euler))", r"e^{-\gamma}"),
    ("Set(Re(z), Im(z), Atan2(y, x))", r"\left\{\operatorname{Re}(z), \operatorname{Im}(z), \operatorname{atan2}(y, x)\right\}"),
    ("Add(NumberE, GoldenRatio, CatalanConstant)", r"e + \varphi + G"),
    ("Add(Sinc(x), Pow(Sinc(x), 2))", r"\operatorname{sinc}(x) + \operatorname{sinc}^{2}\!\left(x\right)"),
    ("AGM(a, b)", r"\operatorname{agm}(a, b)"),
    ("And(Equal(LogBarnesG(z), Log(BarnesG(z))), Equal(LogGamma(z), Log(Gamma(z))))", r"\log G(z) = \log\!\left(G(z)\right) \;\mathbin{\operatorname{and}}\; \log \Gamma(z) = \log\!\left(\Gamma(z)\right)"),
    ("DirichletL(s, chi)", r"L(s, \chi)"),
    ("DirichletLambda(s, chi)", r"\Lambda(s, \chi)"),
    ("Implies(GeneralizedRiemannHypothesis, RiemannHypothesis)", r"\operatorname{GRH} \;\implies\; \operatorname{RH}"),
    ("Set(ModularJ(tau), ModularLambda(tau), JacobiTheta(n, z, tau))", r"\left\{j(\tau), \lambda(\tau), \theta_{n}\!\left(z, \tau\right)\right\}"),
    ("Set(WeierstrassP(z, tau), WeierstrassSigma(z, tau), WeierstrassZeta(z, tau))", r"\left\{\wp(z, \tau), \sigma(z, \tau), \zeta(z, \tau)\right\}"),
    ("Mul(ChebyshevT(n, x), ChebyshevU(n, x))", r"T_{n}\!\left(x\right) U_{n}\!\left(x\right)"),
    ("Add(FresnelC(z), FresnelS(z))", r"C(z) + S(z)"),
    ("Div(EisensteinE(Mul(2, n), tau), EisensteinG(Mul(2, n), tau))", r"\frac{E_{2 n}\!\left(\tau\right)}{G_{2 n}\!\left(\tau\right)}"),
    ("Equal(Div(IncompleteBeta(z, a, b), IncompleteBetaRegularized(z, a, b)), BetaFunction(a, b))", r"\frac{\mathrm{B}_{z}\!\left(a, b\right)}{I_{z}\!\left(a, b\right)} = \mathrm{B}(a, b)"),
    ("Set(PolyLog(s, z), HurwitzZeta(s, z), LerchPhi(z, s, a))", r"\left\{\operatorname{Li}_{s}\!\left(z\right), \zeta(s, z), \Phi(z, s, a)\right\}"),
    ("Equal(PartitionsP(n), Mul(Div(1, n), Sum(Mul(DivisorSigma(1, Sub(n, k)), PartitionsP(k)), For(k, 0, Sub(n, 1)))))", r"p(n) = \frac{1}{n} \sum_{k=0}^{n - 1} \sigma_{1}\!\left(n - k\right) p(k)"),
    ("MultiZetaValue(a, b, c)", r"\zeta(a, b, c)"),
    ("RiemannXi(s)", r"\xi(s)"),
    ("Mul(LiouvilleLambda(n), EulerPhi(n), MoebiusMu(n))", r"\lambda(n) \varphi(n) \mu(n)"),
    ("BetaFunction(a, b)", r"\mathrm{B}(a, b)"),
    ("PrimePi(x)", r"\pi(x)"),
    ("Equal(Min(a, b), Neg(Max(Neg(a), Neg(b))))", r"\min(a, b) = -\max\!\left(-a, -b\right)"),
    ("Equal(Arg(z), Div(Pi, 2))", r"\arg(z) = \frac{\pi}{2}"),
    ("NotEqual(Csgn(z), Sign(z))", r"\operatorname{csgn}(z) \ne \operatorname{sgn}(z)"),
    ("Add(Factorial(0), Factorial(1), Div(1, Factorial(-3)), Factorial(Div(1, 2)), Factorial(Factorial(n)), DoubleFactorial(n))", r"0! + 1! + \frac{1}{\left(-3\right)!} + \left(\frac{1}{2}\right)! + \left(n!\right)! + n!!"),
    ("List(Binomial(x, n), RisingFactorial(x, n), FallingFactorial(x, n), StirlingCycle(x, n), StirlingS1(x, n), StirlingS2(x, n))", r"\left[{x \choose n}, \left(x\right)_{n}, \left(x\right)^{\underline{n}}, \left[{x \atop n}\right], s\!\left(x, n\right), \left\{{x \atop n}\right\}\right]"),
    ("Add(BellNumber(5), BernoulliB(5), EulerE(5), Fibonacci(5), HarmonicNumber(5), Prime(5), RiemannZetaZero(5))", r"\operatorname{B}_{5} + B_{5} + E_{5} + F_{5} + H_{5} + p_{5} + \rho_{5}"),
    ("List(LegendreSymbol(p, q), JacobiSymbol(p, q), KroneckerSymbol(p, q))", r"\left[\left(\frac{p}{q}\right), \left(\frac{p}{q}\right), \left(\frac{p}{q}\right)\right]"),
    ("Add(ExpIntegralEi(x), ExpIntegralE(n, x), SinIntegral(x), SinhIntegral(x), CosIntegral(x), CoshIntegral(x), LogIntegral(x))", r"\operatorname{Ei}(x) + E_{n}\!\left(x\right) + \operatorname{Si}(x) + \operatorname{Shi}(x) + \operatorname{Ci}(x) + \operatorname{Chi}(x) + \operatorname{li}(x)"),
    ("Mul(BesselJ(nu, z), BesselI(nu, z), BesselY(nu, z), BesselK(nu, z))", r"J_{\nu}\!\left(z\right) I_{\nu}\!\left(z\right) Y_{\nu}\!\left(z\right) K_{\nu}\!\left(z\right)"),
    ("Equal(AiryAi(AiryAiZero(n)), AiryBi(AiryBiZero(n)), 0)", r"\operatorname{Ai}\!\left(a_{n}\right) = \operatorname{Bi}\!\left(b_{n}\right) = 0"),
    ("Equal(BesselJ(nu, BesselJZero(nu, n)), BesselY(nu, BesselYZero(nu, n)), 0)", r"J_{\nu}\!\left(j_{\nu, n}\right) = Y_{\nu}\!\left(y_{\nu, n}\right) = 0"),
    ("Equal(RiemannZeta(s), Mul(Mul(Mul(Mul(2, Pow(Mul(2, Pi), Sub(s, 1))), Sin(Div(Mul(Pi, s), 2))), Gamma(Sub(1, s))), RiemannZeta(Sub(1, s))))", r"\zeta(s) = 2 {\left(2 \pi\right)}^{s - 1} \sin\!\left(\frac{\pi s}{2}\right) \Gamma\!\left(1 - s\right) \zeta\!\left(1 - s\right)"),
    ("Pow(Div(Pow(DedekindEta(Mul(2, tau)), 2), Mul(DedekindEta(tau), DedekindEta(Mul(4, tau)))), 24)", r"{\left(\frac{\eta^{2}\!\left(2 \tau\right)}{\eta(\tau) \eta\!\left(4 \tau\right)}\right)}^{24}"),
    ("Mul(Mul(Erf(z), Erfc(z)), Erfi(z))", r"\operatorname{erf}(z) \operatorname{erfc}(z) \operatorname{erfi}(z)"),
    ("Mul(EllipticK(m), EllipticE(m), EllipticPi(n, m))", r"K(m) E(m) \Pi(n, m)"),
    ("Mul(IncompleteEllipticE(z, m), IncompleteEllipticF(z, m), IncompleteEllipticPi(n, z, m))", r"E(z, m) F(z, m) \Pi(n, z, m)"),
    ("Add(CarlsonRF(x, y, z), CarlsonRG(x, y, z), CarlsonRJ(x, y, z, w), CarlsonRD(x, y, z), CarlsonRC(x, y))", r"R_F(x, y, z) + R_G(x, y, z) + R_J(x, y, z, w) + R_D(x, y, z) + R_C(x, y)"),
    ("Mul(Hypergeometric0F1(b, z), Hypergeometric0F1Regularized(b, z))", r"\,{}_0F_1(b, z) \,{}_0{\textbf F}_1(b, z)"),
    ("Mul(Hypergeometric1F1(a, b, z), Hypergeometric1F1Regularized(a, b, z))", r"\,{}_1F_1(a, b, z) \,{}_1{\textbf F}_1(a, b, z)"),
    ("Hypergeometric2F0(a, b, z)", r"\,{}_2F_0(a, b, z)"),
    ("Mul(HypergeometricU(a, b, z), HypergeometricUStar(a, b, z))", r"U(a, b, z) U^{*}(a, b, z)"),
    ("Mul(Hypergeometric2F1(a, b, c, z), Hypergeometric2F1Regularized(a, b, c, z))", r"\,{}_2F_1(a, b, c, z) \,{}_2{\textbf F}_1(a, b, c, z)"),
    ("Mul(Hypergeometric1F2(a, b, c, z), Hypergeometric1F2Regularized(a, b, c, z))", r"\,{}_1F_2(a, b, c, z) \,{}_1{\textbf F}_2(a, b, c, z)"),
    ("Mul(Hypergeometric2F2(a, b, c, d, z), Hypergeometric2F2Regularized(a, b, c, d, z))", r"\,{}_2F_2(a, b, c, d, z) \,{}_2{\textbf F}_2(a, b, c, d, z)"),
    ("Mul(Hypergeometric3F2(a, b, c, d, e, z), Hypergeometric3F2Regularized(a, b, c, d, e, z))", r"\,{}_3F_2(a, b, c, d, e, z) \,{}_3{\textbf F}_2(a, b, c, d, e, z)"),
    ("(Hypergeometric2F1Regularized(Div(-1,4),Div(1,4),1/2, (x-1)/2)**2)", r"{\left(\,{}_2{\textbf F}_1\!\left(-\frac{1}{4}, \frac{1}{4}, \frac{1}{2}, \frac{x - 1}{2}\right)\right)}^{2}"),
    ("Add(ZeroMatrix(2), IdentityMatrix(2), HilbertMatrix(2))", r"0_{2} + I_{2} + H_{2}"),
    ("Set(SpecialLinearGroup(n, ZZ), GeneralLinearGroup(n, ZZ))", r"\left\{\operatorname{SL}_{n}\!\left(\mathbb{Z}\right), \operatorname{GL}_{n}\!\left(\mathbb{Z}\right)\right\}"),
    ("Equal(One(QQ), 1)", r"1_{\mathbb{Q}} = 1"),
    ("Equal(Zero(QQ), 0)", r"0_{\mathbb{Q}} = 0"),
    ("List(Polynomials(QQ, x), Polynomials(QQ, x, y), Polynomials(QQ, Tuple()), Polynomials(QQ, Tuple(x)), Polynomials(QQ, Tuple(x, y)))", r"\left[\mathbb{Q}[x], \mathbb{Q}[x, y], \mathbb{Q}[], \mathbb{Q}[x], \mathbb{Q}[x, y]\right]"),
    ("List(Polynomials(QQ, x), PolynomialFractions(QQ, x), FormalPowerSeries(QQ, x), FormalLaurentSeries(QQ, x), FormalPuiseuxSeries(QQ, x))", r"\left[\mathbb{Q}[x], \mathbb{Q}(x), \mathbb{Q}[[x]], \mathbb{Q}(\!(x)\!), \mathbb{Q}\!\left\langle\!\left\langle x \right\rangle\!\right\rangle\right]"),
]



def test_latex(fexpr):
    print("test latex!")
    fexpr.inject(vars=True)
    for formula, expected in latex_test_cases:
        expr = eval(formula)
        latex = expr.latex()
        if latex != expected:
            raise AssertionError("%s:  got '%s', expected '%s'" % (formula, latex, expected))

def latex_report(fexpr):
    fexpr.inject(vars=True)
    formulas = [eval(formula) for formula, expected in latex_test_cases]
    from os.path import expanduser
    from time import clock
    fp = open(expanduser("~/Desktop/latex_report.html"), "w")
    fp.write("""
<!DOCTYPE html>
<html>
<head>
<title>fexpr to LaTeX test sheet</title>
<meta http-equiv="Content-Type" content="text/html;charset=utf-8" >
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
tt { padding: 0.1em; background-color: #f8f8f8; border:1px solid #eee; }
table { border-collapse:collapse; margin: 1em; }
table, th, td { border: 1px solid #aaa; }
th, td { padding:0.1em 0.3em 0.1em 0.3em; }
table { width: 95%; }
.katex { font-size: 1.1em !important; } 
.katex-display { margin:0.1em; padding:0.1em; }
</style>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.css" integrity="sha384-AfEj0r4/OFrOo5t7NnNe46zW/tFgW6x/bCJG8FqQCEo3+Aro6EYUG4+cU+KJWu/X" crossorigin="anonymous">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.js" integrity="sha384-g7c+Jr9ZivxKLnZTDUhnkOnsh30B4H0rpLUpJ4jAIKs4fnJI+sEnkvrMWph2EDg4" crossorigin="anonymous"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/contrib/auto-render.min.js" integrity="sha384-mll67QQFJfxn0IYznZYonOWZ644AWYC+Pt2cHqMaRhXVrursRwvLnLaebdGIlYNa" crossorigin="anonymous"
    onload="renderMathInElement(document.body);"></script>
<script>
  document.addEventListener("DOMContentLoaded", function() {
      renderMathInElement(document.body, {
          delimiters: [
            {left: "$$", right: "$$", display: true},
            {left: "$", right: "$", display: false}
          ]
      });
  });
</script>
</head>
<body>
""")
    output = [formula.latex() for formula in formulas]
    one_big = fexpr("BigLatex")(*formulas)
    t1 = clock()
    one_big_latex = one_big.latex()
    t2 = clock()
    fp.write("""<h1>fexpr to LaTeX test sheet</h1>""")
    fp.write("""<p>Converted %i formulas (%i leaves, %i bytes) to LaTeX in %f seconds.</p>""" % (len(formulas), one_big.num_leaves(), one_big.size_bytes(), (t2-t1)))
    fp.write("""<table>""")
    fp.write("""<tr><th>fexpr</th> <th>Generated LaTeX</th> <th>KaTeX display</th>""")
    for formula, latex in zip(formulas, output):
        fp.write("""<tr>""")
        fp.write("""<td><tt>%s</tt></td>""" % formula)
        fp.write("""<td><tt>%s</tt></td>""" % latex)
        fp.write("""<td>$$%s$$</td>""" % latex)
        fp.write("""</tr>""")
    fp.write("""</table>""")

    fp.write("""<br/><p>Untested builtins:</p> <p><tt>""")
    s = str(one_big)
    for c in '-+()_,"':
        s = s.replace(c, " ")
    used = set(s.split())
    builtins = [name.strip("_") for name in fexpr.builtins()]
    unused = [name for name in builtins if name not in used]
    for name in unused:
        fp.write(name)
        fp.write(" ")
    fp.write("""</tt></p>""")
    fp.write("""</body></html>""")
    fp.close()
