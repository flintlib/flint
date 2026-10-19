"""
Comparison of the output of an example program with the expected output.

By default the output must match exactly (up to trailing whitespace at the
end). Flags on the ```output fence relax this for outputs that legitimately
differ between platforms or runs (similar to the "# random", "# abs tol"
and "# rel tol" markers of SageMath's doctests):

    ```output radius-tol         radii of printed balls (the number after
                                 "+/-") may differ by a factor up to 10;
    ```output radius-tol=4       ... up to the given factor
    ```output rel-tol=1e-6       every number may differ by this relative
    ```output abs-tol=1e-10      ... or absolute amount (both may be given)
    ```output random             the output is not compared (the program
                                 must still compile and run successfully)

With a tolerance, the text between the numbers must still match exactly,
and so must the number of lines and of numbers. Ball radii are typically
not reproducible bit for bit: they are rigorous bounds, but how tight they
are can depend on the compiler (for example on whether it contracts
floating-point expressions into fused multiply-adds) and on the algorithms
selected for the machine.
"""

import re

FLAGS = ("radius-tol", "rel-tol", "abs-tol", "random")
DEFAULT_RADIUS_FACTOR = 10.0

RE_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")
RE_RADIUS_BEFORE = re.compile(r"\+/-\s*$")


def check_flags(flags):
    """A list of problems with the flags of an output block."""
    problems = []
    for k, v in flags.items():
        if k not in FLAGS:
            problems.append("unknown flag '%s' on an output block (known: %s)" % (k, ", ".join(FLAGS)))
            continue
        if k == "random":
            if v is not True:
                problems.append("'random' takes no value")
            continue
        if v is True:
            if k != "radius-tol":
                problems.append("'%s' needs a value, as in %s=1e-10" % (k, k))
            continue
        try:
            x = float(v)
        except ValueError:
            problems.append("'%s=%s' is not a number" % (k, v))
            continue
        if x < 0 or (k == "radius-tol" and x < 1):
            problems.append("'%s=%s' is out of range" % (k, v))
    return problems


def describe(flags):
    """A short note shown next to the expected output, or ''."""
    if not flags:
        return ""
    if flags.get("random"):
        return "may differ from run to run"
    if "rel-tol" in flags or "abs-tol" in flags:
        return "numbers may differ slightly between platforms"
    if "radius-tol" in flags:
        return "radii may differ slightly between platforms"
    return ""


def _split(line):
    """Split a line into text pieces and numbers: [t0, n1, t1, n2, t2, ...]."""
    parts = []
    pos = 0
    for m in RE_NUMBER.finditer(line):
        parts.append(line[pos:m.start()])
        parts.append(m.group(0))
        pos = m.end()
    parts.append(line[pos:])
    return parts


def compare(want, got, flags=None):
    """Whether `got` matches the expected output `want` under `flags`."""
    flags = flags or {}
    want = want.rstrip()
    got = got.rstrip()
    if flags.get("random"):
        return True
    if want == got:
        return True
    tol = {k: v for k, v in flags.items() if k in ("radius-tol", "rel-tol", "abs-tol")}
    if not tol:
        return False
    factor = None
    if "radius-tol" in tol:
        factor = DEFAULT_RADIUS_FACTOR if tol["radius-tol"] is True else float(tol["radius-tol"])
    rel = float(tol.get("rel-tol", 0))
    ab = float(tol.get("abs-tol", 0))
    wl = want.split("\n")
    gl = got.split("\n")
    if len(wl) != len(gl):
        return False
    for w, g in zip(wl, gl):
        wp = _split(w.rstrip())
        gp = _split(g.rstrip())
        if len(wp) != len(gp):
            return False
        for i, (a, b) in enumerate(zip(wp, gp)):
            if i % 2 == 0:
                if a != b:
                    return False
                continue
            if a == b:
                continue
            x, y = float(a), float(b)
            if factor is not None and RE_RADIUS_BEFORE.search(wp[i - 1]):
                if x == 0 or y == 0:
                    if x != y:
                        return False
                elif not (x / factor <= y <= x * factor):
                    return False
                continue
            if (rel or ab) and abs(x - y) <= ab + rel * abs(x):
                continue
            return False
    return True
