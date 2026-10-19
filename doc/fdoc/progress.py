"""
Progress reporting on standard error.

On a terminal, a stage with a known amount of work is shown as a counter
that is updated in place ("rendering pages  4123/8595  48%"); elsewhere
(make logs, CI) each stage prints one line when it finishes. Set
FDOC_QUIET=1 (or pass --quiet) to silence it.
"""

import os
import sys
import time

QUIET = bool(os.environ.get("FDOC_QUIET"))


class Progress:
    def __init__(self, label, total=None):
        self.label = label
        self.total = total
        self.count = 0
        self.t0 = time.time()
        self.last = 0.0
        self.tty = sys.stderr.isatty()
        if not QUIET and self.tty:
            self._show()

    def step(self, n=1):
        self.count += n
        if QUIET or not self.tty:
            return
        now = time.time()
        if now - self.last >= 0.1 or (self.total and self.count >= self.total):
            self.last = now
            self._show()

    def _show(self):
        if self.total:
            text = "%-30s %*d/%d  %3d%%" % (self.label, len(str(self.total)), self.count, self.total,
                                          100 * self.count // max(self.total, 1))
        else:
            text = "%-30s ..." % self.label
        sys.stderr.write("\r\x1b[K" + text)
        sys.stderr.flush()

    def done(self, note=""):
        if QUIET:
            return
        dt = time.time() - self.t0
        what = ("%d " % self.count if self.total else "")
        line = "%-30s %sdone in %.1f s%s" % (self.label, what, dt, (" (%s)" % note) if note else "")
        if self.tty:
            sys.stderr.write("\r\x1b[K" + line + "\n")
        else:
            sys.stderr.write(line + "\n")
        sys.stderr.flush()


def stage(label):
    """A stage without a counter: Progress(label) with done() at the end."""
    return Progress(label)
