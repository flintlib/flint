import flint

RR = flint.RR

__exclude = set(["prec"])

__globals = globals()
for __obj in set(dir(RR)) - __exclude:
    if __obj not in __exclude:
        __globals[__obj] = getattr(RR, __obj)

