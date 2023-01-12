import flint2

CC = flint2.CC

__exclude = set(["prec"])

__globals = globals()
for __obj in set(dir(CC)) - __exclude:
    if __obj not in __exclude:
        __globals[__obj] = getattr(CC, __obj)

