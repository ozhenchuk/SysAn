from functools import lru_cache

@lru_cache(maxsize=256)
def Chebyshev(p, x):
    if p == 0:
        return 1.
    elif p == 1:
        return 2*x - 1
    else: 
        return 2*(2*x - 1)*Chebyshev(p - 1, x) - Chebyshev(p - 2, x)

@lru_cache(maxsize=256)
def Legender(p, x):
    if p == 0:
        return 1.
    elif p == 1:
        return x
    else: 
        return ((2*p - 1)*x*Legender(p - 1, x) - (p - 1)*Legender(p - 2, x))/p

@lru_cache(maxsize=256)
def Lagger(p, x):
    if p == 0:
        return 1.
    elif p == 1:
        return 1 - x
    else: 
        return (2*p - 1 - x)*Lagger(p - 1, x) - (p - 1)*(p - 1)*Lagger(p - 2, x)

@lru_cache(maxsize=256)
def Hermit(p, x):
    if p == 0:
        return 1.
    elif p == 1:
        return x
    else: 
        return (2*x)*Hermit(p - 1, x) - 2*(p - 1)*Hermit(p - 2, x)

