from itertools import count, chain
def isprime(n):
    for p in primes():
        if n % p == 0:
            return False
        if p * p > n:
            return True


def primes(_cache=[2, 3]):
    yield from _cache
    for n in count(_cache[-1] + 2, 2):
        if isprime(n):
            _cache.append(n)
            yield n

for k in primes():
    print(k)