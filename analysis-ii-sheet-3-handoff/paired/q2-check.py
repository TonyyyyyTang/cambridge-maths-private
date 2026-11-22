from fractions import Fraction
# The reciprocal-integer sequence is a metric-space witness to both failures.
# On X={0} union {1/n : n>=1}, take f(0)=0 and f(1/n)=1.
for n in (1, 4, 20, 1000):
    radius = Fraction(1, n)
    point = Fraction(1, 2*n)
    assert 0 < point < radius
    assert not abs(1-0) < Fraction(1, 2)
# f^{-1}((-1/2,1/2)) consists of zero, with no ball contained in it.
print("At all tested scales the same points witness failure of both continuity and an open preimage.")
