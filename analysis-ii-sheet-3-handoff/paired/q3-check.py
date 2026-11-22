from fractions import Fraction
# In an incomplete interval distance to its missing endpoint gives the explicit unbounded map.
for n in (2, 5, 16, 125):
    x = Fraction(1,n)
    assert 1/x == n
# Evaluate separated bump functions on an unbounded discrete metric space.
for centre in range(1,35):
    bump_sum = sum(n*max(0,1-4*int(centre!=n)) for n in range(1,35))
    assert bump_sum == centre
print("The missing-endpoint and separated-sequence witnesses attain their prescribed unbounded values.")
