from fractions import Fraction
# Realize the construction on the first coordinates of finitely supported sequences.
for n in range(1, 65):
    coefficient = Fraction(1, n)
    old_norm = abs(coefficient)
    phi = n * coefficient
    new_norm = old_norm + abs(phi)
    assert phi == 1 and new_norm / old_norm == n + 1
print("The finite restrictions exhibit unbounded norm ratios and a sequence separating the two norms.")
