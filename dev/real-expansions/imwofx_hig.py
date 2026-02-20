#!/usr/bin/env python3

"""
Compute the first 11 terms of the sequence:
  a_0 = 1/sqrt(pi)
  a_N = (2N-1)!! / 2^N / sqrt(pi)  for N >= 1

where (2N-1)!! = 1 * 3 * 5 * ... * (2N-1) is the double factorial.

Using arb/flint for arbitrary precision arithmetic.
"""

from flint import arb, ctx

# Set precision for 18 decimal digits
# 18 decimal digits ≈ 60 bits (log2(10^18) ≈ 59.8)
ctx.prec = 80  # Use extra bits for safety

def compute_a_terms(n_terms: int) -> list[float]:
    """
    Compute the first n_terms of a_N.

    Returns:
        List of a_N values as floats with 18 digits precision
    """
    results = []

    # Compute 1/sqrt(pi) once
    sqrt_pi = arb.pi().sqrt()
    inv_sqrt_pi = 1 / sqrt_pi

    # a_0 = 1/sqrt(pi)
    a0 = inv_sqrt_pi
    results.append(float(a0))
    print(f"a_0  = {float(a0):.18e}")

    # For N >= 1: a_N = (2N-1)!! / 2^N / sqrt(pi)
    # We can compute recursively:
    # a_N = a_{N-1} * (2N-1) / 2

    a_prev = a0
    for N in range(1, n_terms):
        # (2N-1)!! / 2^N = (2N-1)!! / 2^{N-1} * (2N-1) / 2
        #                = (2(N-1)-1)!! / 2^{N-1} * (2N-1) / 2
        factor = arb(2 * N - 1) / arb(2)
        a_N = a_prev * factor

        results.append(float(a_N))
        print(f"a_{N:<2} = {float(a_N):.18e}")

        a_prev = a_N

    return results


if __name__ == "__main__":
    print("Computing first 11 terms of a_N:\n")
    terms = compute_a_terms(11)

    print("\nVerification:")
    print(f"Number of terms computed: {len(terms)}")
    print(f"All values positive: {all(x > 0 for x in terms)}")
    print(f"Sequence is decreasing: {all(terms[i] > terms[i+1] for i in range(len(terms)-1))}")
