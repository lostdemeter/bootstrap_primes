"""
Mutual Bootstrap Algorithm for Approximating Primes Using Riemann Zeta Zeros

This module implements an approximation to the nth prime by inverting an
explicit formula for the prime counting function π(x), incorporating
oscillatory terms from Riemann zeta zeros. It includes functions for the
logarithmic integral, π(x) approximation, and bisection-based inversion.

Note: Previous versions had a typo in the 75000th prime (909091 instead of 951161),
which inflated the reported error. This is corrected here.

Theoretical basis: The explicit formula relates π(x) to the logarithmic
integral li(x), a correction term, and a sum over zeta zeros.

Author: lostdemeter
Date: July 26, 2025
"""

import math

ZETA_ZEROS = [  # First 20 known zeros for improved accuracy
    14.1347251417, 21.0220396388, 25.0108575801, 30.4248761259, 32.9350615877,
    37.5861781588, 40.9187190121, 43.3270732809, 48.0051508812, 49.7738324777,
    52.9703214777, 56.4462476971, 59.3470440026, 60.8317785246, 65.1125440481,
    67.0798105295, 69.5464017112, 72.0671576745, 75.7046906991, 77.1448400689
]

def approx_li(x, terms=100):
    """
    Approximate the logarithmic integral li(x) using a series expansion for Ei(log x).

    Args:
        x (float): The input value (> 2).
        terms (int): Number of series terms (default: 100).

    Returns:
        float: Approximation to li(x).
    """
    if x <= 2:
        return 0
    z = math.log(x)
    gamma = 0.5772156649015329  # Euler-Mascheroni constant
    res = gamma + math.log(z)
    zk = z
    fact = 1
    for k in range(1, terms):
        fact *= k
        term = zk / (k * fact)
        if term < 1e-12:
            break
        res += term
        zk *= z
    return res

def compute_pi_improved(x, zeros):
    """
    Approximate the prime counting function π(x) using li(x), sqrt(x)/log(x),
    and an oscillatory sum over zeta zeros with damping for stability.

    Args:
        x (float): The input value (> 2).
        zeros (list[float]): Imaginary parts of zeta zeros.

    Returns:
        float: Approximation to π(x).
    """
    if x <= 2:
        return 0
    
    li = approx_li(x)
    logx = math.log(x)
    sqrtx = math.sqrt(x)
    
    osc = 0
    for gamma in zeros:
        if gamma > 0:
            phase = gamma * logx - math.pi / 4
            damping = math.exp(-gamma / (logx + 20))  # Reduced damping
            osc += damping * math.cos(phase) / gamma
    
    osc *= 2 * sqrtx
    return li - sqrtx / logx - osc

def prime_spigot_improved(n, zeros, bisection_iters=40):
    """
    Approximate the nth prime by inverting the π(x) approximation using bisection.

    Args:
        n (int): The index of the prime (> 0).
        zeros (list[float]): Imaginary parts of zeta zeros.
        bisection_iters (int): Number of bisection iterations (default: 40).

    Returns:
        float: Approximation to the nth prime.
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    if n == 1:
        return 2.0
    if n == 2:
        return 3.0
    
    logn = math.log(n)
    loglogn = math.log(logn) if logn > 1 else 0
    
    if n >= 6:
        guess = n * (logn + loglogn - 1 + (loglogn - 2)/logn)
    else:
        guess = n * logn
    
    low = max(2.0, guess * 0.8)
    high = guess * 1.5
    
    while compute_pi_improved(high, zeros) < n:
        high *= 1.5
        if high > 1e18:
            break
    
    for _ in range(bisection_iters):
        mid = (low + high) / 2
        pi_mid = compute_pi_improved(mid, zeros)
        if pi_mid < n:
            low = mid
        else:
            high = mid
    
    return (low + high) / 2

if __name__ == "__main__":
    zeros = ZETA_ZEROS

    # Corrected benchmark test cases (n, actual nth prime)
    test_cases = [
        (15000, 163841),
        (25000, 287117),
        (50000, 611953),
        (75000, 951161),  # Corrected from 909091
        (100000, 1299709),
        (200000, 2750159)
    ]

    print("Benchmark Results")
    print("-" * 60)
    print(f"{'n':<10} {'Actual':<12} {'Estimate':<12} {'Abs Error':<12} {'Rel Error (%)':<14}")
    print("-" * 60)

    errors = []
    for n, actual in test_cases:
        estimate = prime_spigot_improved(n, zeros)
        abs_error = abs(estimate - actual)
        rel_error = (abs_error / actual) * 100
        errors.append(rel_error)
        
        print(f"{n:<10} {actual:<12} {estimate:<12.0f} {abs_error:<12.0f} {rel_error:<14.3f}")

    avg_error = sum(errors) / len(errors)
    print("-" * 60)
    print(f"Average relative error: {avg_error:.3f}%")
