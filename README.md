# Prime Spigot Algorithm
A novel BBP-like spigot algorithm for approximating the nth prime number using the Riemann zeta function's non-trivial zeros and the prime counting function π(x). This implementation inverts an explicit formula for π(x) with oscillatory terms from zeta zeros, achieving reasonable accuracy for large n.

## Overview
This project implements a mutual bootstrap algorithm to approximate the nth prime by inverting the prime counting function π(x). It leverages the logarithmic integral li(x), a correction term, and an oscillatory sum over the first 20 known non-trivial zeros of the Riemann zeta function. The approach is inspired by the Bailey-Borwein-Plouffe (BBP) formula, adapting its series-based methodology to prime number approximation.

### Key Features
- **Novel Approach**: Adapts BBP-like techniques to approximate the nth prime via π(x) inversion.
- **Zeta Zeros**: Incorporates the first 20 non-trivial Riemann zeta zeros for oscillatory terms.
- **Bisection Method**: Uses bisection to invert the π(x) approximation efficiently.
- **Error Analysis**: Includes benchmark tests to evaluate accuracy against known primes.

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/prime-spigot.git
   ```
2. Ensure Python 3.6+ is installed.
3. No external dependencies are required beyond the Python standard library (`math`).

## Usage
Run the script to compute approximations for the nth prime and view benchmark results:

```bash
python bootstrap_primes.py
```

The script outputs a table comparing estimated primes against actual values for specific indices, including absolute and relative errors.

### Example Output
Below is the benchmark output from the script, demonstrating the algorithm's accuracy for various n:

```
Benchmark Results
------------------------------------------------------------
n          Actual       Estimate     Abs Error    Rel Error (%) 
------------------------------------------------------------
15000      163841       164360       519          0.317         
25000      287117       287812       695          0.242         
50000      611953       611517       436          0.071         
75000      951161       952259       1098         0.115         
100000     1299709      1299233      476          0.037         
200000     2750159      2748209      1950         0.071         
------------------------------------------------------------
Average relative error: 0.142%
```

The average relative error is approximately 0.142%, indicating high accuracy for the tested range.

## Methodology
The algorithm approximates the nth prime by:
1. **Computing π(x)**: Uses the logarithmic integral li(x), a correction term (sqrt(x)/log(x)), and an oscillatory sum over zeta zeros with damping for stability.
2. **Inverting π(x)**: Employs bisection to find x such that π(x) ≈ n.
3. **Initial Guess**: Starts with an asymptotic approximation n * (log(n) + log(log(n)) - 1 + ...).
4. **Zeta Zeros**: Incorporates the first 20 non-trivial zeros to capture oscillatory behavior.

### Theoretical Basis
The prime counting function π(x) is approximated using the explicit formula:

π(x) ≈ li(x) - sqrt(x)/log(x) - 2 * sqrt(x) * Σ (cos(γ * log(x) - π/4) * e^(-γ/(log(x) + 20)) / γ)

where γ are the imaginary parts of the zeta zeros, and li(x) is the logarithmic integral.

## Notes
- The algorithm corrects a previous error in the 75000th prime (951161, not 909091).
- Accuracy improves with more zeta zeros, but computational cost increases.
- The damping factor in the oscillatory sum stabilizes the approximation for large x.

## Author
lostdemeter, July 26, 2025
