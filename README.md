# arbcmath
Extra C99 double complex transcendental functions with Arb.

Author: Fredrik Johansson (<fredrik.johansson@gmail.com>)

This code is public domain.

## Description

C99 provides support for complex numbers via the standard library header `complex.h`, but only includes a handful of transcendental functions.

This is a simple wrapper of Arb (https://github.com/fredrik-johansson/arb/), exposing many more useful complex transcendental functions (Riemann zeta, polylogarithm, Bessel, incomplete gamma, hypergeometric, Jacobi theta, etc.) in double precision.

All function arguments as well as the return value have the C99 `double complex` type. The output, if finite, is guaranteed to have a relative error that is at most a small multiple of `2^-53` (todo: specify an explicit multiple here), unless the function value is so small that this value underflows to zero. An infinite/NaN return value indicates either that the evaluation has failed to convergence (possibly due to trying to evaluate the function at a point where it is undefined).

Currently, Arb 2.8.0 (or the git master of Arb) is required.

## Provided functions

    ac_exp(z)               Exponential function
    ac_expm1(z)             Accurate exp(z)-1
    ac_log(z)               Natural logarithm
    ac_log1p(z)             Accurate log(1+z)
    ac_sqrt(z)              Square root
    ac_rsqrt(z)             Reciprocal square root
    ac_cbrt(z)              Cube root
    ac_pow(a,b)             Power a^b
  
    ac_sin(z)               Trigonometric functions
    ac_cos(z)
    ac_tan(z)
    ac_cot(z)
    ac_sinpi(z)             Trigonometric functions, argument multiplied by pi
    ac_cospi(z)
    ac_tanpi(z)
    ac_cotpi(z)
    ac_asin(z)              Inverse trigonometric functions
    ac_acos(z)
    ac_atan(z)
    ac_sinh(z)              Hyperbolic functions
    ac_cosh(z)
    ac_tanh(z)
    ac_coth(z)
    ac_asinh(z)             Inverse hyperbolic functions
    ac_acosh(z)
    ac_atanh(z)
    
    ac_gamma(z)             Gamma function
    ac_rgamma(z)            Reciprocal gamma function
    ac_lgamma(z)            Logarithmic gamma function
    ac_digamma(z)           Digamma function
    ac_zeta(s)              Riemann zeta function
    ac_zeta2(s,a)           Hurwitz zeta function
    ac_polygamma(s,z)       Polygamma function
    ac_polylog(s,z)         Polylogarithm
    ac_barnesg(s)           Barnes G-function
    ac_lbarnesg(s)          Logarithmic Barnes G-function
    
    ac_erf(s)               Error function
    ac_erfc(s)              Complementary error function
    ac_erfi(s)              Imaginary error function
    ac_gammaup(s,z)         Upper incomplete gamma function
    ac_expint(s,z)          Generalized exponential integral E
    ac_ei(z)                Exponential integral Ei
    ac_si(z)                Sine integral
    ac_ci(z)                Cosine integral
    ac_shi(z)               Hyperbolic sine integral
    ac_chi(z)               Hyperbolic cosine integral
    ac_li(z)                Logarithmic integral
    ac_lioffset(z)          Offset logarithmic integral
    
    ac_besselj(v,z)         Bessel function J
    ac_bessely(v,z)         Bessel function Y
    ac_besseli(v,z)         Bessel function I
    ac_besselk(v,z)         Bessel function K
    ac_ai(z)                Airy function Ai
    ac_aiprime(z)           Airy function derivative Ai'
    ac_bi(z)                Airy function Bi
    ac_biprime(z)           Airy function derivative Bi'

    ac_hyperu(a,b,z)        Confluent hypergeometric function U
    ac_hyp0f1(a,z)          Confluent hypergeometric function 0F1
    ac_hyp0f1r(a,z)         Regularized confluent hypergeometric function 0F1
    ac_hyp1f1(a,b,z)        Confluent hypergeometric function 1F1
    ac_hyp1f1r(a,b,z)       Regularized confluent hypergeometric function 1F1
    ac_hyp2f1(a,b,c,z)      Hypergeometric function 2F1
    ac_hyp2f1r(a,b,c,z)     Regularized hypergeometric function 2F1
  
    ac_chebyt(n,z)          Chebyshev polynomial/function T
    ac_chebyu(n,z)          Chebyshev polynomial/function U
    ac_jacobip(n,a,b,z)     Jacobi polynomial/function P
    ac_gegenbauerc(n,m,z)   Gegenbauer polynomial/function C
    ac_laguerrel(n,m,z)     Laguerre polynomial/function L
    ac_hermiteh(n,z)        Hermite polynomial/function H
    ac_legenp(n,m,z)        Associated Legendre polynomial/function P
    ac_legenpv(n,m,z)       Associated Legendre polynomial/function P (alt. branch)
    ac_legenq(n,m,z)        Associated Legendre polynomial/function Q
    ac_legenqv(n,m,z)       Associated Legendre polynomial/function Q (alt. branch)
  
    ac_modeta(tau)          Dedekind eta function
    ac_modj(tau)            Modular j-invariant
    ac_modlambda(tau)       Modular lambda function
    ac_moddelta(tau)        Modular delta function
    ac_agm1(z)              Arithmetic-geometric mean of 1 and z
    ac_ellipk(z)            Complete elliptic integral K
    ac_ellipe(z)            Complete elliptic integral E
    ac_ellipp(z,tau)        Weierstrass elliptic function P
    ac_theta1(z,tau)        Jacobi theta function theta1
    ac_theta2(z,tau)        Jacobi theta function theta2
    ac_theta3(z,tau)        Jacobi theta function theta3
    ac_theta4(z,tau)        Jacobi theta function theta4

## Development ideas

Make it a library (don't define everything as inline functions in a header file).

Write some test code.

Write a similar module for real-valued functions.

Write a similar module for long double, GCC quadruple precision, etc.

Allow the user to set different output tolerances.
For example, one could guarantee correct rounding (requires knowing
the exact input-output pairs of the function), or guarantee that
both the real and imaginary parts are accurate separately
(requires knowing where the function is exactly real/imaginary).
The user might also want to lower the precision (say, guaranteeing
only a relative error of 1e-8) as a tradeoff for speed.

Allow more fine-grained control over what to do when convergence fails,
when the final conversion overflows/underflows the exponent range of a
double, etc.

Wrap the Arb functions that compute several functions or function derivatives
simultaneously, e.g. sin+cos, Bessel J+Y, Jacobi theta 1+2+3+4.

Wrap Arb functions that take integer parameters as inputs (e.g. nth root).
