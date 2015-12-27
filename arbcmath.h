#include <complex.h>
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_modular.h"
#include "double_extras.h"

#define MAKE_COMPLEX(re,im) ((re) + (im)*I)
#define COMPLEX_NAN MAKE_COMPLEX(D_NAN, D_NAN)

#define CHECK_FINITE(x)                                             \
    if (!isfinite(creal(x)) || !isfinite(cimag(x)))                 \
        return COMPLEX_NAN;                                         \

static __inline__ double complex _acb_init_set_c(acb_t z, double complex v)
{
    acb_init(z);
    acb_set_d_d(z, creal(v), cimag(v));
}

static __inline__ double complex _acb_get_c(const acb_t z)
{
    double re, im;
    re = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    im = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    return MAKE_COMPLEX(re, im);
}

#define ADAPTIVE_REFINEMENT(ACBCALL) \
    {  slong wp; \
    for (wp = 80; ; wp *= 2) { \
        ACBCALL; \
        if (acb_rel_accuracy_bits(out1) >= 53) { \
            result = _acb_get_c(out1); \
            break; \
        } \
        if (wp > 10000) { \
            result = COMPLEX_NAN; \
            break; \
        } \
    } } \

#define DEF_FUNCTION_1(cfun, ACBCALL) \
static __inline__ double complex ac_ ## cfun \
(double complex z1) \
{ \
    acb_t in1, out1; \
    double complex result; \
    CHECK_FINITE(z1) \
    _acb_init_set_c(in1, z1); \
    acb_init(out1); \
    ADAPTIVE_REFINEMENT(ACBCALL) \
    acb_clear(in1); \
    acb_clear(out1); \
    return result; \
}

#define DEF_FUNCTION_2(cfun, ACBCALL) \
static __inline__ double complex ac_ ## cfun \
(double complex z1, double complex z2) \
{ \
    acb_t in1, in2, out1; \
    double complex result; \
    CHECK_FINITE(z1) \
    CHECK_FINITE(z2) \
    _acb_init_set_c(in1, z1); \
    _acb_init_set_c(in2, z2); \
    acb_init(out1); \
    ADAPTIVE_REFINEMENT(ACBCALL) \
    acb_clear(in1); \
    acb_clear(in2); \
    acb_clear(out1); \
    return result; \
}

#define DEF_FUNCTION_3(cfun, ACBCALL) \
static __inline__ double complex ac_ ## cfun \
(double complex z1, double complex z2, double complex z3) \
{ \
    acb_t in1, in2, in3, out1; \
    double complex result; \
    CHECK_FINITE(z1) \
    CHECK_FINITE(z2) \
    CHECK_FINITE(z3) \
    _acb_init_set_c(in1, z1); \
    _acb_init_set_c(in2, z2); \
    _acb_init_set_c(in3, z3); \
    acb_init(out1); \
    ADAPTIVE_REFINEMENT(ACBCALL) \
    acb_clear(in1); \
    acb_clear(in2); \
    acb_clear(in3); \
    acb_clear(out1); \
    return result; \
}

#define DEF_FUNCTION_4(cfun, ACBCALL) \
static __inline__ double complex ac_ ## cfun \
(double complex z1, double complex z2, double complex z3, double complex z4) \
{ \
    acb_t in1, in2, in3, in4, out1; \
    double complex result; \
    CHECK_FINITE(z1) \
    CHECK_FINITE(z2) \
    CHECK_FINITE(z3) \
    CHECK_FINITE(z4) \
    _acb_init_set_c(in1, z1); \
    _acb_init_set_c(in2, z2); \
    _acb_init_set_c(in3, z3); \
    _acb_init_set_c(in4, z4); \
    acb_init(out1); \
    ADAPTIVE_REFINEMENT(ACBCALL) \
    acb_clear(in1); \
    acb_clear(in2); \
    acb_clear(in3); \
    acb_clear(in4); \
    acb_clear(out1); \
    return result; \
}

/* placeholder implementations */
static void
_acb_expm1(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_real(z))
    {
        arb_expm1(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
    }
    else
    {
        acb_exp(res, z, prec);
        acb_sub_ui(res, res, 1, prec);
    }
}

DEF_FUNCTION_1(exp, acb_exp(out1, in1, wp))
DEF_FUNCTION_1(expm1, _acb_expm1(out1, in1, wp))
DEF_FUNCTION_1(log, acb_log(out1, in1, wp))
DEF_FUNCTION_1(log1p, acb_log1p(out1, in1, wp))
DEF_FUNCTION_2(pow, acb_pow(out1, in1, in2, wp))

DEF_FUNCTION_1(sqrt, acb_sqrt(out1, in1, wp))
DEF_FUNCTION_1(rsqrt, acb_rsqrt(out1, in1, wp))
DEF_FUNCTION_1(cbrt, acb_root_ui(out1, in1, 3, wp))

DEF_FUNCTION_1(sin, acb_sin(out1, in1, wp))
DEF_FUNCTION_1(cos, acb_cos(out1, in1, wp))
DEF_FUNCTION_1(tan, acb_tan(out1, in1, wp))
DEF_FUNCTION_1(cot, acb_cot(out1, in1, wp))

DEF_FUNCTION_1(sinpi, acb_sin_pi(out1, in1, wp))
DEF_FUNCTION_1(cospi, acb_cos_pi(out1, in1, wp))
DEF_FUNCTION_1(tanpi, acb_tan_pi(out1, in1, wp))
DEF_FUNCTION_1(cotpi, acb_cot_pi(out1, in1, wp))

DEF_FUNCTION_1(asin, acb_asin(out1, in1, wp))
DEF_FUNCTION_1(acos, acb_acos(out1, in1, wp))
DEF_FUNCTION_1(atan, acb_atan(out1, in1, wp))

DEF_FUNCTION_1(sinh, acb_sinh(out1, in1, wp))
DEF_FUNCTION_1(cosh, acb_cosh(out1, in1, wp))
DEF_FUNCTION_1(tanh, acb_tanh(out1, in1, wp))
DEF_FUNCTION_1(coth, acb_coth(out1, in1, wp))

DEF_FUNCTION_1(asinh, acb_asinh(out1, in1, wp))
DEF_FUNCTION_1(acosh, acb_acosh(out1, in1, wp))
DEF_FUNCTION_1(atanh, acb_atanh(out1, in1, wp))

DEF_FUNCTION_1(gamma, acb_gamma(out1, in1, wp))
DEF_FUNCTION_1(rgamma, acb_rgamma(out1, in1, wp))
DEF_FUNCTION_1(lgamma, acb_lgamma(out1, in1, wp))
DEF_FUNCTION_1(digamma, acb_digamma(out1, in1, wp))
DEF_FUNCTION_1(zeta, acb_zeta(out1, in1, wp))
DEF_FUNCTION_2(zeta2, acb_hurwitz_zeta(out1, in1, in2, wp))
DEF_FUNCTION_2(polygamma, acb_polygamma(out1, in1, in2, wp))
DEF_FUNCTION_2(polylog, acb_polylog(out1, in1, in2, wp))
DEF_FUNCTION_1(barnesg, acb_barnes_g(out1, in1, wp))
DEF_FUNCTION_1(lbarnesg, acb_log_barnes_g(out1, in1, wp))

DEF_FUNCTION_1(erf, acb_hypgeom_erf(out1, in1, wp))
DEF_FUNCTION_1(erfc, acb_hypgeom_erfc(out1, in1, wp))
DEF_FUNCTION_1(erfi, acb_hypgeom_erfi(out1, in1, wp))
DEF_FUNCTION_2(gammaup, acb_hypgeom_gamma_upper(out1, in1, in2, 0, wp))
DEF_FUNCTION_2(expint, acb_hypgeom_expint(out1, in1, in2, wp))
DEF_FUNCTION_1(ei, acb_hypgeom_ei(out1, in1, wp))
DEF_FUNCTION_1(si, acb_hypgeom_si(out1, in1, wp))
DEF_FUNCTION_1(ci, acb_hypgeom_ci(out1, in1, wp))
DEF_FUNCTION_1(shi, acb_hypgeom_shi(out1, in1, wp))
DEF_FUNCTION_1(chi, acb_hypgeom_chi(out1, in1, wp))
DEF_FUNCTION_1(li, acb_hypgeom_li(out1, in1, 0, wp))
DEF_FUNCTION_1(lioffset, acb_hypgeom_li(out1, in1, 1, wp))

DEF_FUNCTION_2(besselj, acb_hypgeom_bessel_j(out1, in1, in2, wp))
DEF_FUNCTION_2(bessely, acb_hypgeom_bessel_y(out1, in1, in2, wp))
DEF_FUNCTION_2(besseli, acb_hypgeom_bessel_i(out1, in1, in2, wp))
DEF_FUNCTION_2(besselk, acb_hypgeom_bessel_k(out1, in1, in2, wp))

DEF_FUNCTION_1(ai, acb_hypgeom_airy(out1, NULL, NULL, NULL, in1, wp))
DEF_FUNCTION_1(aiprime, acb_hypgeom_airy(NULL, out1, NULL, NULL, in1, wp))
DEF_FUNCTION_1(bi, acb_hypgeom_airy(NULL, NULL, out1, NULL, in1, wp))
DEF_FUNCTION_1(biprime, acb_hypgeom_airy(NULL, NULL, NULL, out1, in1, wp))

DEF_FUNCTION_2(hyp0f1, acb_hypgeom_0f1(out1, in1, in2, 0, wp))
DEF_FUNCTION_2(hyp0f1r, acb_hypgeom_0f1(out1, in1, in2, 1, wp))
DEF_FUNCTION_3(hyp1f1, acb_hypgeom_m(out1, in1, in2, in3, 0, wp))
DEF_FUNCTION_3(hyp1f1r, acb_hypgeom_m(out1, in1, in2, in3, 1, wp))
DEF_FUNCTION_4(hyp2f1, acb_hypgeom_2f1(out1, in1, in2, in3, in4, 0, wp))
DEF_FUNCTION_4(hyp2f1r, acb_hypgeom_2f1(out1, in1, in2, in3, in4, 1, wp))

DEF_FUNCTION_2(chebyt, acb_hypgeom_chebyshev_t(out1, in1, in2, wp))
DEF_FUNCTION_2(chebyu, acb_hypgeom_chebyshev_u(out1, in1, in2, wp))
DEF_FUNCTION_4(jacobip, acb_hypgeom_jacobi_p(out1, in1, in2, in3, in4, wp))
DEF_FUNCTION_3(gegenbauerc, acb_hypgeom_gegenbauer_c(out1, in1, in2, in3, wp))
DEF_FUNCTION_3(laguerrel, acb_hypgeom_laguerre_l(out1, in1, in2, in3, wp))
DEF_FUNCTION_2(hermiteh, acb_hypgeom_hermite_h(out1, in1, in2, wp))
DEF_FUNCTION_3(legenp, acb_hypgeom_legendre_p(out1, in1, in2, in3, 0, wp))
DEF_FUNCTION_3(legenpv, acb_hypgeom_legendre_p(out1, in1, in2, in3, 1, wp))
DEF_FUNCTION_3(legenq, acb_hypgeom_legendre_q(out1, in1, in2, in3, 0, wp))
DEF_FUNCTION_3(legenqv, acb_hypgeom_legendre_q(out1, in1, in2, in3, 1, wp))

DEF_FUNCTION_1(modeta, acb_modular_eta(out1, in1, wp))
DEF_FUNCTION_1(modj, acb_modular_j(out1, in1, wp))
DEF_FUNCTION_1(modlambda, acb_modular_lambda(out1, in1, wp))
DEF_FUNCTION_1(moddelta, acb_modular_delta(out1, in1, wp))

DEF_FUNCTION_1(agm1, acb_agm1(out1, in1, wp))
DEF_FUNCTION_1(ellipk, acb_modular_elliptic_k(out1, in1, wp))
DEF_FUNCTION_1(ellipe, acb_modular_elliptic_e(out1, in1, wp))
DEF_FUNCTION_2(ellipp, acb_modular_elliptic_p(out1, in1, in2, wp))

static void
_acb_theta1(acb_t res, const acb_t z, const acb_t tau, slong prec)
{
    acb_t a, b, c; acb_init(a); acb_init(b); acb_init(c);
    acb_modular_theta(res, a, b, c, z, tau, prec);
    acb_clear(a); acb_clear(b); acb_clear(c);
}

static void
_acb_theta2(acb_t res, const acb_t z, const acb_t tau, slong prec)
{
    acb_t a, b, c; acb_init(a); acb_init(b); acb_init(c);
    acb_modular_theta(a, res, b, c, z, tau, prec);
    acb_clear(a); acb_clear(b); acb_clear(c);
}

static void
_acb_theta3(acb_t res, const acb_t z, const acb_t tau, slong prec)
{
    acb_t a, b, c; acb_init(a); acb_init(b); acb_init(c);
    acb_modular_theta(a, b, res, c, z, tau, prec);
    acb_clear(a); acb_clear(b); acb_clear(c);
}

static void
_acb_theta4(acb_t res, const acb_t z, const acb_t tau, slong prec)
{
    acb_t a, b, c; acb_init(a); acb_init(b); acb_init(c);
    acb_modular_theta(a, b, c, res, z, tau, prec);
    acb_clear(a); acb_clear(b); acb_clear(c);
}

DEF_FUNCTION_2(theta1, _acb_theta1(out1, in1, in2, wp))
DEF_FUNCTION_2(theta2, _acb_theta2(out1, in1, in2, wp))
DEF_FUNCTION_2(theta3, _acb_theta3(out1, in1, in2, wp))
DEF_FUNCTION_2(theta4, _acb_theta4(out1, in1, in2, wp))

