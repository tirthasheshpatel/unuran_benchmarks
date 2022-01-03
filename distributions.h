#include <unuran.h>

/* For some reason, UNU.RAN `undef`s M_PI. Define it again, if it is `undef`ed. */
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef ERROR_HANDLER
#undef ERROR_HANDLER
#endif

#define ERROR_HANDLER(args...) do { fprintf(stderr, "Error in file '%s' at line %d: ", __FILE__, __LINE__); \
                                    fprintf(stderr, args); exit(EXIT_FAILURE); } while (0)
#define WARNING_HANDLER(args...) do { fprintf(stderr, "Warning: " args); } while (0)


typedef UNUR_PAR *(*UNUR_METHOD_SIGNATURE)(const UNUR_DISTR *);


// Some functions from the CEPHES library.
extern double ndtr(double);                         // CDF of the normal distribution
extern double igam(double, double);                 // Incomplete gamma integral: 0 ... x
extern double igamc(double, double);                // Incomplete gamma integral: x ... inf
extern double btdtr(double, double, double);        // CDF of beta distribution

// Modified Bessel function of the second kind (from AMOS library)
extern void zbesk_(double *Zr, double *Zi, double *FNU, int *KODE, int *N, double *CYr, double *CYi, int *NZ, int *IERR);

// Numerical integration routines from QUADPACK
// dqagse_ => numerical integration with finite bounds
// dqagie_ => numerical integration with infinite bounds
typedef double quadpack_f_t(double *);

extern void dqagse_(quadpack_f_t f, double *a, double *b, double *epsabs, double *epsrel, int *limit, double *result,
                    double *abserr, int *neval, int *ier, double *alist, double *blist, double *rlist, double *elist,
                    int *iord, int *last);

extern void dqagie_(quadpack_f_t f, double *bound, int *inf, double *epsabs, double *epsrel, int *limit,
                    double *result, double *abserr, int *neval, int *ier, double *alist, double *blist,
                    double *rlist, double *elist, int *iord, int *last);


typedef struct {
    double (*pdf)(double, const UNUR_DISTR *);
    double (*cdf)(double, const UNUR_DISTR *);
    double (*dpdf)(double, const UNUR_DISTR *);
    double (*mode)(const UNUR_DISTR *);
    double *support;
    double *params;
    int num_params;
} cont_dist;

typedef struct {
    double (*pmf)(int, const UNUR_DISTR *);
    double (*cdf)(int, const UNUR_DISTR *);
    double (*mode)(const UNUR_DISTR *);
    double *support;
    double *params;
    int num_params;
} discr_dist;


double _quadpack_test(double *x)
{
    return 1/sqrt(2*M_PI) * exp(-0.5 * (*x)*(*x));
}


/**
 * @brief Modified Bessel function of the second kind of real order v
 *
 * @param v Order of Bessel functions
 * @param z Argument at which to evaluate the Bessel functions
 * @return double The result
 */
double kv(double v, double z)
{
    double Zr = z, Zi = 0, FNU = v;
    int KODE = 1, N = 1;
    double CYr, CYi;
    int NZ, IERR;

    /* K_v == K_{-v} even for non-integer v */
    if (v < 0) v = -v;

    zbesk_(&Zr, &Zi, &FNU, &KODE, &N, &CYr, &CYi, &NZ, &IERR);

    if (IERR == 2) { /* overflow */
        return INFINITY;
    }

    else if (IERR) {
        switch (IERR) {
        case 1:
            ERROR_HANDLER("input error: v=%lf, z=%lf\n", v, z);
        case 3:
        case 4:
            WARNING_HANDLER("inaccurate result: v=%lf, z=%lf\n", v, z);
            break;
        case 5:
            ERROR_HANDLER("termination condition not met: v=%lf, z=%lf\n", v, z);
        default:
            ERROR_HANDLER("internal error: this should not happen, please report!\n");
        }
    }

    return CYr;
}

/**
 * @brief Integrate a function.
 *
 * @param func Function with signature `double f(double *)`
 * @param lb Lower bound of the integral
 * @param ub Upper bound of the integral
 * @return double Area under the curve
 */
double quad(quadpack_f_t func, double lb, double ub)
{
    int      limit = 50, infbounds = 0;
    double   bound = 0, epsabs = 1.49e-8, epsrel = 1.49e-8, a = lb, b = ub;

    if (b != INFINITY && a != -INFINITY) {}
    else if (b == INFINITY && a != -INFINITY) {
        infbounds = 1;
        bound = a;
    }
    else if (b == INFINITY && a == -INFINITY) {
        infbounds = 2;
        bound = 0;
    }
    else if (b != INFINITY && a == -INFINITY) {
        infbounds = -1;
        bound = b;
    }
    else {
        ERROR_HANDLER("Infinity comparisons don't work for you.\n");
    }

    int      neval, ier, last, iord[limit];
    double   result, abserr, alist[limit], blist[limit], rlist[limit], elist[limit];

    if (infbounds) {
        dqagie_(func, &bound, &infbounds, &epsabs, &epsrel, &limit,
                &result, &abserr, &neval, &ier, alist, blist, rlist,
                elist, iord, &last);
    }
    else {
        dqagse_(func, &a, &b, &epsabs, &epsrel, &limit, &result,
                &abserr, &neval, &ier, alist, blist, rlist, elist,
                iord, &last);
    }

    if (ier) {
        switch (ier) {
        case 1:
            ERROR_HANDLER("maximum number of subintervals exceeded\n");
        case 2:
            WARNING_HANDLER("round-off errors high\n");
            break;
        case 3:
            ERROR_HANDLER("extremely bad integrad\n");
        case 4:
            ERROR_HANDLER("the algorithm does not converge\n");
        case 5:
            ERROR_HANDLER("the integral is probably divergent\n");
        case 6:
            ERROR_HANDLER("the input is invalid. Check function signature, lb=%lf, ub=%lf\n", a, b);
        default:
            ERROR_HANDLER("internal error: this should not happen, please report!\n");
        }
    }

    return result;
}

double stdnormal_pdf(double x, const UNUR_DISTR *distr)
{
    return 1/sqrt(2*M_PI) * exp(-0.5 *x*x);
}

double stdnormal_dpdf(double x, const UNUR_DISTR *distr)
{
    return -1/sqrt(2*M_PI) * x * exp(-0.5 *x*x);
}

double stdnormal_cdf(double x, const UNUR_DISTR *distr)
{
    return (1 + erf(x/sqrt(2))) / 2;
}

double stdnormal_mode(const UNUR_DISTR *mode)
{
    return 0;
}

double stdnormal_support[] = {-INFINITY, INFINITY};


double gamma_pdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double a = params[0];

    return pow(x, a-1) * exp(-x) / gamma(a);
}

double gamma_dpdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double a = params[0];

    const double first = (a-1) * pow(x, a-2) * exp(-x) / gamma(a),
                 second = -pow(x, a-1) * exp(-x) / gamma(a);
    return first + second;
}

double gamma_cdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double a = params[0];

    return igam(a, x);
}

double gamma_mode(const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double a = params[0];

    if (a < 1) return 0;

    return a-1;
}

double gamma_support[] = {0, INFINITY};


double beta_pdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double a = params[0], b = params[1];

    return gamma(a+b) * pow(x, a-1) * pow(1-x, b-1) / gamma(a) / gamma(b);
}

double beta_dpdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double a = params[0], b = params[1];

    const double z = gamma(a) * gamma(b) / gamma(a+b);

    return ((a-1) * pow(x, a-2) * pow(1-x, b-1) - pow(x, a-1) * (b-1) * pow(1-x, b-2)) / z;
}

double beta_cdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double a = params[0], b = params[1];

    return btdtr(a, b, x);
}

double beta_mode(const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double a = params[0], b = params[1];

    if (a > 1 && b > 1) return (a-1) / (a+b-2);
    if (a == 1 && b == 1) return 0.5;
    if (a <= 1 && b > 1) return 0;
    return 1;
}

double beta_support[] = {0, 1};


double gennorm_pdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double beta = params[0];
    const double z = 2 * gamma(1/beta) / beta;
    return exp(-pow(fabs(x), beta)) / z;
}

double gennorm_dpdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double beta = params[0];
    const double z = 2 * gamma(1/beta) / beta;
    if (x > 0) {
        return ((exp(-pow(x, beta)) / z)
                * (-beta * pow(x, beta-1)));
    }
    return ((exp(-pow(-x, beta)) / z)
            * (beta * pow(-x, beta-1)));
}

double gennorm_cdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double beta = params[0];
    const double c = 0.5 * (x < 0 ? -1: 1);
    return (0.5 + c) - c * igamc(1/beta, pow(fabs(x), beta));
}

double gennorm_mode(const UNUR_DISTR *distr)
{
    return 0;
}

double gennorm_support[] = {-INFINITY, INFINITY};


double nakagami_pdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double nu = params[0];

    const double z = gamma(nu) / 2 / pow(nu, nu);

    return (pow(x, 2*nu-1) * exp(-nu * x*x)) / z;
}

double nakagami_dpdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double nu = params[0];

    const double z = gamma(nu) / 2 / pow(nu, nu);

    const double first = (2*nu-1) * pow(x, 2*nu-2) * exp(-nu * x*x),
                 second = (pow(x, 2*nu-1) * exp(-nu * x*x) * (-nu * 2*x));
    return (first + second) / z;
}

double nakagami_cdf(double x, const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double nu = params[0];

    return igam(nu, nu*x*x);
}

double nakagami_mode(const UNUR_DISTR *distr)
{
    const double *params;
    unur_distr_cont_get_pdfparams(distr, &params);
    const double nu = params[0];
    return sqrt(2) / 2 * sqrt((2*nu - 1)/nu);
}

double nakagami_support[] = {0, INFINITY};


/**
 * @brief Testing the CEPHES, QUADPACK, and AMOS routines
 */
void _test()
{
    printf("Bessel's modified function of order %lf at %lf: %lf\n", 1., 1., kv(1, 1));
    printf("CDF of normal distribution between -1. and 1.: %lf\n", quad(_quadpack_test, -1, 1));
    printf("Area under the normal distribution: %lf\n", quad(_quadpack_test, -INFINITY, INFINITY));
}
