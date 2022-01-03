/*

Author: Tirth Patel

Full UNU.RAN Benchmark Suite for TDR, HINV, PINV, and SROU.

Distributions benchmarked:

- stdnormal     : Standard Normal distribution                (continuous)
- gammma        : Gamma distribution                          (continuous)
- beta          : Beta distribution                           (continuous)
- gennorm       : Generalized normal distribution             (continuous)
- nakagami      : Nakagami distribution                       (continuous)
- hypergeom     : Hypergeometric distribution                 (discrete)         TODO

*/

#define _ISOC99_SOURCE 1
#define _XOPEN_SOURCE
#include <string.h>
#include "FLA_Clock.h"
#include "distributions.h"
#ifdef USE_RNGSTREAMS
    #include <unuran_urng_rngstreams.h>
#else
    #include "mt19937.h"
#endif



/**
 * @brief Kahan-Babushka-Neumaier Sum
 * 
 * Reference: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 * 
 * @param input Input list of doubles
 * @param n Size of the input
 * @return double Sum
 */
double fsum(double *input, const size_t n)
{
    double sum = 0.0, c = 0.0;

    for(size_t i = 0 ; i < n ; ++i) {
        double t = sum + input[i];
        if (fabs(sum) > fabs(input[i])) c += (sum - t) + input[i];
        else                            c += (input[i] - t) + sum;
        sum = t;
    }

    return sum + c;
}


double fmean(double *input, const size_t n)
{
    return fsum(input, n) / n;
}

double fstddev(double *input, const size_t n)
{
    double loc = fmean(input, n);
    double *stddev_input = (double *) malloc(n * sizeof(double));

    if (!stddev_input) {
        fprintf(stderr, "fatal error: can't allocate %lu bytes\n", n * sizeof(double));
        return -1.0;
    }

    for (size_t i = 0 ; i < n ; ++i) stddev_input[i] = (input[i] - loc) * (input[i] - loc);

    double stddev = sqrt(fsum(stddev_input, n) / n);

    free(stddev_input);

    return stddev;
}


void setup_timeit(cont_dist distribution, UNUR_METHOD_SIGNATURE unur_method_new)
{
    UNUR_DISTR *distr;
    UNUR_URNG *urng;

    double start, end, time_taken;

    #ifndef USE_RNGSTREAMS
        mt19937_state state;
        mt19937_seed(&state, 0x0);
    #endif

    distr = unur_distr_cont_new();

    if (!distr) {
        fprintf(stderr, "fatal error: failed to initialize the distribution object\n");
        exit(EXIT_FAILURE);
    }

    unur_distr_cont_set_pdfparams(distr, distribution.params, distribution.num_params);
    unur_distr_cont_set_pdf(distr, distribution.pdf);
    unur_distr_cont_set_dpdf(distr, distribution.dpdf);
    unur_distr_cont_set_cdf(distr, distribution.cdf);
    unur_distr_cont_set_domain(distr, distribution.support[0], distribution.support[1]);
    if (unur_method_new == unur_srou_new) {
        unur_distr_cont_set_pdfarea(distr, 1);
        unur_distr_cont_set_mode(distr, distribution.mode(distr));
    }

    #ifdef USE_RNGSTREAMS
        urng = unur_urng_rngstream_new("urng");
    #else
        // use NumPy's MT19937
        urng = unur_urng_new(mt19937_next_double, &state);
    #endif

    if (!urng) {
        fprintf(stderr, "fatal error: urng object cannot be intialized!\n");
        exit(EXIT_FAILURE);
    }

    size_t number = 0;

    time_taken = 0;

    do {
        start = FLA_Clock();

        UNUR_PAR *par;
        UNUR_GEN *gen;

        par = unur_method_new(distr);

        if (!par) {
            fprintf(stderr, "fatal error: parameter object cannot be intialized!\n");
            exit(EXIT_FAILURE);
        }

        unur_set_urng(par, urng);

        gen = unur_init(par);

        if (!gen) {
            fprintf(stdout, "Failed!\n");
            fflush(stdout);
            unur_distr_free(distr);
            unur_urng_free(urng);
            return;
        }

        unur_free(gen);

        end = FLA_Clock();

        time_taken += end - start;

        ++number;
    } while (time_taken < 0.2);

    number = (size_t)pow(10, ceil(log10((double)number)));

    double time_taken_array[7];

    for (size_t rep = 0 ; rep < 7 ; ++rep) {
        start = FLA_Clock();
        for (int n = 0 ; n < number ; ++n) {
            UNUR_PAR *par;
            UNUR_GEN *gen;

            par = unur_method_new(distr);

            if (!par) {
                fprintf(stderr, "fatal error: parameter object cannot be intialized!\n");
                exit(EXIT_FAILURE);
            }

            unur_set_urng(par, urng);

            gen = unur_init(par);

            if (!gen) {
                fprintf(stdout, "Failed!\n");
                fflush(stdout);
                unur_distr_free(distr);
                unur_urng_free(urng);
                return;
            }

            unur_free(gen);
        }
        end = FLA_Clock();
        time_taken_array[rep] = (end - start) / number;
    }

    double mean = fmean(time_taken_array, 7);
    double stddev = fstddev(time_taken_array, 7);

    const double scaling[] = {1, 1e3, 1e6, 1e9};

    const char *units[] = {"s", "ms", "µs", "ns"};

    // order = min(-int(math.floor(math.log10(timespan)) // 3), 3)
    int mean_order = 0;
    while (mean * scaling[mean_order] < 1) ++mean_order;
    if (mean_order > 3) mean_order = 3;

    fprintf(stdout, "%.3g %s", mean * scaling[mean_order], units[mean_order]);

    fprintf(stdout, " ± ");

    int stddev_order = 0;
    while (stddev * scaling[stddev_order] < 1) ++stddev_order;
    if (stddev_order > 3) stddev_order = 3;

    fprintf(stdout, "%.3g %s", stddev * scaling[stddev_order], units[stddev_order]);

    fprintf(stdout, " per loop (mean ± std. dev. of %d runs, %lu loop%c each)\n",
                    7, number, (number <= 1 ? '\0' : 's'));

    fflush(stdout);

    unur_distr_free(distr);
    unur_urng_free(urng);
}

void sampling_timeit(cont_dist distribution, UNUR_METHOD_SIGNATURE unur_method_new)
{
    size_t num_samples = 1000000;

    UNUR_DISTR *distr;
    UNUR_URNG *urng;
    UNUR_PAR *par;
    UNUR_GEN *gen;

    double start, end, time_taken;

    #ifndef USE_RNGSTREAMS
        mt19937_state state;
        mt19937_seed(&state, 0x0);
    #endif

    distr = unur_distr_cont_new();

    if (!distr) {
        fprintf(stderr, "fatal error: failed to initialize the distribution object\n");
        exit(EXIT_FAILURE);
    }

    unur_distr_cont_set_pdfparams(distr, distribution.params, distribution.num_params);
    unur_distr_cont_set_pdf(distr, distribution.pdf);
    unur_distr_cont_set_dpdf(distr, distribution.dpdf);
    unur_distr_cont_set_cdf(distr, distribution.cdf);
    unur_distr_cont_set_domain(distr, distribution.support[0], distribution.support[1]);
    if (unur_method_new == unur_srou_new) {
        unur_distr_cont_set_pdfarea(distr, 1);
        unur_distr_cont_set_mode(distr, distribution.mode(distr));
    }

    #ifdef USE_RNGSTREAMS
        urng = unur_urng_rngstream_new("urng");
    #else
        // use NumPy's MT19937
        urng = unur_urng_new(mt19937_next_double, &state);
    #endif

    if (!urng) {
        fprintf(stderr, "fatal error: urng object cannot be intialized!\n");
        exit(EXIT_FAILURE);
    }

    par = unur_method_new(distr);

    if (!par) {
        fprintf(stderr, "fatal error: parameter object cannot be intialized!\n");
        exit(EXIT_FAILURE);
    }

    unur_set_urng(par, urng);

    gen = unur_init(par);

    if (!gen) {
        fprintf(stdout, "Failed!\n");
        fflush(stdout);
        unur_distr_free(distr);
        unur_urng_free(urng);
        return;
    }

    unur_distr_free(distr);

    size_t number = 0;

    time_taken = 0;

    double *rvs = (double *) malloc(num_samples * sizeof(double));

    if (!rvs) {
        fprintf(stderr, "fatal error: can't allocate %lu bytes\n", num_samples * sizeof(double));
        exit(EXIT_FAILURE);
    }

    do {
        start = FLA_Clock();

        for (size_t i = 0 ; i < num_samples ; ++i) rvs[i] = unur_sample_cont(gen);

        end = FLA_Clock();

        time_taken += end - start;

        ++number;
    } while (time_taken < 0.2);

    number = (size_t)pow(10, ceil(log10((double)number)));

    double time_taken_array[7];

    for (size_t rep = 0 ; rep < 7 ; ++rep) {
        start = FLA_Clock();
        for (int n = 0 ; n < number ; ++n) {
            for (size_t i = 0 ; i < num_samples ; ++i) rvs[i] = unur_sample_cont(gen);
        }
        end = FLA_Clock();
        time_taken_array[rep] = (end - start) / number;
    }

    double mean = fmean(time_taken_array, 7);
    double stddev = fstddev(time_taken_array, 7);

    const double scaling[] = {1, 1e3, 1e6, 1e9};

    const char *units[] = {"s", "ms", "µs", "ns"};

    // order = min(-int(math.floor(math.log10(timespan)) // 3), 3)
    int mean_order = 0;
    while (mean * scaling[mean_order] < 1) ++mean_order;
    if (mean_order > 3) mean_order = 3;

    fprintf(stdout, "%.3g %s", mean * scaling[mean_order], units[mean_order]);

    fprintf(stdout, " ± ");

    int stddev_order = 0;
    while (stddev * scaling[stddev_order] < 1) ++stddev_order;
    if (stddev_order > 3) stddev_order = 3;

    fprintf(stdout, "%.3g %s", stddev * scaling[stddev_order], units[stddev_order]);

    fprintf(stdout, " per loop (mean ± std. dev. of %d runs, %lu loop%c each)\n",
                    7, number, (number <= 1 ? '\0' : 's'));

    fflush(stdout);

    unur_free(gen);
    unur_urng_free(urng);
    free(rvs);
}


int main()
{
    UNUR_METHOD_SIGNATURE methods[] = {unur_pinv_new, unur_hinv_new, unur_tdr_new, unur_srou_new};

    const char *method_names[] = {"PINV", "HINV", "TDR", "SROU"};

    FILE *_error_log = fopen("_unuran_errors.log", "w");
    if (!_error_log) {
        fprintf(stderr, "fatal error: can't create/open file '_unuran_errors.log' to log UNU.RAN errors\n");
        exit(EXIT_FAILURE);
    }
    unur_set_stream(_error_log);

    for (int i=0 ; i<4 ; ++i) {
        cont_dist distribution = {
            .pdf = stdnormal_pdf,
            .cdf = stdnormal_cdf,
            .dpdf = stdnormal_dpdf,
            .mode = stdnormal_mode,
            .support = stdnormal_support,
            .params = NULL,
            .num_params = 0
        };
        printf("%s, stdnormal(), [setup]    : ", method_names[i]);
        setup_timeit(distribution, methods[i]);
        printf("%s, stdnormal() [sampling]  : ", method_names[i]);
        sampling_timeit(distribution, methods[i]);
    }

    double gamma_params[] = {0.05, 0.5, 3.0};

    for (int i=0 ; i<sizeof(gamma_params)/sizeof(gamma_params[0]) ; ++i) {
        cont_dist distribution = {
            .pdf = gamma_pdf,
            .cdf = gamma_cdf,
            .dpdf = gamma_dpdf,
            .mode = gamma_mode,
            .support = gamma_support,
            .params = gamma_params + i,
            .num_params = 1
        };
        for (int j=0 ; j<4 ; ++j) {
            printf("%s, gamma(%lf), [setup]    : ", method_names[j], gamma_params[i]);
            setup_timeit(distribution, methods[j]);
            printf("%s, gamma(%lf) [sampling]  : ", method_names[j], gamma_params[i]);
            sampling_timeit(distribution, methods[j]);
        }
    }

    double beta_params[] = {
        0.5, 0.5,
        0.5, 1.0,
        1.3, 1.2,
        3.0, 2.0
    };

    for (int i=0 ; i<sizeof(beta_params)/sizeof(beta_params[0]) ; i+=2) {
        cont_dist distribution = {
            .pdf = beta_pdf,
            .cdf = beta_cdf,
            .dpdf = beta_dpdf,
            .mode = beta_mode,
            .support = beta_support,
            .params = beta_params + i,
            .num_params = 2
        };
        for (int j=0 ; j<4 ; ++j) {
            printf("%s, beta(%lf, %lf), [setup]    : ", method_names[j], beta_params[i], beta_params[i+1]);
            setup_timeit(distribution, methods[j]);
            printf("%s, beta(%lf, %lf) [sampling]  : ", method_names[j], beta_params[i], beta_params[i+1]);
            sampling_timeit(distribution, methods[j]);
        }
    }

    /* Non-standard distribution ==> generalized-normal benchmark */

    double gennorm_params[] = {0.25, 0.45, 0.75, 1., 1.5, 2, 5, 8};

    for (int i=0 ; i<sizeof(gennorm_params)/sizeof(gennorm_params[0]) ; ++i) {
        cont_dist distribution = {
            .pdf = gennorm_pdf,
            .cdf = gennorm_cdf,
            .dpdf = gennorm_dpdf,
            .mode = gennorm_mode,
            .support = gennorm_support,
            .params = gennorm_params + i,
            .num_params = 1
        };
        printf("%s, gennorm(%lf), [setup]    : ", "PINV", gennorm_params[i]);
        setup_timeit(distribution, unur_pinv_new);
        printf("%s, gennorm(%lf) [sampling]  : ", "PINV", gennorm_params[i]);
        sampling_timeit(distribution, unur_pinv_new);
    }

    fclose(_error_log);
}
