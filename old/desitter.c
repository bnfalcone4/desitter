#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

// Constants
const double Lzero = 0.5;
double cf;  // cf will be calculated at runtime

// Functions
inline double switch_func(double z, double DH) {
    return pow(cos(M_PI * z / (2 * DH)), 4);
}

inline double k(int n, int l, int r) {
    return sqrt(n * n + l * l + r * r);
}

inline double nonunit(double z, double k_val) {
    return exp(-z) / sqrt(k_val);
}

inline double RealDyn(double z, double EA, double k_val, double exp_val) {
    double phase = 2 * M_PI * k_val * exp_val / Lzero;
    return cos(z * EA) * cos(phase) + sin(z * EA) * sin(phase);
}

inline double ImDyn(double z, double EA, double k_val, double exp_val) {
    double phase = 2 * M_PI * k_val * exp_val / Lzero;
    return cos(z * EA) * sin(phase) - sin(z * EA) * cos(phase);
}

double fReal(double z, double DH, double EA, double k_val, double exp_val) {
    return switch_func(z, DH) * nonunit(z, k_val) * RealDyn(z, EA, k_val, exp_val);
}

double fImag(double z, double DH, double EA, double k_val, double exp_val) {
    return switch_func(z, DH) * nonunit(z, k_val) * ImDyn(z, EA, k_val, exp_val);
}

// Trapezoidal integration
double integrate(const double* y, const double* x, int size) {
    double sum = 0.0;
    for (int i = 1; i < size; ++i) {
        sum += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    }
    return sum;
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        printf("Usage: %s <N> <L> <R> <n_ea> <n_dh>\n", argv[0]);
        return 1;
    }

    // Calculate cf at runtime
    cf = 1.0 / (4.0 * M_PI * pow(Lzero, 2));

    int N = atoi(argv[1]);
    int L = atoi(argv[2]);
    int R = atoi(argv[3]);
    int n_ea = atoi(argv[4]);
    int n_dh = atoi(argv[5]);

    // Define parameter grids
    double *EA_array = (double*) malloc(n_ea * sizeof(double));
    double *DH_array = (double*) malloc(n_dh * sizeof(double));

    for (int i = 0; i < n_ea; ++i) {
        EA_array[i] = -12 + (24.0 / (n_ea - 1)) * i;
    }
    for (int i = 0; i < n_dh; ++i) {
        DH_array[i] = 0.175 + (0.625 / (n_dh - 1)) * i;
    }

    double **fR_array = (double**) malloc(n_ea * sizeof(double*));
    double **fI_array = (double**) malloc(n_ea * sizeof(double*));
    double complex **ResultadoIntegral = (double complex**) malloc(n_ea * sizeof(double complex*));
    double **F_array = (double**) malloc(n_ea * sizeof(double*));

    for (int i = 0; i < n_ea; ++i) {
        fR_array[i] = (double*) malloc(n_dh * sizeof(double));
        fI_array[i] = (double*) malloc(n_dh * sizeof(double));
        ResultadoIntegral[i] = (double complex*) malloc(n_dh * sizeof(double complex));
        F_array[i] = (double*) malloc(n_dh * sizeof(double));
    }

    double *z_list = (double*) malloc(100000 * sizeof(double));

    for (int ind_ea = 0; ind_ea < n_ea; ++ind_ea) {
        double ea = EA_array[ind_ea];
        for (int ind_dh = 0; ind_dh < n_dh; ++ind_dh) {
            double dh = DH_array[ind_dh];

            // Recalculate z_list for the current dh
            for (int i = 0; i < 100000; ++i) {
                z_list[i] = -dh + (2.0 * dh / 99999) * i;
            }

            for (int n = 0; n <= N; n++) {
                for (int l = 0; l <= L; l++) {
                    for (int r = 0; r <= R; r++) {
                        if (n == 0 && l == 0 && r == 0) {
                            continue;
                        }

                        double k_val = k(n, l, r);
                        double exp_val = exp(-dh);  // Assuming exp(-z) where z is dh for simplicity

                        double *fR_values = (double*) malloc(100000 * sizeof(double));
                        double *fI_values = (double*) malloc(100000 * sizeof(double));

                        for (int i = 0; i < 100000; ++i) {
                            fR_values[i] = fReal(z_list[i], dh, ea, k_val, exp_val);
                            fI_values[i] = fImag(z_list[i], dh, ea, k_val, exp_val);
                        }

                        double integrated_fR = integrate(fR_values, z_list, 100000);
                        double integrated_fI = integrate(fI_values, z_list, 100000);

                        fR_array[ind_ea][ind_dh] = integrated_fR;
                        fI_array[ind_ea][ind_dh] = integrated_fI;

                        ResultadoIntegral[ind_ea][ind_dh] += integrated_fR + integrated_fI * I;
                        F_array[ind_ea][ind_dh] += pow(cabs(ResultadoIntegral[ind_ea][ind_dh]), 2) * cf;

                        free(fR_values);
                        free(fI_values);
                    }
                }
            }
        }
    }

    // Save result to a file
    FILE* outfile = fopen("result.txt", "w");
    for (int i = 0; i < n_ea; ++i) {
        for (int j = 0; j < n_dh; ++j) {
            fprintf(outfile, "%f ", F_array[i][j]);
        }
        fprintf(outfile, "\n");
    }
    fclose(outfile);

    printf("Results saved to result.txt\n");

    // Free allocated memory
    free(EA_array);
    free(DH_array);
    for (int i = 0; i < n_ea; ++i) {
        free(fR_array[i]);
        free(fI_array[i]);
        free(ResultadoIntegral[i]);
        free(F_array[i]);
    }
    free(fR_array);
    free(fI_array);
    free(ResultadoIntegral);
    free(F_array);
    free(z_list);

    return 0;
}
