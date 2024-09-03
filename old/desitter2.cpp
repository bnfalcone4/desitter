#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <cstdlib>

// Constants
const double Lzero = 0.5;
const double cf = 1 / (4 * M_PI * std::pow(Lzero, 2));

// Functions
inline double switch_func(double z, double DH) {
    return std::pow(std::cos(M_PI * z / (2 * DH)), 4);
}

inline double k(int n, int l, int r) {
    return std::sqrt(n * n + l * l + r * r);
}

inline double nonunit(double z, double k_val) {
    return std::exp(-z) / std::sqrt(k_val);
}

inline double RealDyn(double z, double EA, double k_val, double exp_val) {
    double phase = 2 * M_PI * k_val * exp_val / Lzero;
    return std::cos(z * EA) * std::cos(phase) + std::sin(z * EA) * std::sin(phase);
}

inline double ImDyn(double z, double EA, double k_val, double exp_val) {
    double phase = 2 * M_PI * k_val * exp_val / Lzero;
    return std::cos(z * EA) * std::sin(phase) - std::sin(z * EA) * std::cos(phase);
}

double fReal(double z, double DH, double EA, double k_val, double exp_val) {
    return switch_func(z, DH) * nonunit(z, k_val) * RealDyn(z, EA, k_val, exp_val);
}

double fImag(double z, double DH, double EA, double k_val, double exp_val) {
    return switch_func(z, DH) * nonunit(z, k_val) * ImDyn(z, EA, k_val, exp_val);
}

// Trapezoidal integration
inline double integrate(const std::vector<double>& y, const std::vector<double>& x) {
    double sum = 0.0;
    for (size_t i = 1; i < x.size(); ++i) {
        sum += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    }
    return sum;
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <N> <L> <R> <n_ea> <n_dh>" << std::endl;
        return 1;
    }

    int N = std::stoi(argv[1]);
    int L = std::stoi(argv[2]);
    int R = std::stoi(argv[3]);
    int n_ea = std::stoi(argv[4]);
    int n_dh = std::stoi(argv[5]);

    // Define parameter grids
    std::vector<double> EA_array(n_ea);
    std::vector<double> DH_array(n_dh);
    for (int i = 0; i < n_ea; ++i) {
        EA_array[i] = -12 + (24.0 / (n_ea - 1)) * i;
    }
    for (int i = 0; i < n_dh; ++i) {
        DH_array[i] = 0.175 + (0.625 / (n_dh - 1)) * i;
    }

    std::vector<std::vector<double> > fR_array(n_ea, std::vector<double>(n_dh, 0.0));
    std::vector<std::vector<double> > fI_array(n_ea, std::vector<double>(n_dh, 0.0));
    std::vector<std::vector<std::complex<double> > > ResultadoIntegral(n_ea, std::vector<std::complex<double> >(n_dh, std::complex<double>(0.0, 0.0)));
    std::vector<std::vector<double> > F_array(n_ea, std::vector<double>(n_dh, 0.0));

    std::vector<double> z_list(100000);
    for (int i = 0; i < 100000; ++i) {
        z_list[i] = -DH_array[0] + (2.0 * DH_array[0] / 99999) * i;  // Precompute z_list assuming DH_array[0] as a base
    }

    for (int ind_ea = 0; ind_ea < n_ea; ++ind_ea) {
        double ea = EA_array[ind_ea];
        for (int ind_dh = 0; ind_dh < n_dh; ++ind_dh) {
            double dh = DH_array[ind_dh];

            // Recalculate z_list for the current dh
            for (int i = 0; i < 100000; ++i) {
                z_list[i] = -dh + (2.0 * dh / 99999) * i;
            }

            for (int n = 0; n <= N; ++n) {
                for (int l = 0; l <= L; ++l) {
                    for (int r = 0; r <= R; ++r) {
                        if (n == 0 && l == 0 && r == 0) {
                            continue;
                        }

                        double k_val = k(n, l, r);
                        double exp_val = std::exp(-dh);  // Assuming exp(-z) where z is dh for simplicity

                        std::vector<double> fR_values(100000);
                        std::vector<double> fI_values(100000);

                        for (int i = 0; i < 100000; ++i) {
                            fR_values[i] = fReal(z_list[i], dh, ea, k_val, exp_val);
                            fI_values[i] = fImag(z_list[i], dh, ea, k_val, exp_val);
                        }

                        double integrated_fR = integrate(fR_values, z_list);
                        double integrated_fI = integrate(fI_values, z_list);

                        fR_array[ind_ea][ind_dh] = integrated_fR;
                        fI_array[ind_ea][ind_dh] = integrated_fI;

                        ResultadoIntegral[ind_ea][ind_dh] += std::complex<double>(integrated_fR, integrated_fI);
                        F_array[ind_ea][ind_dh] += std::norm(ResultadoIntegral[ind_ea][ind_dh]) * cf;
                    }
                }
            }
        }
    }

    // Save result to a file
    std::ofstream outfile("result.txt");
    for (int i = 0; i < n_ea; ++i) {
        for (int j = 0; j < n_dh; ++j) {
            outfile << F_array[i][j] << " ";
        }
        outfile << std::endl;
    }
    outfile.close();

    std::cout << "Results saved to result.txt" << std::endl;

    return 0;
}
