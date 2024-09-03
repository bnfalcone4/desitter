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
double switch_func(double z, double DH) {
    return std::pow(std::cos(M_PI * z / (2 * DH)), 4);
}

double k(int n, int l, int r) {
    return std::sqrt(n * n + l * l + r * r);
}

double nonunit(double z, int n, int l, int r) {
    return std::exp(-z) / std::sqrt(k(n, l, r));
}

double RealDyn(double z, double EA, int n, int l, int r) {
    double k_val = k(n, l, r);
    double exp_val = std::exp(-z);
    return std::cos(z * EA) * std::cos(2 * M_PI * k_val * exp_val / Lzero) +
           std::sin(z * EA) * std::sin(2 * M_PI * k_val * exp_val / Lzero);
}

double ImDyn(double z, double EA, int n, int l, int r) {
    double k_val = k(n, l, r);
    double exp_val = std::exp(-z);
    return std::cos(z * EA) * std::sin(2 * M_PI * k_val * exp_val / Lzero) -
           std::sin(z * EA) * std::cos(2 * M_PI * k_val * exp_val / Lzero);
}

double fReal(double z, double DH, double EA, int n, int l, int r) {
    return switch_func(z, DH) * nonunit(z, n, l, r) * RealDyn(z, EA, n, l, r);
}

double fImag(double z, double DH, double EA, int n, int l, int r) {
    return switch_func(z, DH) * nonunit(z, n, l, r) * ImDyn(z, EA, n, l, r);
}

// Trapezoidal integration
double integrate(const std::vector<double>& y, const std::vector<double>& x) {
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

    // Fix the space between consecutive right angle brackets and correct the initializer list
    std::vector<std::vector<double> > fR_array(n_ea, std::vector<double>(n_dh, 0.0));
    std::vector<std::vector<double> > fI_array(n_ea, std::vector<double>(n_dh, 0.0));
    std::vector<std::vector<std::complex<double> > > ResultadoIntegral(n_ea, std::vector<std::complex<double> >(n_dh, std::complex<double>(0.0, 0.0)));
    std::vector<std::vector<double> > F_array(n_ea, std::vector<double>(n_dh, 0.0));

    for (int ind_ea = 0; ind_ea < n_ea; ++ind_ea) {
        double ea = EA_array[ind_ea];
        for (int ind_dh = 0; ind_dh < n_dh; ++ind_dh) {
            double dh = DH_array[ind_dh];
            std::vector<double> z_list(100000);
            for (int i = 0; i < 100000; ++i) {
                z_list[i] = -dh + (2.0 * dh / 99999) * i;
            }

            for (int n = 0; n <= N; ++n) {
                for (int l = 0; l <= L; ++l) {
                    for (int r = 0; r <= R; ++r) {
                        if (n == 0 && l == 0 && r == 0) {
                            continue;
                        }

                        std::vector<double> fR_values(100000);
                        std::vector<double> fI_values(100000);

                        for (int i = 0; i < 100000; ++i) {
                            fR_values[i] = fReal(z_list[i], ea, dh, n, l, r);
                            fI_values[i] = fImag(z_list[i], ea, dh, n, l, r);
                        }

                        fR_array[ind_ea][ind_dh] = integrate(fR_values, z_list);
                        fI_array[ind_ea][ind_dh] = integrate(fI_values, z_list);

                        ResultadoIntegral[ind_ea][ind_dh] += std::complex<double>(fR_array[ind_ea][ind_dh], fI_array[ind_ea][ind_dh]);
                        F_array[ind_ea][ind_dh] += std::norm(ResultadoIntegral[ind_ea][ind_dh]) * cf;
                    }
                }
            }

            std::cout << ind_ea << std::endl;
        }
    }

    // Save result to a file (as a placeholder for plotting)
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
