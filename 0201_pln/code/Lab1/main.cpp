#include <iostream>
#include <vector>
#include <numbers>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>


// ---------------- EXTRA FUNCTIONS ----------------

double sign(double x)
{
    if (x > 0) return 1.0;
    else if (x == 0) return 0.0;
    else return -1.0;
}

std::vector<double> divided_differences(const std::vector<double>& x, const std::vector<double>& y)
{
    int n = x.size();
    std::vector<double> coef = y;

    for (int j = 1; j < n; j++)
    {
        for (int i = n - 1; i >= j; i--)
        {
            coef[i] = (coef[i] - coef[i - 1]) / (x[i] - x[i - j]);
        }
    }
    return coef;
}


// ---------------- INPUT FUNCTIONS ----------------

double f1(double x)
{
    return log(1 + x * x) + pow(x, 3);
}

double f2(double x)
{
    return pow(x, 5) + 0.4 * sign(x) * pow(x, 4) + 2;
}


// ---------------- GRIDS ----------------

std::vector<double> uniform_grid(double a, double b, int n)
{
    std::vector<double> x(n);
    double h = (b - a) / (n - 1);
    for (int i = 0; i < n; i++)
        x[i] = a + i * h;
    return x;
}

std::vector<double> chebyshev_grid(double a, double b, int n)
{
    std::vector<double> x(n);
    for (int i = 0; i < n; i++)
    {
        double t = cos((2.0 * i + 1) / (2 * n) * std::numbers::pi);
        x[i] = (a + b) / 2 + (b - a) / 2 * t;
    }
    return x;
}


// ---------------- METHOD ----------------

double newton_eval(const std::vector<double>& x, const std::vector<double>& coef, double t)
{
    int n = x.size();
    double result = coef[n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        result = result * (t - x[i]) + coef[i];
    }

    return result;
}


// ---------------- FILES ----------------

void save_data(const std::string& filename,
    const std::vector<double>& x,
    const std::vector<double>& f,
    const std::vector<double>& p)
{
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(10);

    for (size_t i = 0; i < x.size(); i++)
    {
        double err = fabs(f[i] - p[i]);
        file << x[i] << " " << f[i] << " " << p[i] << " " << err << "\n";
    }

    file.close();
}


// ---------------- TESTS ----------------

void run_experiment(double (*func)(double),
    std::string prefix,
    double a, double b,
    std::vector<int> node_counts) {

    int dense_N = 10000;

    std::vector<double> x_dense(dense_N);
    for (int i = 0; i < dense_N; i++)
    {
        x_dense[i] = a + (b - a) * i / (dense_N - 1);
    }

    std::vector<double> f_dense(dense_N);
    for (int i = 0; i < dense_N; i++)
    {
        f_dense[i] = func(x_dense[i]);
    }

    for (int n : node_counts)
    {
        // uniform grid
        
        std::vector<double> x_nodes = uniform_grid(a, b, n);

        std::vector<double> y_nodes(n);
        for (int i = 0; i < n; i++)
            y_nodes[i] = func(x_nodes[i]);

        std::vector<double> coef = divided_differences(x_nodes, y_nodes);

        std::vector<double> p_dense(dense_N);
        for (int i = 0; i < dense_N; i++)
            p_dense[i] = newton_eval(x_nodes, coef, x_dense[i]);

        save_data(prefix + "_uniform_n" + std::to_string(n) + ".txt",
            x_dense, f_dense, p_dense);


        // chebyshev grid

        std::vector<double> x_nodes_c = chebyshev_grid(a, b, n);

        std::vector<double> y_nodes_c(n);
        for (int i = 0; i < n; i++)
            y_nodes_c[i] = func(x_nodes_c[i]);

        std::vector<double> coef_c = divided_differences(x_nodes_c, y_nodes_c);

        std::vector<double> p_dense_c(dense_N);
        for (int i = 0; i < dense_N; i++)
            p_dense_c[i] = newton_eval(x_nodes_c, coef_c, x_dense[i]);

        save_data(prefix + "_cheb_n" + std::to_string(n) + ".txt",
            x_dense, f_dense, p_dense_c);
    }
}

void error_vs_n(double (*func)(double),
    std::string filename,
    double a, double b)
{
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(10);

    int dense_N = 10000;

    std::vector<double> x_dense(dense_N);
    for (int i = 0; i < dense_N; i++)
        x_dense[i] = a + (b - a) * i / (dense_N - 1);

    std::vector<double> f_dense(dense_N);
    for (int i = 0; i < dense_N; i++)
        f_dense[i] = func(x_dense[i]);

    for (int n = 5; n <= 100; n++)
    {
        // uniform grid 

        std::vector<double> x_nodes = uniform_grid(a, b, n);

        std::vector<double> y_nodes(n);
        for (int i = 0; i < n; i++)
            y_nodes[i] = func(x_nodes[i]);

        std::vector<double> coef = divided_differences(x_nodes, y_nodes);

        double max_err_uniform = 0;

        for (int i = 0; i < dense_N; i++) {
            double p = newton_eval(x_nodes, coef, x_dense[i]);
            double err = fabs(f_dense[i] - p);
            if (err > max_err_uniform) max_err_uniform = err;
        }

        // chebyshev grid

        std::vector<double> x_nodes_c = chebyshev_grid(a, b, n);

        std::vector<double> y_nodes_c(n);
        for (int i = 0; i < n; i++)
            y_nodes_c[i] = func(x_nodes_c[i]);

        std::vector<double> coef_c = divided_differences(x_nodes_c, y_nodes_c);

        double max_err_cheb = 0;

        for (int i = 0; i < dense_N; i++) {
            double p = newton_eval(x_nodes_c, coef_c, x_dense[i]);
            double err = fabs(f_dense[i] - p);
            if (err > max_err_cheb) max_err_cheb = err;
        }

        file << n << " " << max_err_uniform << " " << max_err_cheb << "\n";
    }

    file.close();
}



void perturbation_experiment(double (*func)(double),
    std::string filename,
    double a, double b,
    int n,
    bool use_cheb) {

    std::ofstream file(filename);
    file << std::fixed << std::setprecision(10);

    std::mt19937 gen(42);

    int dense_N = 10000;

    std::vector<double> x_dense(dense_N);
    for (int i = 0; i < dense_N; i++)
        x_dense[i] = a + (b - a) * i / (dense_N - 1);

    std::vector<double> f_dense(dense_N);
    for (int i = 0; i < dense_N; i++)
        f_dense[i] = func(x_dense[i]);

    std::vector<double> eps_list = { 0.01, 0.02, 0.03, 0.04, 0.05 };

    for (double eps : eps_list) {

        for (int exp_id = 0; exp_id < 20; exp_id++) {
            std::vector<double> x_nodes = use_cheb ?
                chebyshev_grid(a, b, n) :
                uniform_grid(a, b, n);

            std::vector<double> y_nodes(n);
            std::vector<double> y_perturbed(n);

            double avg_noise = 0;

            for (int i = 0; i < n; i++) {
                y_nodes[i] = func(x_nodes[i]);

                double delta = std::uniform_real_distribution<double>(-eps, eps)(gen);
                avg_noise += fabs(delta);

                y_perturbed[i] = y_nodes[i] * (1 + delta);
            }

            avg_noise /= n;

            std::vector<double> coef = divided_differences(x_nodes, y_perturbed);

            double max_err = 0;
            double max_f = 0;

            for (int i = 0; i < dense_N; i++) {
                double p = newton_eval(x_nodes, coef, x_dense[i]);

                double err = fabs(f_dense[i] - p);
                if (err > max_err) max_err = err;

                if (fabs(f_dense[i]) > max_f) max_f = fabs(f_dense[i]);
            }

            double rel_err = max_err / max_f;

            file << avg_noise << " " << rel_err << "\n";
        }
    }

    file.close();
}


// ---------------- MAIN ----------------

int main() {

    double a = -2.0;
    double b = 2.0;

    std::vector<int> node_counts = { 5, 7, 10 };

    std::string src_prefix = "../../data/";

    run_experiment(f1, src_prefix + "f1", a, b, node_counts);
    error_vs_n(f1, src_prefix + "f1_error.txt", a, b);

    run_experiment(f2, src_prefix + "f2", a, b, node_counts);
    error_vs_n(f2, src_prefix + "f2_error.txt", a, b);

    perturbation_experiment(f1, src_prefix + "f1_perturb_uniform.txt", a, b, 10, false);
    perturbation_experiment(f1, src_prefix + "f1_perturb_cheb.txt", a, b, 10, true);

    perturbation_experiment(f2, src_prefix + "f2_perturb_uniform.txt", a, b, 10, false);
    perturbation_experiment(f2, src_prefix + "f2_perturb_cheb.txt", a, b, 10, true);

    return 0;
}
