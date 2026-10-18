#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>



static size_t g_func_eval_count = 0;

void reset_func_counter() {
    g_func_eval_count = 0;
}

size_t get_func_counter() {
    return g_func_eval_count;
}

double P(double x) {
    return std::pow(x, 5) - 6.2 * std::pow(x, 3) + 3.5 * std::pow(x, 2) - 7.0 * x;
}

double func(double x) {
    ++g_func_eval_count; // Счетчик встроен непосредственно в код функции
    return P(x) * std::cos(2.0 * x);
}

double exact_integral() {
    // I_exact = 348.25 * sin(20) + 35 * cos(20)
    return 348.25 * std::sin(20.0) + 35.0 * std::cos(20.0);
}


struct MethodResult {
    std::string method_name;
    double eps;
    double integral;
    double exact_integral;
    double error;
    size_t func_calls;
    size_t segments;
};

struct StepResult {
    double integral;
    size_t segments;
};


StepResult trapezoidal_adaptive_step(double a, double b, double fa, double fb,
                                     double I_whole, double eps, int depth, int max_depth) {
    double mid = 0.5 * (a + b);
    double fmid = func(mid);

    double h_half = 0.5 * (b - a);
    double I_left = 0.5 * h_half * (fa + fmid);
    double I_right = 0.5 * h_half * (fmid + fb);
    double I_half = I_left + I_right;

    // Оценка погрешности по правилу Рунге для метода 2-го порядка:
    // delta = |I_half - I_whole| / (2^2 - 1) = |I_half - I_whole| / 3
    double delta = std::abs(I_half - I_whole) / 3.0;

    if (delta <= eps || depth >= max_depth) {
        return {I_half, 1};
    }

    auto left_res = trapezoidal_adaptive_step(a, mid, fa, fmid, I_left, eps * 0.5, depth + 1, max_depth);
    auto right_res = trapezoidal_adaptive_step(mid, b, fmid, fb, I_right, eps * 0.5, depth + 1, max_depth);

    return {left_res.integral + right_res.integral, left_res.segments + right_res.segments};
}

MethodResult integrate_trapezoidal_adaptive(double a, double b, double eps) {
    reset_func_counter();

    double fa = func(a);
    double fb = func(b);
    double I_whole = 0.5 * (b - a) * (fa + fb);

    int max_depth = 28;
    StepResult res = trapezoidal_adaptive_step(a, b, fa, fb, I_whole, eps, 0, max_depth);

    size_t calls = get_func_counter();
    double I_exact = exact_integral();
    double err = std::abs(res.integral - I_exact);

    return {"Trapezoidal", eps, res.integral, I_exact, err, calls, res.segments};
}

// Метод 2: Квадратурные формулы Гаусса–Лобатто 4 и 5 узлов с адаптивным разбиением

double transform_canonical(double a, double b, double t) {
    return 0.5 * (b + a) + 0.5 * (b - a) * t;
}

double lobatto_4(double a, double b, double fa, double fb) {
    static const double t[2] = {-std::sqrt(1.0 / 5.0), std::sqrt(1.0 / 5.0)};
    static const double w[2] = {5.0 / 6.0, 5.0 / 6.0};
    static const double w_end = 1.0 / 6.0;

    double sum = w_end * (fa + fb);
    for (int i = 0; i < 2; ++i) {
        double x = transform_canonical(a, b, t[i]);
        sum += w[i] * func(x);
    }
    return 0.5 * (b - a) * sum;
}

double lobatto_5(double a, double b, double fa, double fb) {
    static const double t[3] = {-std::sqrt(3.0 / 7.0), 0.0, std::sqrt(3.0 / 7.0)};
    static const double w[3] = {49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0};
    static const double w_end = 1.0 / 10.0;

    double sum = w_end * (fa + fb);
    for (int i = 0; i < 3; ++i) {
        double x = transform_canonical(a, b, t[i]);
        sum += w[i] * func(x);
    }
    return 0.5 * (b - a) * sum;
}

StepResult lobatto_adaptive_step(double a, double b, double fa, double fb,
                                 double eps, int depth, int max_depth) {
    double I4 = lobatto_4(a, b, fa, fb);
    double I5 = lobatto_5(a, b, fa, fb);

    double delta = std::abs(I5 - I4);

    if (delta <= eps || depth >= max_depth) {
        return {I5, 1};
    }

    double mid = 0.5 * (a + b);
    double fmid = func(mid);

    auto left_res = lobatto_adaptive_step(a, mid, fa, fmid, eps * 0.5, depth + 1, max_depth);
    auto right_res = lobatto_adaptive_step(mid, b, fmid, fb, eps * 0.5, depth + 1, max_depth);

    return {left_res.integral + right_res.integral, left_res.segments + right_res.segments};
}

MethodResult integrate_lobatto_adaptive(double a, double b, double eps) {
    reset_func_counter();

    double fa = func(a);
    double fb = func(b);

    int max_depth = 25;
    StepResult res = lobatto_adaptive_step(a, b, fa, fb, eps, 0, max_depth);

    size_t calls = get_func_counter();
    double I_exact = exact_integral();
    double err = std::abs(res.integral - I_exact);

    return {"Gauss-Lobatto", eps, res.integral, I_exact, err, calls, res.segments};
}


void export_results_csv(const std::string &filepath,
                        const std::vector<MethodResult> &trapezoid_results,
                        const std::vector<MethodResult> &lobatto_results) {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        throw std::runtime_error("Не удалось открыть файл для записи: " + filepath);
    }

    out << "method,eps,integral,exact_integral,error,func_calls,segments\n";

    auto write_row = [&out](const MethodResult &r) {
        out << r.method_name << ","
            << std::scientific << std::setprecision(6) << r.eps << ","
            << std::fixed << std::setprecision(12) << r.integral << ","
            << std::fixed << std::setprecision(12) << r.exact_integral << ","
            << std::scientific << std::setprecision(8) << r.error << ","
            << r.func_calls << ","
            << r.segments << "\n";
    };

    for (const auto &r : trapezoid_results) {
        write_row(r);
    }
    for (const auto &r : lobatto_results) {
        write_row(r);
    }

    out.close();
}


int main() {
    try {
        const double a = -10.0;
        const double b = 10.0;
        const double I_exact = exact_integral();

        std::vector<double> eps_values;
        for (double p = 1.0; p <= 10.0; p += 0.5) {
            eps_values.push_back(std::pow(10.0, -p));
        }

        std::vector<MethodResult> trapezoid_results;
        std::vector<MethodResult> lobatto_results;

        std::cout << std::left << std::setw(12) << "Eps"
                  << std::setw(16) << "Trap Calls"
                  << std::setw(16) << "Trap Error"
                  << std::setw(16) << "Lobatto Calls"
                  << std::setw(16) << "Lobatto Error"
                  << "\n";
        std::cout << std::string(76, '-') << "\n";

        for (double eps : eps_values) {
            MethodResult tr = integrate_trapezoidal_adaptive(a, b, eps);
            MethodResult lob = integrate_lobatto_adaptive(a, b, eps);

            trapezoid_results.push_back(tr);
            lobatto_results.push_back(lob);

            std::cout << std::scientific << std::setprecision(2) << std::setw(12) << eps
                      << std::setw(16) << tr.func_calls
                      << std::scientific << std::setprecision(3) << std::setw(16) << tr.error
                      << std::setw(16) << lob.func_calls
                      << std::scientific << std::setprecision(3) << std::setw(16) << lob.error
                      << "\n";
        }

        std::string csv_path = "../data/comparison_data.csv";
        export_results_csv(csv_path, trapezoid_results, lobatto_results);

    } catch (const std::exception &ex) {
        std::cerr << "ОШИБКА: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
