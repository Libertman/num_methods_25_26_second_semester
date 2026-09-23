#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct Point {
  double x;
  double y;
};

struct ErrorPoint {
  double x;
  double y_exact;
  double error;
};

struct AccuracyResult {
  double eps;
  double max_error;
  size_t steps;
};

struct StepPoint {
  double x;
  double h;
};

struct SweepInterval {
  double x_start;
  double h;
  double alpha_start, beta_start;
  double alpha_mid, beta_mid;
  double alpha_end, beta_end;
};

// p(x) y'' + q(x) y' + r(x) y = f(x)
// p(x) = 1

double q_func(double x) { return 1.0 + std::sin(x) * std::sin(x); }

double r_func(double x) { return std::cos(x) * std::cos(x); }

double f_func(double x) { return 3.0 * std::exp(x); }

double exact_solution(double x) { return std::exp(x); }

void forward_rhs(double x, double alpha, double beta, double &dalpha,
                 double &dbeta) {
  double q = q_func(x);
  double r = r_func(x);
  double f = f_func(x);

  dalpha = -alpha * alpha - q * alpha - r;
  dbeta = -(alpha + q) * beta + f;
}

void rk4_forward_step(double x, double h, double alpha_in, double beta_in,
                      double &alpha_out, double &beta_out) {
  double k1_a, k1_b;
  forward_rhs(x, alpha_in, beta_in, k1_a, k1_b);

  double k2_a, k2_b;
  forward_rhs(x + 0.5 * h, alpha_in + 0.5 * h * k1_a, beta_in + 0.5 * h * k1_b,
              k2_a, k2_b);

  double k3_a, k3_b;
  forward_rhs(x + 0.5 * h, alpha_in + 0.5 * h * k2_a, beta_in + 0.5 * h * k2_b,
              k3_a, k3_b);

  double k4_a, k4_b;
  forward_rhs(x + h, alpha_in + h * k3_a, beta_in + h * k3_b, k4_a, k4_b);

  alpha_out = alpha_in + (h / 6.0) * (k1_a + 2.0 * k2_a + 2.0 * k3_a + k4_a);
  beta_out = beta_in + (h / 6.0) * (k1_b + 2.0 * k2_b + 2.0 * k3_b + k4_b);
}

double rk4_backward_step(double h, double y_curr, double a_start,
                         double b_start, double a_mid, double b_mid,
                         double a_end, double b_end) {
  double k1 = a_end * y_curr + b_end;

  double k2 = a_mid * (y_curr - 0.5 * h * k1) + b_mid;

  double k3 = a_mid * (y_curr - 0.5 * h * k2) + b_mid;

  double k4 = a_start * (y_curr - h * k3) + b_start;

  return y_curr - (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

std::vector<Point> solve_factorization_fixed(double a, double b, double h,
                                             double alpha_0, double alpha_1,
                                             double A, double beta_0,
                                             double beta_1, double B) {
  int N = static_cast<int>(std::round((b - a) / h));
  if (N < 1)
    N = 1;
  double actual_h = (b - a) / N;

  std::vector<SweepInterval> intervals(N);

  double alpha_curr = -alpha_0 / alpha_1;
  double beta_curr = A / alpha_1;

  double x_curr = a;

  for (int i = 0; i < N; ++i) {
    double x_next = x_curr + actual_h;

    double alpha_mid, beta_mid;
    rk4_forward_step(x_curr, 0.5 * actual_h, alpha_curr, beta_curr, alpha_mid,
                     beta_mid);

    double alpha_next, beta_next;
    rk4_forward_step(x_curr + 0.5 * actual_h, 0.5 * actual_h, alpha_mid,
                     beta_mid, alpha_next, beta_next);

    intervals[i] = {x_curr,    actual_h, alpha_curr, beta_curr,
                    alpha_mid, beta_mid, alpha_next, beta_next};

    x_curr = x_next;
    alpha_curr = alpha_next;
    beta_curr = beta_next;
  }

  double denom = beta_0 + beta_1 * alpha_curr;
  if (std::abs(denom) < 1e-14) {
    throw std::runtime_error(
        "Знаменатель в краевом условии на правом конце равен 0");
  }
  double y_right = (B - beta_1 * beta_curr) / denom;

  std::vector<Point> solution(N + 1);
  solution[N] = {b, y_right};

  double y_curr = y_right;
  for (int i = N - 1; i >= 0; --i) {
    const auto &iv = intervals[i];
    y_curr =
        rk4_backward_step(iv.h, y_curr, iv.alpha_start, iv.beta_start,
                          iv.alpha_mid, iv.beta_mid, iv.alpha_end, iv.beta_end);
    solution[i] = {iv.x_start, y_curr};
  }

  return solution;
}

std::vector<Point>
solve_factorization_adaptive(double a, double b, double eps, double alpha_0,
                             double alpha_1, double A, double beta_0,
                             double beta_1, double B,
                             std::vector<StepPoint> &step_history) {
  step_history.clear();

  double x = a;
  double alpha_curr = -alpha_0 / alpha_1;
  double beta_curr = A / alpha_1;

  double h = std::min(0.5, std::max(0.001, std::pow(eps, 0.2)));

  std::vector<SweepInterval> intervals;
  const double safety = 0.9;
  const double min_h = 1e-7;
  const double max_h = 0.5;

  while (x < b - 1e-12) {
    if (x + h > b) {
      h = b - x;
    }

    double alpha_full, beta_full;
    rk4_forward_step(x, h, alpha_curr, beta_curr, alpha_full, beta_full);

    double alpha_mid, beta_mid;
    rk4_forward_step(x, 0.5 * h, alpha_curr, beta_curr, alpha_mid, beta_mid);

    double alpha_half, beta_half;
    rk4_forward_step(x + 0.5 * h, 0.5 * h, alpha_mid, beta_mid, alpha_half,
                     beta_half);

    double err_alpha = std::abs(alpha_half - alpha_full) / 15.0;
    double err_beta = std::abs(beta_half - beta_full) / 15.0;
    double err = std::max(err_alpha, err_beta);

    if (err <= eps || h <= min_h) {
      step_history.push_back({x, h});

      intervals.push_back({x, h, alpha_curr, beta_curr, alpha_mid, beta_mid,
                           alpha_half, beta_half});

      x += h;
      alpha_curr = alpha_half;
      beta_curr = beta_half;

      double factor = safety * std::pow(eps / (err + 1e-15), 0.2);
      factor = std::max(0.5, std::min(2.0, factor));
      h = std::min(max_h, h * factor);
    } else {
      double factor = safety * std::pow(eps / (err + 1e-15), 0.2);
      factor = std::max(0.1, std::min(0.5, factor));
      h = std::max(min_h, h * factor);
    }
  }

  double denom = beta_0 + beta_1 * alpha_curr;
  if (std::abs(denom) < 1e-14) {
    throw std::runtime_error(
        "Знаменатель в краевом условии на правом конце равен 0");
  }
  double y_right = (B - beta_1 * beta_curr) / denom;

  size_t N = intervals.size();
  std::vector<Point> solution(N + 1);
  solution[N] = {b, y_right};

  double y_curr = y_right;
  for (int i = static_cast<int>(N) - 1; i >= 0; --i) {
    const auto &iv = intervals[i];
    y_curr =
        rk4_backward_step(iv.h, y_curr, iv.alpha_start, iv.beta_start,
                          iv.alpha_mid, iv.beta_mid, iv.alpha_end, iv.beta_end);
    solution[i] = {iv.x_start, y_curr};
  }

  return solution;
}

std::vector<ErrorPoint> compute_errors(const std::vector<Point> &solution) {
  std::vector<ErrorPoint> errors;
  errors.reserve(solution.size());

  for (const auto &p : solution) {
    double exact = exact_solution(p.x);
    errors.push_back({p.x, exact, std::abs(p.y - exact)});
  }
  return errors;
}

double compute_max_error(const std::vector<Point> &solution) {
  double max_err = 0.0;
  for (const auto &p : solution) {
    double err = std::abs(p.y - exact_solution(p.x));
    if (err > max_err) {
      max_err = err;
    }
  }
  return max_err;
}

std::vector<AccuracyResult> accuracy_study(double a, double b, double alpha_0,
                                           double alpha_1, double A,
                                           double beta_0, double beta_1,
                                           double B) {
  std::vector<AccuracyResult> results;
  std::vector<double> eps_values = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};

  std::vector<StepPoint> unused_history;

  for (double eps : eps_values) {
    auto sol = solve_factorization_adaptive(a, b, eps, alpha_0, alpha_1, A,
                                            beta_0, beta_1, B, unused_history);

    double max_err = compute_max_error(sol);
    results.push_back({eps, max_err, sol.size() - 1});
  }

  return results;
}

// writes

void write_solution_csv(const std::string &filename,
                        const std::vector<Point> &sol_h1,
                        const std::vector<Point> &sol_h2, double a, double b) {
  std::ofstream file(filename);
  if (!file.is_open())
    throw std::runtime_error("Не удалось открыть файл: " + filename);

  file << "type,x,y_numeric,y_exact\n";

  for (const auto &p : sol_h1) {
    file << "h1," << p.x << "," << p.y << "," << exact_solution(p.x) << "\n";
  }

  for (const auto &p : sol_h2) {
    file << "h2," << p.x << "," << p.y << "," << exact_solution(p.x) << "\n";
  }
}

void write_error_csv(const std::string &filename,
                     const std::vector<ErrorPoint> &err_h1,
                     const std::vector<ErrorPoint> &err_h2) {
  std::ofstream file(filename);
  if (!file.is_open())
    throw std::runtime_error("Не удалось открыть файл: " + filename);

  file << "type,x,y_exact,error\n";

  for (const auto &ep : err_h1) {
    file << "h1," << ep.x << "," << ep.y_exact << "," << ep.error << "\n";
  }

  for (const auto &ep : err_h2) {
    file << "h2," << ep.x << "," << ep.y_exact << "," << ep.error << "\n";
  }
}

void write_accuracy_csv(const std::string &filename,
                        const std::vector<AccuracyResult> &data) {
  std::ofstream file(filename);
  if (!file.is_open())
    throw std::runtime_error("Не удалось открыть файл: " + filename);

  file << "eps,max_error,steps\n";

  for (const auto &r : data) {
    file << r.eps << "," << r.max_error << "," << r.steps << "\n";
  }
}

void write_step_size_csv(const std::string &filename,
                         const std::vector<StepPoint> &data) {
  std::ofstream file(filename);
  if (!file.is_open())
    throw std::runtime_error("Не удалось открыть файл: " + filename);

  file << "x,h\n";

  for (const auto &sp : data) {
    file << sp.x << "," << sp.h << "\n";
  }
}

int main() {
  try {
    double a = 0.0;
    double b = 1.0;

    double alpha_0 = 1.0;
    double alpha_1 = 1.0;
    double A = alpha_0 * exact_solution(a) + alpha_1 * exact_solution(a);

    double beta_0 = 1.0;
    double beta_1 = 0.0;
    double B = beta_0 * exact_solution(b);

    double h1 = 0.1;
    double h2 = 0.05;

    std::cout << "Расчет решения при h1 = " << h1 << " и h2 = " << h2
              << "...\n";
    auto sol_h1 = solve_factorization_fixed(a, b, h1, alpha_0, alpha_1, A,
                                            beta_0, beta_1, B);
    auto sol_h2 = solve_factorization_fixed(a, b, h2, alpha_0, alpha_1, A,
                                            beta_0, beta_1, B);

    auto err_h1 = compute_errors(sol_h1);
    auto err_h2 = compute_errors(sol_h2);

    std::cout << "Максимальная ошибка при h1 = " << h1 << ": "
              << compute_max_error(sol_h1) << "\n";
    std::cout << "Максимальная ошибка при h2 = " << h2 << ": "
              << compute_max_error(sol_h2) << "\n";

    std::cout << "Исследование точности от заданной погрешности eps\n";
    auto accuracy =
        accuracy_study(a, b, alpha_0, alpha_1, A, beta_0, beta_1, B);

    double eps_step = 1e-6;
    std::cout << "Расчет адаптивного шага при eps = " << eps_step << "...\n";
    std::vector<StepPoint> step_history;
    auto sol_adaptive = solve_factorization_adaptive(
        a, b, eps_step, alpha_0, alpha_1, A, beta_0, beta_1, B, step_history);

    std::cout << "Число шагов при eps = 1e-6: " << step_history.size() << "\n";
    std::cout << "Максимальная ошибка адаптивного решения: "
              << compute_max_error(sol_adaptive) << "\n";

    std::string data_dir = "../../data";

    write_solution_csv(data_dir + "solution.csv", sol_h1, sol_h2, a, b);
    write_error_csv(data_dir + "error.csv", err_h1, err_h2);
    write_accuracy_csv(data_dir + "accuracy.csv", accuracy);
    write_step_size_csv(data_dir + "step_size.csv", step_history);

  } catch (const std::exception &ex) {
    std::cerr << "ОШИБКА: " << ex.what() << std::endl;
    return 1;
  }

  return 0;
}
