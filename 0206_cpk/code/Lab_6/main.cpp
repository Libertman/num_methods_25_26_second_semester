#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <numbers>



struct Point
{
    double x;
    double y;
};

struct ErrorPoint
{
    double x;
    double y_exact;
    double error;
};

struct AccuracyResult
{
    double eps;
    double max_error;
};


double f(double x, double y)
{
    return y / x + x * std::cos(x);
}


double exact_solution(double x)
{
    return x * std::sin(x);
}




std::vector<Point> solve_euler_trapezoid(
    double a,
    double b,
    double h,
    double eps,
    double y0
)
{
    const int MAX_ITER = 100;

    std::vector<Point> result;

    double x = a;
    double y = y0;

    Point first;

    first.x = x;
    first.y = y;

    result.push_back(first);

    while (x + h < b + 1e-12)
    {
        double x_next = x + h;

        double f_left = f(x, y);

        double current = y + h * f_left;

        bool converged = false;

        for (int iter = 0; iter < MAX_ITER; iter++)
        {
            double next = y + h * (
                f_left
                + f(x_next, current)
                ) / 2.0;

            if (std::abs(next - current) < eps)
            {
                current = next;
                converged = true;
                break;
            }

            current = next;
        }

        if (!converged)
        {
            throw std::runtime_error(
                "Corrector iteration did not converge"
            );
        }

        y = current;
        x = x_next;

        Point p;

        p.x = x;
        p.y = y;

        result.push_back(p);
    }

    return result;
}

std::vector<ErrorPoint> compute_errors(
    const std::vector<Point>& solution
)
{
    std::vector<ErrorPoint> errors;

    for (const auto& p : solution)
    {
        double exact = exact_solution(p.x);

        errors.push_back(
            {
                p.x,
                exact,
                std::abs(p.y - exact)
            });
    }

    return errors;
}


static double compute_max_error(
    const std::vector<Point>& data
)
{
    double max_error = 0.0;

    for (const auto& p : data)
    {
        double exact = exact_solution(p.x);
        double error = std::abs(p.y - exact);
        if (error > max_error)
        {
            max_error = error;
        }
    }

    return max_error;
}


std::vector<AccuracyResult> accuracy_study(
    double a,
    double b,
    double h,
    double y0
)
{

    std::vector<AccuracyResult> results;
    std::vector<double> eps_values =
    {
        1e-2,
        1e-3,
        1e-4,
        1e-5,
        1e-6,
        1e-7,
        1e-8
    };

    for (double eps : eps_values)
    {
        auto solution = solve_euler_trapezoid(
            a,
            b,
            h,
            eps,
            y0
        );

        double max_error = compute_max_error(solution);

        AccuracyResult r;

        r.eps = eps;
        r.max_error = max_error;

        results.push_back(r);
    }

    return results;
}

void write_solution_csv(
    const std::string& filename,
    const std::vector<Point>& data
)
{
    std::ofstream file(filename);

    file << "x,y_numeric\n";

    for (const auto& p : data)
    {
        file
            << p.x << ","
            << p.y << "\n";
    }
}


void write_error_csv(
    const std::string& filename,
    const std::vector<ErrorPoint>& errors
)
{
    std::ofstream file(filename);

    file << "x,y_exact,error\n";

    for (const auto& p : errors)
    {
        file
            << p.x
            << ","
            << p.y_exact
            << ","
            << p.error
            << "\n";
    }
}


void write_accuracy_csv(
    const std::string& filename,
    const std::vector<AccuracyResult>& data
)
{
    std::ofstream file(filename);

    file << "eps,max_error\n";

    for (const auto& r : data)
    {
        file
            << r.eps << ","
            << r.max_error << "\n";
    }
}


int main()
{
    try
    {
        double a = std::numbers::pi / 2;
        double b = 2 * std::numbers::pi;
        double h = 0.01;
        double eps = pow(10, -8);
        double y0 = std::numbers::pi / 2;

        auto solution = solve_euler_trapezoid(
            a, b, h, eps, y0
        );

        auto errors = compute_errors(solution);

        auto accuracy = accuracy_study(
            a, b, h, y0
        );

        write_solution_csv(
            "../../data/solution.csv",
            solution
        );

        write_accuracy_csv(
            "../../data/accuracy.csv",
            accuracy
        );

        write_error_csv(
            "../../data/error.csv",
            errors
        );

    }
    catch (const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
    }

    return 0;
}