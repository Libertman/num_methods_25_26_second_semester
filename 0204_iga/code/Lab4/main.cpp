#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>


double P(const double& x)
{
	return pow(x, 5) - 6.2 * pow(x, 3) + 3.5 * pow(x, 2) - 7 * x;
}


double func(const double& x)
{
	return P(x) * cos(2 * x);
}


double transform(double a, double b, double t)
{
	return (b + a) / 2 + (b - a) / 2 * t;
}

double richardson(double Ih, double I2h, int p)
{
    return Ih + (Ih - I2h) / (pow(2, p) - 1);
}


double lobatto_4(double a, double b, double fa, double fb)
{
    double t[] = { -sqrt(1.0 / 5), sqrt(1.0 / 5) };
    double w[] = { 5.0 / 6, 5.0 / 6 };
    double w_end = 1.0 / 6;

    double sum = w_end * (fa + fb);

    for (int i = 0; i < 2; i++)
    {
        double x = transform(a, b, t[i]);
        sum += w[i] * func(x);
    }

    return (b - a) / 2 * sum;
}

double lobatto_5(double a, double b, double fa, double fb)
{
    double t[] = { -sqrt(3.0 / 7), 0.0, sqrt(3.0 / 7) };
    double w[] = { 49.0 / 90, 32.0 / 45, 49.0 / 90 };
    double w_end = 1.0 / 10;

    double sum = w_end * (fa + fb);

    for (int i = 0; i < 3; i++)
    {
        double x = transform(a, b, t[i]);
        sum += w[i] * func(x);
    }

    return (b - a) / 2 * sum;
}


struct AdaptiveResult
{
    double integral;
    int segments;
};

AdaptiveResult adaptive_integrate(double a, double b, double fa, double fb, double eps)
{
    double I4 = lobatto_4(a, b, fa, fb);
    double I5 = lobatto_5(a, b, fa, fb);

    double error = fabs(I5 - I4);

    if (error < eps) {
        return {I5, 1};
    }

    double mid = (a + b) / 2;
    double fmid = func(mid);

    AdaptiveResult left = adaptive_integrate(a, mid, fa, fmid, eps/2);
    AdaptiveResult right = adaptive_integrate(mid, b, fmid, fb, eps/2);

    return {left.integral + right.integral,
            left.segments + right.segments};
}


struct Params {
    double a, b;
};

Params readParams() {
    Params p;
    std::cout << "Enter a b: ";
    std::cin >> p.a >> p.b;
    return p;
}

int main()
{
    Params p = readParams();
    std::ofstream out("../../data/results.txt");

    double fa = func(p.a);
    double fb = func(p.b);

    for (double eps = 1e-1; eps >= 1e-8; eps /= 2)
    {
        AdaptiveResult I1 = adaptive_integrate(p.a, p.b, fa, fb, eps);
        AdaptiveResult I2 = adaptive_integrate(p.a, p.b, fa, fb, eps / 2);

        double I_rich = richardson(I1.integral, I2.integral, 6);

        double h = (p.b - p.a) / I1.segments;

        out << std::scientific << std::setprecision(3) << eps << " "
            << std::fixed << std::setprecision(15) << I1.integral << " "
            << std::fixed << std::setprecision(15) << I_rich << " "
            << std::fixed << std::setprecision(3) << h << "\n";
    }

    out.close();
}
