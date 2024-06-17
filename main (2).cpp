#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

struct Spline {
    vector<double> x, y, a, b, c, d;
};

Spline cubicSpline(const vector<double>& x, const vector<double>& y) {
    int n = x.size() - 1;
    vector<double> h(n), alpha(n), l(n + 1), mu(n), z(n + 1);
    vector<double> a = y, b(n), c(n + 1), d(n);

    for (int i = 0; i < n; ++i)
        h[i] = x[i + 1] - x[i];

    for (int i = 1; i < n; ++i)
        alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1]);

    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n; ++i) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for (int j = n - 1; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    return {x, y, a, b, c, d};
}

double evaluateSpline(const Spline& spline, double x) {
    int n = spline.x.size() - 1;
    int i = 0;
    while (i < n && x > spline.x[i + 1]) ++i;

    double dx = x - spline.x[i];
    return spline.a[i] + spline.b[i] * dx + spline.c[i] * dx * dx + spline.d[i] * dx * dx * dx;
}

int main() {
    vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    vector<double> y = {0.0, 0.5, 0.86603, 1.0, 0.86603};

    Spline spline = cubicSpline(x, y);

    double x_star = 1.5;
    double y_star = evaluateSpline(spline, x_star);

    cout << fixed << setprecision(5) << "Значение функции в точке x* = " << x_star << ": " << y_star << endl;

    return 0;
}
