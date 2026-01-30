#include "roots.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

namespace {
constexpr double TOL = 1e-6;
constexpr int MAX_ITERS = 1'000'000;

inline bool finite(double x) { return std::isfinite(x); }
inline bool near_zero(double x) { return std::abs(x) <= TOL; }

inline bool same_sign(double a, double b) {
    return (a > 0.0 && b > 0.0) || (a < 0.0 && b < 0.0);
}
}
 
bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);

    double fa = f(a);
    double fb = f(b);
    if (!finite(fa) || !finite(fb)) return false;

    if (near_zero(fa)) { *root = a; return true; }
    if (near_zero(fb)) { *root = b; return true; }

    // Must bracket a sign change
    if (same_sign(fa, fb)) return false;

    double left = a, right = b;
    double fleft = fa, fright = fb;

    for (int i = 0; i < MAX_ITERS; ++i) {
        double mid = 0.5 * (left + right);
        double fmid = f(mid);
        if (!finite(fmid)) return false;

        if (near_zero(fmid) || std::abs(right - left) <= TOL) {
            *root = mid;
            return true;
        }

        if (same_sign(fleft, fmid)) {
            left = mid;
            fleft = fmid;
        } else {
            right = mid;
            fright = fmid;
        }
    }

    *root = 0.5 * (left + right);
    return true;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);

    double fa = f(a);
    double fb = f(b);
    if (!finite(fa) || !finite(fb)) return false;

    if (near_zero(fa)) { *root = a; return true; }
    if (near_zero(fb)) { *root = b; return true; }

    // Must bracket a sign change
    if (same_sign(fa, fb)) return false;

    double left = a, right = b;
    double fleft = fa, fright = fb;

    for (int i = 0; i < MAX_ITERS; ++i) {
        double denom = (fright - fleft);
        if (std::abs(denom) < 1e-15) {
            // Degenerate slope, fall back to midpoint-ish answer
            double mid = 0.5 * (left + right);
            double fmid = f(mid);
            if (!finite(fmid)) return false;
            *root = mid;
            return near_zero(fmid) || std::abs(right - left) <= TOL;
        }

        // False position (linear interpolation)
        double x = (left * fright - right * fleft) / denom;
        double fx = f(x);
        if (!finite(fx)) return false;

        if (near_zero(fx) || std::abs(right - left) <= TOL) {
            *root = x;
            return true;
        }

        if (same_sign(fleft, fx)) {
            left = x;
            fleft = fx;
        } else {
            right = x;
            fright = fx;
        }
    }

    // Best guess after max iters
    double denom = (fright - fleft);
    if (std::abs(denom) < 1e-15) {
        *root = 0.5 * (left + right);
    } else {
        *root = (left * fright - right * fleft) / denom;
    }
    return true;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);

    // Start guess must be in [a,b]
    if (c < a || c > b) return false;

    double x = c;

    for (int i = 0; i < MAX_ITERS; ++i) {
        double fx = f(x);
        double gx = g(x);
        if (!finite(fx) || !finite(gx)) return false;

        if (near_zero(fx)) { *root = x; return true; }

        // derivative too small => unstable
        if (std::abs(gx) < 1e-15) return false;

        double x_next = x - fx / gx;
        if (!finite(x_next)) return false;

        // fail if iteration leaves interval
        if (x_next < a || x_next > b) return false;

        if (std::abs(x_next - x) <= TOL) {
            *root = x_next;
            return true;
        }

        x = x_next;
    }

    *root = x;
    return true; // best effort
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);

    // We’ll use a and b as the two initial points (classic secant),
    // and require c to be in interval as per your signature (starting guess).
    if (c < a || c > b) return false;

    double x0 = a;
    double x1 = b;

    double f0 = f(x0);
    double f1 = f(x1);
    if (!finite(f0) || !finite(f1)) return false;

    if (near_zero(f0)) { *root = x0; return true; }
    if (near_zero(f1)) { *root = x1; return true; }

    for (int i = 0; i < MAX_ITERS; ++i) {
        double denom = (f1 - f0);
        if (std::abs(denom) < 1e-15) return false;

        double x2 = x1 - f1 * (x1 - x0) / denom;
        if (!finite(x2)) return false;

        // Lab note: fail if iteration leaves interval
        if (x2 < a || x2 > b) return false;

        if (std::abs(x2 - x1) <= TOL) {
            *root = x2;
            return true;
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;

        f1 = f(x1);
        if (!finite(f1)) return false;
        if (near_zero(f1)) { *root = x1; return true; }
    }

    *root = x1;
    return true;
}


