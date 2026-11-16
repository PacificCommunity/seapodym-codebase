#ifndef SEAPODYM_UPWIND_SPLIT_H
#define SEAPODYM_UPWIND_SPLIT_H

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace seapodym_upwind {

inline double parse_epsilon(const char* text, const char* setting)
{
    char* end = nullptr;
    errno = 0;
    const double value = std::strtod(text, &end);
    if (end == text || *end || errno == ERANGE || !std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(std::string(setting) + " must be finite and nonnegative");
    }
    return value;
}

inline bool in_band(double velocity, double width)
{
    return std::abs(velocity) < width;
}

struct Split {
    double positive;
    double negative;
    double d_positive;
    double d_negative;
};

// Precondition: |velocity| < width. The caller retains the native stencil elsewhere.
inline Split inside_split(double velocity, double width)
{
    const double relative = velocity / width;
    const double d_positive = (1.0 + relative) * 0.5;
    const double d_negative = (1.0 - relative) * 0.5;
    return {
        (width * d_positive) * d_positive,
        -(width * d_negative) * d_negative,
        d_positive,
        d_negative,
    };
}

} // namespace seapodym_upwind

#endif
