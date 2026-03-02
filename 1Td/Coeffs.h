#pragma once

#include <vector>

// Container for Aubin coefficients
struct AubinCoeffValues {
    double A0;
    double A1;
    double A2;
    double A3;
    double C;
    double eta;
};

// Species-specific coefficient functions
AubinCoeffValues AubinCoeffs_Argon(
    double eps,
    const std::vector<double>& a
);

AubinCoeffValues AubinCoeffs_Nitrogen(
    double eps,
    const std::vector<double>& a
);

AubinCoeffValues AubinCoeffs_N2(
    double eps,
    const std::vector<double>& a
);



