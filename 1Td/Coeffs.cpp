#include "Coeffs.h"
#include <cmath>

AubinCoeffValues AubinCoeffs_Argon(
    double eps,
    const std::vector<double>& a
) {
    AubinCoeffValues v;

    if (eps <= 0.0) eps = 1e-6;

    double x = std::log(eps);

    v.C = a[0]*std::pow(x,4.0)-
        + a[1]*std::pow(x,3.0)
        + a[2]*std::pow(x,2.0)
        + a[3]*x
        + a[4];

    v.eta = a[5]*std::exp(-eps/a[6]) + a[7];

    if (eps <= 15.0) {
        v.A1 = a[8]
             + std::pow(eps, a[11]) *
               (a[9]*std::exp(-std::pow(eps/a[12], a[13]))
              -a[10]*std::exp(-std::pow(eps/a[14], a[15])))
             + a[16]*eps
             + a[17]*std::pow(eps, 2.0);
    } else {
        double xh = std::log(eps);
        v.A1 = (a[18]*xh + a[19])*xh
             + a[20]
             + a[21]*std::pow(xh,3.0)
             + a[22]*std::pow(xh,4.0);
    }

    v.A2 = (a[23]*std::pow(a[24]/2.0,2.0)
           /(std::pow(eps - a[25],2.0)
           + std::pow(a[24]/2.0,2.0))
           + a[26]) + a[27]*eps;

    v.A3 = a[28]*std::pow(eps, a[29])
         * std::exp(-std::pow(eps/a[30], a[31]))
         + (a[32]*eps)/(1.0 + a[33]*eps);

    v.A0 = 1.0;

    return v;
}

AubinCoeffValues AubinCoeffs_Nitrogen(
    double eps,
    const std::vector<double>& a
) {
    AubinCoeffValues v;

    if (eps <= 0.0) eps = 1e-6;

    double x = std::log(eps);

    v.C = a[0]
        + a[1]*x
        + a[2]*std::pow(x, 2.0)
        + a[3]*std::pow(x, 3.0);

    v.eta = a[4]
        + a[5]*x
        + a[6]*std::pow(x, 2.0)
        + a[7]*std::pow(x, 3.0);

    if (eps <= 20.0) {
        v.A1 = a[8]
            + a[9]*x
            + a[10]*std::pow(x, 2.0)
            + a[11]*std::pow(x, 3.0)
            + a[12]*x*std::exp(-x/0.25);

    } else {
        v.A1 = a[13]
            + a[14]*x
            + a[15]*std::pow(x, 2.0)
            + a[16]*std::pow(x, 3.0);
    }

    v.A2 = a[17]
     + a[18]*x
     + a[19]*std::pow(x, 2.0)
     + a[20]*std::pow(x, 3.0)
     + a[21]*std::pow(x, 4.0);

    v.A3 = a[22]
     + a[23]*eps
     + a[24]*std::pow(eps, 2.0)
     + a[25]*std::pow(eps, 3.0);


    v.A0 = 1.0;

    return v;
}

AubinCoeffValues AubinCoeffs_N2(
    double eps,
    const std::vector<double>& a
) {
    AubinCoeffValues v;

    if (eps <= 0.0) eps = 1e-6;

    double x = std::log(eps);

    v.C = a[0]
        + a[1]*x;

    v.eta = a[2]
        + a[3]*x;

    if (eps <= 20.0) {
        v.A1 = a[4]
            + a[5]*x
            + a[6]*std::pow(x, 2.0);

    } else {
        v.A1 = a[7]
            + a[8]*x;
    }

    v.A2 =
        ( a[9]
        + a[10]*eps
        + a[11]*std::pow(eps, 2.0) )
    / ( 1.0
        + a[12]*eps
        + a[13]*std::pow(eps, 2.0)
        + a[14]*std::pow(eps, 3.0) );


    v.A3 =
        ( a[15]
        + a[16]*eps
        + a[17]*std::pow(eps, 2.0) )
    / ( 1.0
        + a[18]*eps
        + a[19]*std::pow(eps, 2.0)
        + a[20]*std::pow(eps, 3.0) );

    v.A0 = 1.0;

    return v;
}

