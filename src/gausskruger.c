/* Gauss Kruger.
 * Version: 1.2
 *
 * Authors: Erik Lundin [https://github.com/f03el] 2016.
 * Authors: zvezdochiot [https://github.com/zvezdochiot] 2020.
 *
 * Unlicense.
 * This is free and unencumbered software released into the public domain.
 * For more information, please refer to <https://unlicense.org>
*/

#include "gausskruger.h"

GKcoord gausskruger(GKcoord source, GKparams ellipsoid, int freverse)
{
    double e2, n, rectifyingRadius;
    double A, B, C, D;
    double phi, phis, phic, lambda, lambda0, xi, eta;
    double deltaLambda, phiStar, xiPrim, etaPrim;
    double beta1, beta2, beta3, beta4;
    GKcoord target;

    e2 = ellipsoid.flattening * (2.0 - ellipsoid.flattening); // e2: first eccentricity squared
    n = ellipsoid.flattening / (2.0 - ellipsoid.flattening); // n: 3rd flattening
    rectifyingRadius = ellipsoid.equatorialRadius / (1.0 + n) * (1.0 + 0.25 * n * n + 0.015625 * n * n * n * n);

    if (freverse == 0)
    {
        A = e2;
        B = (5.0 * e2 * e2 - e2 * e2 * e2) / 6.0;
        C = (104.0 * e2 * e2 * e2 - 45.0 * e2 * e2 * e2 * e2) / 120.0;
        D = (1237.0 * e2 * e2 * e2 * e2) / 1260.0;

        // Latitude and longitude are expected to be given in degrees
        // phi and lambda: latitude and longitude in radians
        phi = source.coord[0] * M_PI / 180.0;
        lambda = source.coord[1] * M_PI / 180.0;
        lambda0 = ellipsoid.centralMeridian * M_PI / 180.0;

        // deltaLambda: longitude relative to the central meridian
        deltaLambda = lambda - lambda0;

        // phiStar: conformal latitude
        phis = sin(phi);
        phic = cos(phi);
        phiStar = phi - phis * phic *
                (A + B * phis * phis + C * phis * phis * phis * phis + D * phis * phis * phis * phis * phis * phis);

        xiPrim = atan(tan(phiStar) / cos(deltaLambda));
        etaPrim = atanh(cos(phiStar) * sin(deltaLambda));

        beta1 = 1.0/2.0 * n - 2.0/3.0 * n * n + 5.0/16.0 * n * n * n + 41.0/180.0 * n * n * n * n;
        beta2 = 13.0/48.0 * n * n - 3.0/5.0 * n * n * n + 557.0/1440.0 * n * n * n * n;
        beta3 = 61.0/240.0 * n * n * n - 103.0/140.0 * n * n * n * n;
        beta4 = 49561.0/161280.0 * n * n * n * n;

        target.coord[0] = ellipsoid.falseNorthing
                          + ellipsoid.scale * rectifyingRadius * (xiPrim
                                                + beta1 * sin(2.0*xiPrim) * cosh(2.0*etaPrim)
                                                + beta2 * sin(4.0*xiPrim) * cosh(4.0*etaPrim)
                                                + beta3 * sin(6.0*xiPrim) * cosh(6.0*etaPrim)
                                                + beta4 * sin(8.0*xiPrim) * cosh(8.0*etaPrim));
        target.coord[1] = ellipsoid.falseEasting
                          + ellipsoid.scale * rectifyingRadius * (etaPrim
                                                + beta1 * cos(2.0*xiPrim) * sinh(2.0*etaPrim)
                                                + beta2 * cos(4.0*xiPrim) * sinh(4.0*etaPrim)
                                                + beta3 * cos(6.0*xiPrim) * sinh(6.0*etaPrim)
                                                + beta4 * cos(8.0*xiPrim) * sinh(8.0*etaPrim));
    } else {
        xi = (source.coord[0] - ellipsoid.falseNorthing) / (ellipsoid.scale * rectifyingRadius);
        eta = (source.coord[1] - ellipsoid.falseEasting) / (ellipsoid.scale * rectifyingRadius);

        beta1 = 1.0/2.0 * n - 2.0/3.0 * n * n + 37.0/96.0 * n * n * n - 1.0/360.0 * n * n * n * n;
        beta2 = 1.0/48.0 * n * n + 1.0/15.0 * n * n * n - 437.0/1440.0 * n * n * n * n;
        beta3 = 17.0/480.0 * n * n * n - 37.0/840.0 * n * n * n * n;
        beta4 = 4397.0/161280.0 * n * n * n * n;

        xiPrim = xi
                - beta1 * sin(2.0*xi) * cosh(2.0*eta)
                - beta2 * sin(4.0*xi) * cosh(4.0*eta)
                - beta3 * sin(6.0*xi) * cosh(6.0*eta)
                - beta4 * sin(8.0*xi) * cosh(8.0*eta);
        etaPrim = eta
                - beta1 * cos(2.0*xi) * sinh(2.0*eta)
                - beta2 * cos(4.0*xi) * sinh(4.0*eta)
                - beta3 * cos(6.0*xi) * sinh(6.0*eta)
                - beta4 * cos(8.0*xi) * sinh(8.0*eta);

        phiStar = asin(sin(xiPrim) / cosh(etaPrim)); // Conformal latitude
        deltaLambda = atan(sinh(etaPrim) / cos(xiPrim));

        A = e2 + e2 * e2 + e2 * e2 * e2 + e2 * e2 * e2 * e2;
        B = -(7.0 * e2 * e2 + 17.0 * e2 * e2 * e2 + 30.0 * e2 * e2 * e2 * e2) / 6.0;
        C = (224.0 * e2 * e2 * e2 + 889.0 * e2 * e2 * e2 * e2) / 120.0;
        D = -(4279.0 * e2 * e2 * e2 * e2) / 1260.0;

        phis = sin(phiStar);
        phic = cos(phiStar);

        phi = phiStar
                + phis * phic * ( A
                                + B * phis * phis
                                + C * phis * phis * phis * phis
                                + D * phis * phis * phis * phis * phis * phis);

        // phi: latitude in radians, lambda: longitude in radians
        // Return latitude and longitude as degrees
        target.coord[0] = phi * 180.0 / M_PI;
        target.coord[1] = ellipsoid.centralMeridian + deltaLambda * 180.0 / M_PI;
    }

    return target;
}
