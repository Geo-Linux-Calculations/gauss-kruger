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

#include <math.h>

#ifndef GAUSSKRUGER_H
#define GAUSSKRUGER_H

#ifndef M_PI
#define M_PI 3.1415926535
#endif

typedef struct
{
    double coord[3];
}
GKcoord;

typedef struct
{
    double inverseFlattening;
    double flattening;
    double equatorialRadius;
    double centralMeridian;
    double scale;
    double falseNorthing;
    double falseEasting;
}
GKparams;

GKcoord gausskruger(GKcoord source, GKparams ellipsoid, int freverse);

#endif // GAUSSKRUGER_H
