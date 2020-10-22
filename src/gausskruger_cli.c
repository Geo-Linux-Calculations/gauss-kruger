/* Gauss Kruger.
 * Version: 1.0
 *
 * Authors: Erik Lundin [https://github.com/f03el] 2016.
 * Authors: zvezdochiot [https://github.com/zvezdochiot] 2020.
 *
 * Public Domain Mark 1.0
 * No Copyright
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "gausskruger.h"

#define INVERSEFLATTERING 298.257222101
#define EQUATORIALRADIUS 6378137.0
#define CENTRALMERIDIAN 11.30625
#define SCALE 1.000006
#define FALSENORTHING -667.282
#define FALSEEASTING 1500025.141

void gausskrugerhelp(char *pname)
{
    printf("Gauss Kruger.\n");
    printf("Homepage: https://github.com/Geo-Linux-Calculations/gauss-kruger\n\n");
    printf("Usage:\n");
    printf("%s [options] latitude longitude\n", pname);
    printf("Options:\n");
    printf("    -i      Inverse flattening of the ellipsoid (default = %f)\n", INVERSEFLATTERING);
    printf("    -a      Equatorial radius (default = %f)\n", EQUATORIALRADIUS);
    printf("    -m      Longitude of the central meridian (default = %f)\n", CENTRALMERIDIAN);
    printf("    -s      Scale factor along the central meridian (default = %f)\n", SCALE);
    printf("    -n      False northing (default = %f)\n", FALSENORTHING);
    printf("    -e      False easting (default = %f)\n", FALSEEASTING);
    printf("    -r      Reverse transformation (default = FALSE)\n");
    printf("    -h      this help\n");
}

int main(int argc, char *argv[])
{
    GKparams ellipsoid;
    ellipsoid.inverseFlattening = INVERSEFLATTERING;
    ellipsoid.equatorialRadius = EQUATORIALRADIUS;
    ellipsoid.centralMeridian = CENTRALMERIDIAN;
    ellipsoid.scale = SCALE;
    ellipsoid.falseNorthing = FALSENORTHING;
    ellipsoid.falseEasting = FALSEEASTING;
    GKcoord coords;
    int opt, fReverse = 0, fhelp = 0;

    while ((opt = getopt(argc, argv, ":i:a:m:s:n:e:rh")) != -1)
    {
        switch(opt)
        {
            case 'i':
                ellipsoid.inverseFlattening = atof(optarg);
                break;
            case 'a':
                ellipsoid.equatorialRadius = atof(optarg);
                break;
            case 'm':
                ellipsoid.centralMeridian = atof(optarg);
                break;
            case 's':
                ellipsoid.scale = atof(optarg);
                break;
            case 'n':
                ellipsoid.falseNorthing = atof(optarg);
                break;
            case 'e':
                ellipsoid.falseEasting = atof(optarg);
                break;
            case 'r':
                fReverse = 1;
                break;
            case 'h':
                fhelp = 1;
                break;
            case ':':
                printf("option needs a value\n");
                break;
            case '?':
                printf("unknown option: %c\n", optopt);
                break;
        }
    }

    if(optind + 2 > argc || fhelp)
    {
        gausskrugerhelp(argv[0]);
        return 0;
    }
    ellipsoid.flattening = (ellipsoid.inverseFlattening > 0.0) ? (1.0 / ellipsoid.inverseFlattening) : 0.0;
    coords.coord[0] = atof(argv[optind]);
    coords.coord[1] = atof(argv[optind + 1]);

    printf("Ellipsoid:\n");
    printf("  Inverse flattening = %f\n", ellipsoid.inverseFlattening);
    printf("  Equatorial radius = %f\n", ellipsoid.equatorialRadius);
    printf("  Central meridian = %f\n", ellipsoid.centralMeridian);
    printf("  Scale factor = %f\n", ellipsoid.scale);
    printf("  False northing = %f\n", ellipsoid.falseNorthing);
    printf("  False easting = %f\n", ellipsoid.falseEasting);

    printf("Point:\n");
    if (fReverse == 1) {
        printf("  Northing = %.3f\n", coords.coord[0]);
        printf("  Easting = %.3f\n", coords.coord[1]);
        coords = gausskruger_gridtogeodetic(coords, ellipsoid);
        printf("  Latitude = %.10f\n", coords.coord[0]);
        printf("  Longitude = %.10f\n", coords.coord[1]);
    } else {
        printf("  Latitude = %.10f\n", coords.coord[0]);
        printf("  Longitude = %.10f\n", coords.coord[1]);
        coords = gausskruger_geodetictogrid(coords, ellipsoid);
        printf("  Northing = %.3f\n", coords.coord[0]);
        printf("  Easting = %.3f\n", coords.coord[1]);
    }

    return 0;
}
