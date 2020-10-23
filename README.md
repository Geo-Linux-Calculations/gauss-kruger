# Simple Gauss-Krüger transformation for C

## Description

This project provides a C implementation of coordinate transformation between geodetic
coordinates and grid coordinates of Gauss (Transverse Mercator) projections, using the Krüger-n series developed to fourth order.

This is a simplification of the C++ implementation: https://github.com/f03el/gauss-kruger-cpp .

Except the library, a simple command line tool is included so you can play around with the transformations.

## Dependencies

The library itself requires no more than standard library components.

## Build

To build, type:
```sh
$ make
```

## Example

Example of direct transformation from SWEREF 99 to RT90 (same values as from the library usage example)::
```sh
$ ./gausskruger -i 298.257222101 -a 6378137 -m 11.30625 -s 1.000006 -n -667.282 -e 1500025.141 67.877566667 21.06025
Ellipsoid:
  Inverse flattening = 298.257222
  Equatorial radius = 6378137.000000
  Central meridian = 11.306250
  Scale factor = 1.000006
  False northing = -667.282000
  False easting = 1500025.141000
Point:
  Latitude = 67.8775666670
  Longitude = 21.0602500000
  Northing = 7563929.5303
  Easting = 1908686.7147
```

To see all options, run the tool without options or with `-h`:
```sh
$ ./gausskruger -h
Usage:
./gausskruger [options] latitude longitude
Options:
    -i      Inverse flattening of the ellipsoid (default = 298.257222)
    -a      Equatorial radius (default = 6378137.000000)
    -m      Longitude of the central meridian (default = 11.306250)
    -s      Scale factor along the central meridian (default = 1.000006)
    -n      False northing (default = -667.282000)
    -e      False easting (default = 1500025.141000)
    -r      Reverse transformation (default = FALSE)
    -h      this help
```

## License

Public Domain Mark 1.0  
 No Copyright

---

https://github.com/Geo-Linux-Calculations/gauss-kruger
