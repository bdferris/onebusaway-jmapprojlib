package com.jhlabs.map.proj;

import java.awt.geom.Point2D;

import com.jhlabs.map.Ellipsoid;

public class CoordinateSystemToCoordinateSystem {

  private final static double PI_OVER_2 = Math.PI / 2;

  /**
   * cosine of 67.5 degrees
   */
  private final static double COS_67P5 = 0.38268343236508977;

  /**
   * Toms region 1 constant
   */
  private final static double AD_C = 1.0026000;

  private static final double genau = 1.E-12;
  private static final double genau2 = (genau * genau);
  private static final int maxiter = 30;

  public static Point2D.Double transform(Projection fromProjection,
      Projection toProjection, Point2D.Double from, Point2D.Double to) {

    Point2D.Double intermediate = new Point2D.Double();
    fromProjection.inverseTransformRadians(from, intermediate);

    transformDatum(fromProjection, toProjection, intermediate, intermediate);

    return toProjection.transformRadians(intermediate, to);
  }

  public static void transformDatum(Projection fromProjection,
      Projection toProjection, Point2D.Double from, Point2D.Double to) {

    double src_a = fromProjection.a;
    double src_es = fromProjection.es;

    double dst_a = toProjection.a;
    double dst_es = toProjection.es;

    if (src_a == dst_a && src_es == dst_es) {
      to.x = from.x;
      to.y = from.y;
      return;
    }

    Point3D p3 = new Point3D();
    ConvertGeodeticToGeocentric(src_a, src_es, from, p3);
    GeocentricToWgs84(fromProjection.ellipsoid, p3);
    ConvertGeocentricToGeodeticIterative(dst_a, dst_es, p3, to);
  }

  /**
   * Converts geodetic coordinates (latitude, longitude) to geocentric
   * coordinates (X, Y), according to the current ellipsoid parameters.
   * 
   * @param a Semi-major axis of ellipsoid in meters
   * @param es Eccentricity squared
   * @param from - geodetic latitude and longitude in radians
   * @param to - resulting geocentric x and y
   */
  private static void ConvertGeodeticToGeocentric(double a, double es,
      Point2D.Double from, Point3D to) {

    double lat = from.y;
    double lon = from.x;
    double height = 0;

    /**
     * Don't blow up if Latitude is just a little out of the value* range as it
     * may just be a rounding issue. Also removed longitude* test, it should be
     * wrapped by cos() and sin(). NFW for PROJ.4, Sep/2001.
     */
    if (lat < -PI_OVER_2 && lat > -1.001 * PI_OVER_2) {
      lat = -PI_OVER_2;
    } else if (lat > PI_OVER_2 && lat < 1.001 * PI_OVER_2) {
      lat = PI_OVER_2;
    } else if ((lat < -PI_OVER_2) || (lat > PI_OVER_2)) {
      /*
       * Latitude out of range
       */
      throw new ProjectionException("latitude is out of range " + lat);
    }

    if (lon > Math.PI) {
      lon -= (2 * Math.PI);
    }

    double sinLat = Math.sin(lat);
    double cosLat = Math.cos(lat);
    double sin2Lat = sinLat * sinLat;
    // Earth radius at location
    double radiusOfEarthAtLocation = a / (Math.sqrt(1.0e0 - es * sin2Lat));
    to.x = (radiusOfEarthAtLocation + height) * cosLat * Math.cos(lon);
    to.y = (radiusOfEarthAtLocation + height) * cosLat * Math.sin(lon);
    to.z = ((radiusOfEarthAtLocation * (1 - es)) + height) * sinLat;
  }

  private static void GeocentricToWgs84(Ellipsoid ellipsoid, Point3D point) {
    double[] params = ellipsoid.datumParams;
    if (params == null) {
      return;
    }
    if (params.length == 3) {
      point.x += params[0];
      point.y += params[1];
      point.z += params[2];
    } else if (params.length == 7) {
      double x = point.x;
      double y = point.y;
      double z = point.z;
      point.x = params[6] * (x - params[5] * y + params[4] * z) + params[0];
      point.y = params[6] * (params[5] * x + y - params[3] * z) + params[1];
      point.z = params[6] * (-params[4] * x + params[3] * y + z) + params[2];
    }
  }

  /**
   * The method used here is derived from 'An Improved Algorithm for Geocentric
   * to Geodetic Coordinate Conversion', by Ralph Toms, Feb 1996
   */
  private static void ConvertGeocentricToGeodeticNonIterative(double a,
      double es2, Point3D from, Point2D.Double to) {

    double b = (es2 == 0.0) ? a : a * Math.sqrt(1 - Math.sqrt(es2));

    // double S1;
    // double Sin_B0;
    double Sin3_B0; /* cube of sin(B0) */
    double Cos_B0; /* cos(B0) */
    double Sin_p1; /* sin(phi1), phi1 is estimated latitude */
    double Cos_p1; /* cos(phi1) */
    double Rn; /* Earth radius at location */
    double Sum; /* numerator of cos(phi1) */

    double X = from.x;
    double Y = from.y;
    double Z = from.z;

    /* indicates location is in polar region */
    boolean At_Pole = false;

    if (X != 0.0) {
      to.y = Math.atan2(Y, X);
    } else {
      if (Y > 0) {
        to.x = PI_OVER_2;
      } else if (Y < 0) {
        to.x = -PI_OVER_2;
      } else {
        At_Pole = true;
        to.x = 0.0;
        if (Z > 0.0) { /* north pole */
          to.y = PI_OVER_2;
        } else if (Z < 0.0) { /* south pole */
          to.y = -PI_OVER_2;
        } else { /* center of earth */
          to.y = PI_OVER_2;
          // height = -Geocent_b;
          return;
        }
      }
    }

    double squareOfDistanceFromZAxis = X * X + Y * Y;
    double distanceFromZAxis = Math.sqrt(squareOfDistanceFromZAxis);
    double initialEstimateOfVerticalComponent = Z * AD_C;
    double initialEstimateOfHorizontalComponent = Math.sqrt(initialEstimateOfVerticalComponent
        * initialEstimateOfVerticalComponent + squareOfDistanceFromZAxis);
    /* sin(B0), B0 is estimate of Bowring aux variable */
    double sin_B0 = initialEstimateOfVerticalComponent
        / initialEstimateOfHorizontalComponent;
    Cos_B0 = distanceFromZAxis / initialEstimateOfHorizontalComponent;
    Sin3_B0 = sin_B0 * sin_B0 * sin_B0;
    double correctedEstimateOfVerticalComponent = Z + b * es2 * Sin3_B0;
    Sum = distanceFromZAxis - a * es2 * Cos_B0 * Cos_B0 * Cos_B0;
    /* corrected estimate of horizontal component */
    double correctedEstimateOfHorizontalComponent = Math.sqrt(correctedEstimateOfVerticalComponent
        * correctedEstimateOfVerticalComponent + Sum * Sum);
    Sin_p1 = correctedEstimateOfVerticalComponent
        / correctedEstimateOfHorizontalComponent;
    Cos_p1 = Sum / correctedEstimateOfHorizontalComponent;
    Rn = a / Math.sqrt(1.0 - es2 * Sin_p1 * Sin_p1);
    if (Cos_p1 >= COS_67P5) {
      // height = W / Cos_p1 - Rn;
    } else if (Cos_p1 <= -COS_67P5) {
      // height = W / -Cos_p1 - Rn;
    } else {
      // height = Z / Sin_p1 + Rn * (es2 - 1.0);
    }
    if (At_Pole == false) {
      to.y = Math.atan(Sin_p1 / Cos_p1);
    }
  }

  private static void ConvertGeocentricToGeodeticIterative(double a, double es,
      Point3D from, Point2D.Double to) {

    double P; /* distance between semi-minor axis and location */
    double RR; /* distance between center and location */
    double CT; /* sin of geocentric latitude */
    double ST; /* cos of geocentric latitude */
    double RX;
    double RK;
    double RN; /* Earth radius at location */
    double CPHI0; /* cos of start or old geodetic latitude in iterations */
    double SPHI0; /* sin of start or old geodetic latitude in iterations */
    double CPHI; /* cos of searched geodetic latitude */
    double SPHI; /* sin of searched geodetic latitude */
    double SDPHI; /*
                   * end-criterium: addition-theorem of
                   * sin(Latitude(iter)-Latitude(iter-1))
                   */
    boolean At_Pole; /* indicates location is in polar region */
    int iter; /* # of continous iteration, max. 30 is always enough (s.a.) */

    double X = from.x;
    double Y = from.y;
    double Z = from.z;

    double b = (es == 0.0) ? a : a * Math.sqrt(1 - Math.sqrt(es));
    double height = 0;

    At_Pole = false;
    P = Math.sqrt(X * X + Y * Y);
    RR = Math.sqrt(X * X + Y * Y + Z * Z);

    /* special cases for latitude and longitude */
    if (P / a < genau) {

      /* special case, if P=0. (X=0., Y=0.) */
      At_Pole = true;
      to.x = 0.;

      /*
       * if (X,Y,Z)=(0.,0.,0.) then Height becomes semi-minor axis of ellipsoid
       * (=center of mass), Latitude becomes PI/2
       */
      if (RR / a < genau) {
        to.y = PI_OVER_2;
        height = b;
        return;

      }
    } else {
      /*
       * ellipsoidal (geodetic) longitude interval: -PI < Longitude <= +PI
       */
      to.x = Math.atan2(Y, X);
    }

    /*
     * -------------------------------------------------------------- Following
     * iterative algorithm was developped by "Institut für Erdmessung",
     * University of Hannover, July 1988. Internet: www.ife.uni-hannover.de
     * Iterative computation of CPHI,SPHI and Height. Iteration of CPHI and SPHI
     * to 10**-12 radian resp. 2*10**-7 arcsec.
     * --------------------------------------------------------------
     */
    CT = Z / RR;
    ST = P / RR;
    RX = 1.0 / Math.sqrt(1.0 - es * (2.0 - es) * ST * ST);
    CPHI0 = ST * (1.0 - es) * RX;
    SPHI0 = CT * RX;
    iter = 0;

    /*
     * loop to find sin(Latitude) resp. Latitude until
     * |sin(Latitude(iter)-Latitude(iter-1))| < genau
     */
    do {
      iter++;
      RN = a / Math.sqrt(1.0 - es * SPHI0 * SPHI0);

      /* ellipsoidal (geodetic) height */
      height = P * CPHI0 + Z * SPHI0 - RN * (1.0 - es * SPHI0 * SPHI0);

      RK = es * RN / (RN + height);
      RX = 1.0 / Math.sqrt(1.0 - RK * (2.0 - RK) * ST * ST);
      CPHI = ST * (1.0 - RK) * RX;
      SPHI = CT * RX;
      SDPHI = SPHI * CPHI0 - CPHI * SPHI0;
      CPHI0 = CPHI;
      SPHI0 = SPHI;
    } while (SDPHI * SDPHI > genau2 && iter < maxiter);

    /* ellipsoidal (geodetic) latitude */
    to.y = Math.atan(SPHI / Math.abs(CPHI));

    return;

  }

  private static class Point3D {
    public double x;
    public double y;
    public double z;
  }
}
