package com.jhlabs.map.proj;

import java.awt.geom.Point2D;

public class LatLongProjection extends Projection {

  private static final long serialVersionUID = 1L;

  @Override
  public Point2D.Double transformRadians(Point2D.Double src, Point2D.Double dst) {
    dst.x = src.x * RTD;
    dst.y = src.y * RTD;
    return dst;
  }

  @Override
  public Point2D.Double inverseTransformRadians(Point2D.Double src,
      Point2D.Double dst) {
    dst.x = src.x * DTR;
    dst.y = src.y * DTR;
    return dst;
  }
}
