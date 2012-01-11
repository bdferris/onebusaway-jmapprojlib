/*
Copyright 2006 Jerry Huxtable

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
 */

/*
 * This file was semi-automatically converted from the public-domain USGS PROJ source.
 *
 * Bernhard Jenny, 19 September 2010: fixed inverse spherical.
 * 23. September 2010: Change super class to ConicProjection.
 */
package com.jhlabs.map.proj;

import java.awt.geom.*;
import com.jhlabs.map.*;

public class BonneProjection extends ConicProjection {

    private double phi1;
    private double cphi1;
    private double am1;
    private double m1;
    private double[] en;

    public Point2D.Double project(double lplam, double lpphi, Point2D.Double out) {
        if (spherical) {
            double E, rh;

            rh = cphi1 + phi1 - lpphi;
            if (Math.abs(rh) > EPS10) {
                out.x = rh * Math.sin(E = lplam * Math.cos(lpphi) / rh);
                out.y = cphi1 - rh * Math.cos(E);
            } else {
                out.x = out.y = 0.;
            }
        } else {
            double rh, E, c;

            rh = am1 + m1 - MapMath.mlfn(lpphi, E = Math.sin(lpphi), c = Math.cos(lpphi), en);
            E = c * lplam / (rh * Math.sqrt(1. - es * E * E));
            out.x = rh * Math.sin(E);
            out.y = am1 - rh * Math.cos(E);
        }
        return out;
    }

    public Point2D.Double projectInverse(double xyx, double xyy, Point2D.Double out) {
        if (spherical) {
            double rh;

            rh = MapMath.distance(xyx, xyy = cphi1 - xyy);
            out.y = cphi1 + phi1 - rh;
            if (Math.abs(out.y) > MapMath.HALFPI) {
                throw new ProjectionException("I");
            }
            if (Math.abs(Math.abs(out.y) - MapMath.HALFPI) <= EPS10) {
                out.x = 0.;
            } else {
                out.x = rh * Math.atan2(xyx, xyy) / Math.cos(out.y);
            }
        } else {
            double s, rh;

            rh = MapMath.distance(xyx, out.y = am1 - xyy);
            out.y = MapMath.inv_mlfn(am1 + m1 - rh, es, en);
            if ((s = Math.abs(out.y)) < MapMath.HALFPI) {
                s = Math.sin(out.y);
                out.x = rh * Math.atan2(xyx, xyy)
                        * Math.sqrt(1. - es * s * s) / Math.cos(out.y);
            } else if (Math.abs(s - MapMath.HALFPI) <= EPS10) {
                out.x = 0.;
            } else {
                throw new ProjectionException("I");
            }
        }
        return out;
    }

    /**
     * Returns true if this projection is equal area
     */
    public boolean isEqualArea() {
        return true;
    }

    public boolean hasInverse() {
        return true;
    }

    public void initialize() {
        super.initialize();

        double c;

//		phi1 = pj_param(params, "rlat_1").f;
        phi1 = MapMath.HALFPI;
        if (Math.abs(phi1) < EPS10) {
            throw new ProjectionException("-23");
        }
        if (!spherical) {
            en = MapMath.enfn(es);
            m1 = MapMath.mlfn(phi1, am1 = Math.sin(phi1),
                    c = Math.cos(phi1), en);
            am1 = c / (Math.sqrt(1. - es * am1 * am1) * am1);
        } else {
            if (Math.abs(phi1) + EPS10 >= MapMath.HALFPI) {
                cphi1 = 0.;
            } else {
                cphi1 = 1. / Math.tan(phi1);
            }
        }
    }

    public String toString() {
        return "Bonne";
    }
}
