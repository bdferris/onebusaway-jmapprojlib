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
 * Bernhard Jenny, 21 September 2010: Changed name from MBTFPSProjection to
 * McBrydeThomasFlatPolarSine2Projection.
 * 23 September 2010: changed super class to PseudoCylindricalProjection, added
 * isEqualArea.
 */
package com.jhlabs.map.proj;

import java.awt.geom.*;
import com.jhlabs.map.*;

public class McBrydeThomasFlatPolarSine2Projection extends PseudoCylindricalProjection {

    private final static int MAX_ITER = 10;
    private final static double LOOP_TOL = 1e-7;
    private final static double C1 = 0.45503;
    private final static double C2 = 1.36509;
    private final static double C3 = 1.41546;
    private final static double C_x = 0.22248;
    private final static double C_y = 1.44492;
    private final static double C1_2 = 0.33333333333333333333333333;

    public Point2D.Double project(double lplam, double lpphi, Point2D.Double out) {
        double k, V, t;
        int i;

        k = C3 * Math.sin(lpphi);
        for (i = MAX_ITER; i > 0; i--) {
            t = lpphi / C2;
            out.y -= V = (C1 * Math.sin(t) + Math.sin(lpphi) - k)
                    / (C1_2 * Math.cos(t) + Math.cos(lpphi));
            if (Math.abs(V) < LOOP_TOL) {
                break;
            }
        }
        t = lpphi / C2;
        out.x = C_x * lplam * (1. + 3. * Math.cos(lpphi) / Math.cos(t));
        out.y = C_y * Math.sin(t);
        return out;
    }

    public Point2D.Double projectInverse(double xyx, double xyy, Point2D.Double out) {
        double t, s;

        out.y = C2 * (t = MapMath.asin(xyy / C_y));
        out.x = xyx / (C_x * (1. + 3. * Math.cos(out.y) / Math.cos(t)));
        out.y = MapMath.asin((C1 * Math.sin(t) + Math.sin(out.y)) / C3);
        return out;
    }

    public boolean hasInverse() {
        return true;
    }

    public boolean isEqualArea() {
        return true; // FIXME verify if correct.
    }

    public String toString() {
        return "McBryde-Thomas Flat-Pole Sine (No. 2)";
    }
}
