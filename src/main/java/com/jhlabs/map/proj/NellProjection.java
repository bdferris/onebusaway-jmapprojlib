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
 * Bernhard Jenny, 19 September 2010: fixed forward projection.
 * 23 September 2010: changed super class to PseudoCylindricalProjection.
 */
package com.jhlabs.map.proj;

import java.awt.geom.*;
import com.jhlabs.map.*;

public class NellProjection extends PseudoCylindricalProjection {

    private final static int MAX_ITER = 10;
    private final static double LOOP_TOL = 1e-7;

    public Point2D.Double project(double lplam, double lpphi, Point2D.Double out) {
        double k, V;
        int i;

        k = 2. * Math.sin(lpphi);
        V = lpphi * lpphi;
        lpphi *= 1.00371 + V * (-0.0935382 + V * -0.011412);
        for (i = MAX_ITER; i > 0; --i) {
            lpphi -= V = (lpphi + Math.sin(lpphi) - k) / (1. + Math.cos(lpphi));
            if (Math.abs(V) < LOOP_TOL) {
                break;
            }
        }
        out.x = 0.5 * lplam * (1. + Math.cos(lpphi));
        out.y = lpphi;
        return out;
    }

    public Point2D.Double projectInverse(double xyx, double xyy, Point2D.Double out) {
        out.x = 2. * xyx / (1. + Math.cos(xyy));
        out.y = MapMath.asin(0.5 * (xyy + Math.sin(xyy)));
        return out;
    }

    public boolean hasInverse() {
        return true;
    }

    public String toString() {
        return "Nell";
    }
}
