/*
Copyright 2009 Bernhard Jenny

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

package com.jhlabs.map.proj;

/**
 * Kavraisky 7 Projection. Based on Proj4.
 * @author Bernhard Jenny, Institute of Cartography, ETH Zurich
 */
public class Kavraisky7Projection extends EllipticalPseudoCylindrical {

    /** Creates a new instance of Kavraisky7Projection */
    public Kavraisky7Projection() {
        super(Math.sqrt(3.) / 2.,
                1.,
                0.,
                3. / Math.PI / Math.PI);
    }

    public String toString() {
        return "Kavraisky VII";
    }
}
