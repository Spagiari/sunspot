/*
 * Copyright 2010-2011
 * Adilson Eduardo Spagiari (e.spagiari@gmail.com)
 * Israel Florentino (learsi@gmail.com)
 * Wander Lairson Costa aka walac (wander.lairson@gmail.com)
 *
 * This is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * version 2 of the License.
 * 
 * This software is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * source code. If not, see http://www.gnu.org/licenses/.
 */

#ifndef sunspot_SolarAnalyzerGroup_hpp
#define sunspot_SolarAnalyzerGroup_hpp

#include "SolarAnalyzer.hpp"

#define SUNRADIUS_6GRAUS 12
#define DELTA 3
#define THRESHOLD 2
#define PI 3.14159265

namespace SolarAnalyzer
{
    
    IplImage *circle_sunspots(IplImage *img, blob_collection &blobs);
    
    SolarStructure center_sun(IplImage *img, int debug);
    
    int count_groups(const SolarStructure &sun, blob_collection &blobs, IplImage *img, vSunspotGroup &vGroups);
    
    void put_group (vSunspotGroup &groups, unsigned int blobida, unsigned int blobidb);

}

#endif // UTILS_HPP