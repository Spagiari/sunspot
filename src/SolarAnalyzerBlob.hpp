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

#ifndef sunspot_SolarAnalyzerBlob_hpp
#define sunspot_SolarAnalyzerBlob_hpp

#include <map>
#include <cv.h>

#include "SolarAnalyzer.hpp"

namespace SolarAnalyzer
{
    int DetectBlobs(IplImage* frame, SolarAnalyzer::blob_collection & blobs);
}

#endif // BLOB_HPP
