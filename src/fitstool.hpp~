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

#ifndef FITSTOOL_HPP
#define FITSTOOL_HPP

#define _USE_MATH_DEFINES

#include "fitsio.h"
#include "blob.hpp"
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <memory>

struct struct_sun
{
	coordinate center;
	int radius;
};

typedef std::map<int, double *> MAPFITS;

class fitsTool
{
protected:
    fitsfile *fptr = NULL;

    MAPFITS mapFitLines;

    long imgSize[3] = {1,1,1};

public:

    int openFit(const char * sFitFile);

    int sunPhysics(struct_sun * sun);

    int dateObs(std::string & sdate);

    int obsDist(double * dDist);

    int imScale(double * dImScale);

    bool validPoint(coordinate point);

    int readFitFile(coordinate pt);

    int readFitPoint(coordinate pt, double *dValue);

    int reposCenterBlob(blob *b);

    int closeFit();

    int IntensityAverage(const std::vector<coordinate> & vPoints, double *dIntAverage);
};

class fitsToolMDIContinum : public fitsTool
{
public:
    int openFit(const char * sFitFile);
    int regionGrowing(coordinate point, int *npoints, double *dIntensity, IplImage *src, const double dMin, std::vector<coordinate> &vInPoints);
    int sunMeanIntensity(double * dMean, double * dMin, double * dMax);

    int Area(const blob *b, double *AreaM2);
};

class fitsToolMag : public fitsTool
{
public:
    int openFit(const char * sFitFile);
    int magFieldAverage(const std::vector<coordinate> & vPoints, double *dMagAverage, double *dMagAveragePos, double *dMagAverageNeg);
};

#endif //FITSTOOL_HPP
