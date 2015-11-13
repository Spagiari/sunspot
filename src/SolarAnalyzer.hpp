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

#ifndef SolarAnalyzer_hpp
#define SolarAnalyzer_hpp

#include <string>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <memory>
#include <fitsio.h>
#include <cv.h>
#include <boost/shared_array.hpp>
#include <boost/noncopyable.hpp>

#define QTYMEANPOINTS 20
#define DELTAINTENSITY 1.000001

namespace SolarAnalyzer
{
    typedef std::map<int, boost::shared_array<double> > MAPFITS;
    typedef std::pair<int, boost::shared_array<double> > MAPFITSLINE;

    class coordinate
    {
    public:
        void * data;
        unsigned int x, y;
        coordinate()
        :x(0)
        ,y(0)
        {
        };
        coordinate(int ix, int iy)
        :x(ix)
        ,y(iy)
        {
        };
        coordinate(const coordinate & cc)
        :x(cc.x)
        ,y(cc.y)
        {
        };
        bool operator==(const coordinate &other) const
        {
            return this->x == other.x && this->y == other.y;
        };
    };

    struct blob
    {
        coordinate min, max;

        coordinate center;

        double lat;
        double lon;

        std::vector<coordinate> listPoints;

        blob(const coordinate & cmin, const coordinate & cmax)
        :min(cmin)
        ,max(cmax)
        ,center(0,0)
        ,lat(0)
        ,lon(0)
        {
        };

        blob()
        {
        };

    };

    struct lineBlob
    {
        unsigned int min, max;
        unsigned int blobId;
        bool attached;
    };

    typedef std::map<unsigned int, blob> blob_collection;

    typedef enum {  E_FIT_FILE_ALREADY_OPEN = 998,
                    E_FIT_FILE_NOT_OPEN = 999,
                    E_READ_PIX = 1001,
                    E_JPG_FILE_ALREADY_OPEN = 1002,
                    E_JPG_FILE_NOT_OPEN = 1003,
                    E_SUNSPOT_DETECT = 1004,
                    E_SUNSPOT_BLOB  = 1005,
                    E_SUNSPOT_MEAN_INTENSITY = 1006,
                    E_OK = 0
    } SAReturn;

    typedef enum {  X = 0,
                    Y = 1,
                    Z = 2
    } SANaxis;

    struct SolarStructure
    {
        coordinate center;
        int radius;
    };

    struct SunspotGroup {
        int id_group;
        std::vector<unsigned int>  blobids;
    };

    typedef std::vector<SunspotGroup> vSunspotGroup;

    class cSunImageIMPL
    {
    public:
        std::string sFileName;

        std::shared_ptr<fitsfile>   ptrFitSource;
        std::shared_ptr<IplImage>   ptrJpgSource;
        std::shared_ptr<IplImage>   ptrJpgCustom;
        std::shared_ptr<IplImage>   ptrJpgTarget;

        MAPFITS                     mapFitLines;
        SolarStructure              sSun;
        long                        imgSize[3];

        blob_collection             bBlobs;
        double                      dMeanIntensity;

        int                         iDebug;

        cSunImageIMPL(const std::string & sFn, int Debug)
        :sFileName(sFn)
        ,iDebug(Debug)
        //,imgSize({1,1,1})
        {
            imgSize[0] = 1;
            imgSize[1] = 1;
            imgSize[2] = 1;
        };
    };

    class cSolarAnalyzer : public boost::noncopyable
    {
    private:
        std::shared_ptr<SolarAnalyzer::cSunImageIMPL> pSunInfo;
        int reposCenterBlob(blob *b);
        int openFit();
        int openJpg();
        int MeanIntensity   (double * dMean);
        int readFitPoint    (const coordinate & pt, double *dValue);
        int readJpgPoint    (const coordinate & pt, double *dValue);
        int readFitFile     (const coordinate & pt);
        int regionGrowing   (const coordinate & point, int *npoints, double *dIntensity, IplImage *src);
        bool validPoint     (const coordinate & point) const;
        static int closeFit(fitsfile *ptrFitFile);
        static int closeJpg(IplImage *ptrJpgFile);

    public:
        explicit cSolarAnalyzer(const std::string & sFileName, int Debug) throw(int);
        ~cSolarAnalyzer(){};
        int Physics(SolarStructure & sun);
        int ProcessImages();
        int GetWolfNumber();
    };

    static const std::shared_ptr<fitsfile> & null_fitfile_ptr()
    {
        static std::shared_ptr<fitsfile> onull;
        return onull;
    }
    static const std::shared_ptr<IplImage> & null_Jpgfile_ptr()
    {
        static std::shared_ptr<IplImage> onull;
        return onull;
    }
}

#endif
