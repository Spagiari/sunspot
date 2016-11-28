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

#include "detect.hpp"
#include "blob.hpp"
#include "group.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <getopt.h>
#include <unistd.h>
#include "fitstool.hpp"
#include "highgui.h"
#include "cv.h"
#include "utils.hpp"
#include "fitstool.hpp"

int debug = 0;

int main(int argc, char **argv)
{
    create_window("Original");
    create_window("Result");

    cvMoveWindow("Original", 10, 200);
    cvMoveWindow("Result", 20 + WND_WIDTH, 200);

    const char *optstr = "d";
    int ch;

    while ((ch = getopt(argc, argv, optstr)) != EOF) {
        switch (ch)
        {
            case 'd':
                debug = 1;
                break;
        }
    }

    if ((argc - optind) > 2) {

        std::ofstream fWolf;
        std::ofstream fSpotCatalog;

        fWolf.open ("Wolfs.txt");
        fSpotCatalog.open ("SpotCatalog.txt");

        fSpotCatalog << "Image;Date;Spot Id;Lat;Lon;Area;Intensity Relative;Mag Avg;Mag Avg Pos;Mag Avg Neg; Mag Min; Mag Max; Sun Intensity Avg; Sun Intensity Max; Sun Intensity Min" << std::endl;

        for (int i = optind; i < argc; ++i) {
            IplImage *src = cvLoadImage(argv[i]);
            if (debug)
                cvShowImage("Original", src);

            std::cout << "Detect sun" << std::endl;
            IplImage *dst = detect_sunspots(src);

            std::cout << "Detect Blob" << std::endl;
            blob_collection b = detectBlobs(dst);

            std::cout << "Found " << b.size() << " sunspots in image " << argv[i] << std::endl;

            //struct_sun sun = center_sun(src, debug);

            fitsToolMDIContinum oFitMDIcontinum;
            fitsToolMag oFitMag;

            oFitMDIcontinum.openFit(argv[i]);

            struct_sun sun;

            oFitMDIcontinum.sunPhysics(&sun);

            std::string date;

            oFitMDIcontinum.dateObs(date);

            std::cout << "SD X0 = " << sun.center.x << " Y0 = " << sun.center.y << " Radius = " << sun.radius  << " Date:" << date << std::endl;

            group_sunspot_vector groups = count_groups(sun, b, debug? dst : NULL);

            std::cout << "groups found:" << groups.size() << std::endl;

            std::cout << "Wolf number " << b.size() + (groups.size() * 10) << " sunspots in image " << argv[i] << std::endl;

            double dIntensity;
            double dMin;
            double dMax;

            oFitMDIcontinum.sunMeanIntensity(&dIntensity, &dMin, &dMax);

            fWolf << argv[i] << ";"<< b.size() + (groups.size() * 10) << ";" << dMax << ";" << date << std::endl;

            oFitMag.openFit(argv[i]);

            blob_collection::iterator iter;

            int icnt = 0;

            for (iter = b.begin(); iter != b.end(); iter++)
            {
                double dMagAverage, dMagAveragePos, dMagAverageNeg, dArea, dIntensityAvg, dMagMin, dMagMax;

                oFitMag.magFieldAverage(iter->second.listPoints, &dMagAverage, &dMagAveragePos, &dMagAverageNeg, &dMagMin, &dMagMax);
                oFitMDIcontinum.Area(&(iter->second), &dArea);
                oFitMDIcontinum.IntensityAverage(iter->second.listPoints, &dIntensityAvg);

                std::cout  << "/t" << iter->first << "| Center:" << iter->second.getLatDeg() << "," << iter->second.getLonDeg() << "| Points : " << iter->second.listPoints.size() << "| AREA : " << dArea << "| Intensidade " << dIntensityAvg / dIntensity << "| Mag:" << dMagAverage << " : " << dMagAveragePos << " : " << dMagAverageNeg << std::endl;

                fSpotCatalog << argv[i] << ";" << date << ";" << icnt << ";" << iter->second.getLatDeg() << ";" << iter->second.getLonDeg() << ";" << dArea << ";" <<  dIntensityAvg / dIntensity << ";" << dMagAverage << ";" << dMagAveragePos << ";" << dMagAverageNeg << ";"<< dMagMin << ";" << dMagMax << ";" << dIntensity << ";" << dMax << ";" << dMin  << std::endl;;

                ++icnt;


            }

            oFitMDIcontinum.closeFit();
            oFitMag.closeFit();

            cvShowImage("Original", src);
            cvShowImage("Result", dst);

			/*double dIntensity;

			openFit(argv[i]);
			sunMeanIntensity(&dIntensity);

			std::cout << "Sun Intensity: " << dIntensity;
			*/

            oFitMDIcontinum.closeFit();

            cvWaitKey(5);
            cvReleaseImage(&src);
            cvReleaseImage(&dst);
        }
        fWolf.close();
        fSpotCatalog.close();

    } else {
        IplImage *src = cvLoadImage(argv[optind]);
        if (debug)
            cvShowImage("Original", src);
        IplImage *dst = detect_sunspots(src);

        blob_collection b = detectBlobs(dst);

        std::cout << "Found " << b.size() << " sunspots in image " << argv[optind] << std::endl;

        fitsToolMDIContinum oFitMDIcontinum;
        fitsToolMag oFitMag;

        //struct_sun sun = center_sun(src, debug);

		oFitMDIcontinum.openFit(argv[optind]);

		struct_sun sun;

		oFitMDIcontinum.sunPhysics(&sun);

        std::string date;

        oFitMDIcontinum.dateObs(date);

        std::cout << "SDO X0 = " << sun.center.x << " Y0 = " << sun.center.y << " Radius = " << sun.radius << " Date:" << date << std::endl;

        group_sunspot_vector groups = count_groups(sun, b, debug? dst : NULL);

        std::cout << "groups found:" << groups.size() << std::endl;

        std::cout << "Wolf number: " << b.size() + (groups.size() * 10) << " sunspots in image " << argv[optind] << std::endl;

        if (!debug) {
            cvShowImage("Original", src);
            cvShowImage("Result", dst);

			std::cout << "Src Width:" << src->width << "\t" << "Height:"<< src->height << std::endl;
			std::cout << "Dst Width:" << dst->width << "\t" << "Height:"<< dst->height << std::endl;

        }

		double dIntensity;
        double dMin;
        double dMax;

		oFitMDIcontinum.sunMeanIntensity(&dIntensity, &dMin, &dMax);

        dMin = dIntensity * 0.75;

		std::cout << "Sun Intensity: " << dIntensity << std::endl;

		blob_collection::iterator iter;
/*
		int npoints=0;
		double dspotInt=0;
        double dMagAverage = 0;
*/
        oFitMag.openFit(argv[optind]);


		for (iter = b.begin(); iter != b.end(); iter++)
		{
            double dMagAverage, dMagAveragePos, dMagAverageNeg, dArea, dIntensityAvg, dMagMin, dMagMax;

            oFitMag.magFieldAverage(iter->second.listPoints, &dMagAverage, &dMagAveragePos, &dMagAverageNeg, &dMagMin, &dMagMax);
            oFitMDIcontinum.Area(&(iter->second), &dArea);
            oFitMDIcontinum.IntensityAverage(iter->second.listPoints, &dIntensityAvg);

            std::cout << iter->first << "| Center:" << iter->second.getLatDeg() << "," << iter->second.getLonDeg() << "| Points : " << iter->second.listPoints.size() << "| AREA : " << dArea << "| Intensidade " << dIntensityAvg / dIntensity << "| Mag:" << dMagAverage << " : " << dMagAveragePos << " : " << dMagAverageNeg << std::endl;

            // intensidade mancha // center em graus // depois do primeiro teste 40 ou 50 // B/Int B/A Int/A | int média mancha / int central sol


/*
            iter->second.listPoints.clear();
			std::cout << iter->first << "| Center:" << iter->second.center.x << "," << iter->second.center.y;
			oFitMDIcontinum.reposCenterBlob(&(iter->second));
			std::cout << "| Re-Center:" << iter->second.center.x << "," << iter->second.center.y;
			oFitMDIcontinum.regionGrowing(iter->second.center, &npoints, &dspotInt, src, dMin, iter->second.listPoints);
			std::cout << " | Points:" << npoints << "\t| Inyensity: " << dspotInt << std::endl;
 */
		}

		oFitMDIcontinum.closeFit();
        oFitMag.closeFit();

        while (27 == (cvWaitKey(0) & 0xff));

        cvDestroyAllWindows();
        cvReleaseImage(&src);
        cvReleaseImage(&dst);

    }

    return 0;
}

