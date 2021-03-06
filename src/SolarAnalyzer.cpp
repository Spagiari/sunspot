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

#include "SolarAnalyzer.hpp"
#include "SolarAnalyzerDetect.hpp"
#include "SolarAnalyzerBlob.hpp"

#include <cv.h>
#include <highgui.h>

namespace SolarAnalyzer
{
    cSolarAnalyzer::cSolarAnalyzer(const std::string & sFileName, int Debug) throw(int)
    :pSunInfo(new cSunImageIMPL(sFileName, Debug))
    {
        int iFitRet = openFit();
        if (iFitRet != E_OK)
            throw iFitRet;

        int iJpgRet = openJpg();
        if (iJpgRet != E_OK)
            throw iJpgRet;
    }

    int cSolarAnalyzer::ProcessImages()
    {
        Physics(pSunInfo->sSun);

        IplImage *src = detect_sunspots(pSunInfo->ptrJpgCustom.get(), pSunInfo->iDebug);

        if(src == nullptr)
            return E_SUNSPOT_DETECT;

        pSunInfo->ptrJpgCustom = std::shared_ptr<IplImage>(src,closeJpg);

        if(DetectBlobs(pSunInfo->ptrJpgCustom.get(), pSunInfo->bBlobs) != E_OK)
            return E_SUNSPOT_BLOB;

        if(MeanIntensity(&pSunInfo->dMeanIntensity) != E_OK)
           return E_SUNSPOT_MEAN_INTENSITY;

        return E_OK;

    }

    int cSolarAnalyzer::openFit()
    {
        int status = 0;  /* MUST initialize status */

        if (pSunInfo->ptrFitSource == null_fitfile_ptr())
        {
            std::string sFit(pSunInfo->sFileName);
            sFit.replace(sFit.length()-3,3,"fits");

            fitsfile *tmpptr;

            fits_open_file(&tmpptr, sFit.c_str(), READONLY, &status);

            if (!status)
            {
                pSunInfo->ptrFitSource = std::shared_ptr<fitsfile>(tmpptr, closeFit);
                fits_get_img_size(pSunInfo->ptrFitSource.get(), 3, pSunInfo->imgSize, &status);
            }

            return status;

        }
        return SAReturn::E_FIT_FILE_ALREADY_OPEN;
    }

    int cSolarAnalyzer::openJpg()
    {
        if (pSunInfo->ptrJpgSource == null_Jpgfile_ptr())
        {
            std::string sJpg(pSunInfo->sFileName);
            sJpg.replace(sJpg.length()-3,3,"fits");

            IplImage *src = cvLoadImage(sJpg.c_str());

            if (src == nullptr)
            return E_JPG_FILE_NOT_OPEN;

            pSunInfo->ptrJpgSource = std::shared_ptr<IplImage>(src, closeJpg);

            return E_OK;
        }
        return E_JPG_FILE_ALREADY_OPEN;
    }

    int cSolarAnalyzer::closeFit(fitsfile *ptrFitFile)
    {
        int status = 0;

        if (ptrFitFile != nullptr)
        {
            if(fits_close_file(ptrFitFile, &status))
                return status;
            ptrFitFile = nullptr;

            return E_OK;

        }
        return SAReturn::E_FIT_FILE_NOT_OPEN;
    }

    int cSolarAnalyzer::closeJpg(IplImage *ptrJpgFile)
    {
        if (ptrJpgFile != nullptr)
        {
            cvReleaseImage(&ptrJpgFile);
            ptrJpgFile = nullptr;

            return E_OK;
        }
        return SAReturn::E_JPG_FILE_NOT_OPEN;
    }

    int cSolarAnalyzer::Physics(SolarStructure & sun)
    {
        int status = 0;  /* MUST initialize status */
        double center_y, center_x, radius;

        if (pSunInfo->ptrFitSource == null_fitfile_ptr()) /* OPEN file first */
        {
            return SAReturn::E_FIT_FILE_NOT_OPEN;
        }

        if ( fits_read_key(pSunInfo->ptrFitSource.get(), TDOUBLE, "CENTER_Y", &center_y, NULL, &status) )
            return status;
        else
            sun.center.y = (unsigned int) center_y;

        if ( fits_read_key(pSunInfo->ptrFitSource.get(), TDOUBLE, "CENTER_X", &center_x, NULL, &status) )
            return status;
        else
            sun.center.x = (unsigned int) center_x;

        if ( fits_read_key(pSunInfo->ptrFitSource.get(), TDOUBLE, "R_SUN", &radius, NULL, &status) )
            return status;
        else
            sun.radius = (unsigned int) radius;

        return status;

    }

    int cSolarAnalyzer::MeanIntensity(double * dMean)
    {

        double dSum = 0;
        double dTmp = 0;

        for (int i = 0; i < QTYMEANPOINTS; ++i)
        {
            for (int j = 0; j < QTYMEANPOINTS; ++j)
            {
                int err = readFitPoint(coordinate(i,j), &dTmp);
                if (err != E_OK)
                    return err;
                dSum += dTmp;
            }
        }
        *dMean = dSum / (QTYMEANPOINTS * QTYMEANPOINTS);
        return E_OK;

    }

    bool cSolarAnalyzer::validPoint(const coordinate & point) const
    {
        if (pSunInfo->ptrFitSource == null_fitfile_ptr()) /* OPEN file first */
        {
            return false;
        }

        return point.x > 0 && point.x < pSunInfo->imgSize[SANaxis::X] && point.y > 0 && point.y < pSunInfo->imgSize[SANaxis::Y];
    }


    int cSolarAnalyzer::readFitPoint(const coordinate & pt1, double *dValue)
    {
        int status = 0;

        coordinate pt(pt1);

        pt.y = pSunInfo->imgSize[SANaxis::Y] - pt.y;

        status = readFitFile(pt);

        if (status)
            return status;

        *dValue =  pSunInfo->mapFitLines[pt.y][pt.x];

        return 0;
    }

    int cSolarAnalyzer::readJpgPoint(const coordinate & pt, double *dValue)
    {
        if (pSunInfo->ptrJpgSource == null_Jpgfile_ptr())
        {
            return E_JPG_FILE_NOT_OPEN;
        }
        *dValue = (double) pSunInfo->ptrJpgSource.get()->imageData[(pt.x*(pSunInfo->ptrJpgSource.get()->width))+ pt.y];
        return E_OK;
    }

    int cSolarAnalyzer::readFitFile(const coordinate & pt)
    {
        if (pSunInfo->ptrFitSource == null_fitfile_ptr()) /* OPEN file first */
        {
            return SAReturn::E_FIT_FILE_NOT_OPEN;
        }

        MAPFITS::iterator iter;
        iter = pSunInfo->mapFitLines.find(pt.y);

        if (iter == pSunInfo->mapFitLines.end())
        {
            long firstpix[3] = {1,pt.y,1};
            int status = SAReturn::E_OK;

            pSunInfo->mapFitLines.insert(MAPFITSLINE(pt.y, boost::shared_array<double>(new double[pSunInfo->imgSize[SANaxis::X]])));

            fits_read_pix(pSunInfo->ptrFitSource.get(), TDOUBLE, firstpix, pSunInfo->imgSize[SANaxis::X], NULL, pSunInfo->mapFitLines[pt.y].get(),NULL, &status);
            return status;
        }
        return SAReturn::E_OK;
    }

    int cSolarAnalyzer::regionGrowing(const coordinate & point, int *npoints, double *dIntensity, IplImage *src)
    {
        if (pSunInfo->ptrFitSource == null_fitfile_ptr()) /* OPEN file first */
        {
            return SAReturn::E_FIT_FILE_NOT_OPEN;
        }

        std::vector<coordinate> vExploitPoints;
        std::vector<coordinate> vInPoints;
        std::vector<coordinate> vOutPoints;

        int status = 0;
        double dSum = 0;
        coordinate CandidatePoint;

        readFitFile(point);

        double dValueBase = pSunInfo->mapFitLines[point.y][point.x] * DELTAINTENSITY;

        vExploitPoints.push_back(point);

        std::map<int, double *>::iterator iter;

        while (vExploitPoints.size() > 0)
        {

            status = readFitFile(vExploitPoints[0]);

            if (status)
                return status;

            if (pSunInfo->mapFitLines[vExploitPoints[0].y][vExploitPoints[0].x] < dValueBase)
            {
                vInPoints.push_back(vExploitPoints[0]);

                dSum += pSunInfo->mapFitLines[vExploitPoints[0].y][vExploitPoints[0].x];

                if (src != NULL)
                {
                    CvPoint drawP = cvPoint(vExploitPoints[0].x, 1024 - vExploitPoints[0].y);
                    cvCircle(src, drawP, 1, CV_RGB(255,255,255), 1, CV_AA, 0);
                    cvShowImage("Result",src);
                }

                //Generate Candidate Points

                //superior

                CandidatePoint.x = vExploitPoints[0].x - 1; CandidatePoint.y = vExploitPoints[0].y - 1;

                if (validPoint(CandidatePoint))
                {
                    if (std::find(vExploitPoints.begin(), vExploitPoints.end(), CandidatePoint) == vExploitPoints.end() &&
                        std::find(vInPoints.begin(), vInPoints.end(), CandidatePoint) == vInPoints.end() &&
                        std::find(vOutPoints.begin(), vOutPoints.end(), CandidatePoint) == vOutPoints.end())
                    {
                        vExploitPoints.push_back(CandidatePoint);
                    }
                }

                CandidatePoint.x = vExploitPoints[0].x + 0; CandidatePoint.y = vExploitPoints[0].y - 1;

                if (validPoint(CandidatePoint))
                {
                    if (std::find(vExploitPoints.begin(), vExploitPoints.end(), CandidatePoint) == vExploitPoints.end() &&
                        std::find(vInPoints.begin(), vInPoints.end(), CandidatePoint) == vInPoints.end() &&
                        std::find(vOutPoints.begin(), vOutPoints.end(), CandidatePoint) == vOutPoints.end())
                    {
                        vExploitPoints.push_back(CandidatePoint);
                    }
                }

                CandidatePoint.x = vExploitPoints[0].x + 1; CandidatePoint.y = vExploitPoints[0].y - 1;

                if (validPoint(CandidatePoint))
                {
                    if (std::find(vExploitPoints.begin(), vExploitPoints.end(), CandidatePoint) == vExploitPoints.end() &&
                        std::find(vInPoints.begin(), vInPoints.end(), CandidatePoint) == vInPoints.end() &&
                        std::find(vOutPoints.begin(), vOutPoints.end(), CandidatePoint) == vOutPoints.end())
                    {
                        vExploitPoints.push_back(CandidatePoint);
                    }
                }

                //inferior
                CandidatePoint.x = vExploitPoints[0].x - 1; CandidatePoint.y = vExploitPoints[0].y + 1;

                if (validPoint(CandidatePoint))
                {
                    if (std::find(vExploitPoints.begin(), vExploitPoints.end(), CandidatePoint) == vExploitPoints.end() &&
                        std::find(vInPoints.begin(), vInPoints.end(), CandidatePoint) == vInPoints.end() &&
                        std::find(vOutPoints.begin(), vOutPoints.end(), CandidatePoint) == vOutPoints.end())
                    {
                        vExploitPoints.push_back(CandidatePoint);
                    }
                }

                CandidatePoint.x = vExploitPoints[0].x + 0; CandidatePoint.y = vExploitPoints[0].y + 1;

                if (validPoint(CandidatePoint))
                {
                    if (std::find(vExploitPoints.begin(), vExploitPoints.end(), CandidatePoint) == vExploitPoints.end() &&
                        std::find(vInPoints.begin(), vInPoints.end(), CandidatePoint) == vInPoints.end() &&
                        std::find(vOutPoints.begin(), vOutPoints.end(), CandidatePoint) == vOutPoints.end())
                    {
                        vExploitPoints.push_back(CandidatePoint);
                    }
                }

                CandidatePoint.x = vExploitPoints[0].x + 1; CandidatePoint.y = vExploitPoints[0].y + 1;

                if (validPoint(CandidatePoint))
                {
                    if (std::find(vExploitPoints.begin(), vExploitPoints.end(), CandidatePoint) == vExploitPoints.end() &&
                        std::find(vInPoints.begin(), vInPoints.end(), CandidatePoint) == vInPoints.end() &&
                        std::find(vOutPoints.begin(), vOutPoints.end(), CandidatePoint) == vOutPoints.end())
                    {
                        vExploitPoints.push_back(CandidatePoint);
                    }
                }

                //direita - esquerda
                CandidatePoint.x = vExploitPoints[0].x - 1; CandidatePoint.y = vExploitPoints[0].y + 0;

                if (validPoint(CandidatePoint))
                {
                    if (std::find(vExploitPoints.begin(), vExploitPoints.end(), CandidatePoint) == vExploitPoints.end() &&
                        std::find(vInPoints.begin(), vInPoints.end(), CandidatePoint) == vInPoints.end() &&
                        std::find(vOutPoints.begin(), vOutPoints.end(), CandidatePoint) == vOutPoints.end())
                    {
                        vExploitPoints.push_back(CandidatePoint);
                    }
                }

                CandidatePoint.x = vExploitPoints[0].x + 1; CandidatePoint.y = vExploitPoints[0].y + 0;

                if (validPoint(CandidatePoint))
                {
                    if (std::find(vExploitPoints.begin(), vExploitPoints.end(), CandidatePoint) == vExploitPoints.end() &&
                        std::find(vInPoints.begin(), vInPoints.end(), CandidatePoint) == vInPoints.end() &&
                        std::find(vOutPoints.begin(), vOutPoints.end(), CandidatePoint) == vOutPoints.end())
                    {
                        vExploitPoints.push_back(CandidatePoint);
                    }
                }

            }
            else
                vOutPoints.push_back(vExploitPoints[0]);


            vExploitPoints.erase(vExploitPoints.begin());

        }

        *npoints = vInPoints.size();
        *dIntensity = dSum / (*npoints);

        return 0;

    }

    int cSolarAnalyzer::reposCenterBlob(blob *b)
    {
        double dIntensity = 0;
        double dtmp = 0;
        coordinate ptMinInt;

        std::vector<coordinate>::iterator inter;

        readFitPoint(b->listPoints[0], &dIntensity);
        ptMinInt.x = b->listPoints[0].x; ptMinInt.y = b->listPoints[0].y;

        for (inter = b->listPoints.begin(); inter != b->listPoints.end(); inter++)
        {
            int status = readFitPoint((*inter), &dtmp);
            if (status)
                return status;
            if (dtmp < dIntensity)
            {
                dIntensity = dtmp;
                ptMinInt.x = inter->x; ptMinInt.y = inter->y;
            }
        }
        b->center.x = ptMinInt.x; b->center.y = ptMinInt.y;
        return 0;
    }

}
