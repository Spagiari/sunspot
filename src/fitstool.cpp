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

#include "fitsTool.hpp"
#include <iostream>

#define QTYMEANPOINTS 20
#define DELTAINTENSITY 1.09


int fitsTool::openFit(const char * sFitFile)
{
	int status = 0;  /* MUST initialize status */
	
	if (fptr == NULL) 
	{
		fits_open_file(&fptr, sFitFile, READONLY, &status);
		
		if (status) /* ERROR open file */
		{
			FILE *report = fopen ("reperror.txt","w");
			
			fits_report_error(report,status);
			
			fclose(report);
			fptr = NULL;
		}
		else
		{
            //std::cout << "FIT Open:" << sFitFile << std::endl;
			fits_get_img_size(fptr, 3, imgSize, &status);
		}
		
		return status;
		
	}
	return 998;
}

int fitsToolMDIContinum::openFit(const char * sFitFile)
{
    std::string sFit(sFitFile);
    sFit.replace(sFit.length()-3,3,"fits");
    return fitsTool::openFit(sFit.c_str());
}

int fitsToolMag::openFit(const char * sFitFile)
{
    std::string sFit(sFitFile);
    sFit.replace(sFit.length()-4,4,"_mag.fits");
    return fitsTool::openFit(sFit.c_str());
}

int fitsTool::closeFit()
{
	int status = 0;
	
	if (fptr != NULL) 
	{
		if(fits_close_file(fptr, &status))
			return status;
		fptr = NULL;
        //std::cout << "FIT Closed" << std::endl;
		
		MAPFITS::iterator iter;
		
		for (iter = mapFitLines.begin(); iter != mapFitLines.end(); iter++) 
		{
			free(iter->second);
		}
		mapFitLines.clear();
	}
	return status;
}

int fitsToolMDIContinum::sunMeanIntensity(double * dMean, double * dMin, double * dMax)
{
	if (fptr == NULL) /* OPEN file first */ 
	{
		return 999;
	}
	struct_sun sun;
    
    *dMin = 1000000;
    *dMax = 0;
	
	int iRet = sunPhysics(&sun);
	int status = 0;
	
	if (iRet)
		return iRet;
	
	double * ptrPix;
	long firstpix[3] = {sun.center.x - (QTYMEANPOINTS/2),sun.center.y - (QTYMEANPOINTS/2),1};
	
	ptrPix = (double *) malloc(QTYMEANPOINTS * sizeof(double));
	
	double dSum = 0;
	
	for (int i = 0; i < QTYMEANPOINTS; ++i) 
	{
		if (fits_read_pix(fptr, TDOUBLE, firstpix, QTYMEANPOINTS, NULL, ptrPix,NULL, &status))
			printf("Couldn't read pix\n");
		else 
		{
			for (int j = 0; j < QTYMEANPOINTS; ++j) 
			{
				dSum += ptrPix[j];
                if (ptrPix[j] < *dMin)
                {
                    *dMin = ptrPix[j];
                }
                if (ptrPix[j] > *dMax)
                {
                    *dMax = ptrPix[j];
                }
			}
		}
		firstpix[1] = firstpix[1] + 1;
	}
	
	*dMean = dSum / (QTYMEANPOINTS * QTYMEANPOINTS);
	free(ptrPix);
	return status;
}

int fitsTool::sunPhysics(struct_sun * sun)
{
	int status = 0;  /* MUST initialize status */
	double center_y, center_x, radius;
	
	if (fptr == NULL) /* OPEN file first */ 
	{
		return 999;
	}
	/*
     SOHO CENTER_Y CENTER_X
     SDO  Y0_MP    X0_MP
    */
    
	if ( fits_read_key(fptr, TDOUBLE, "CENTER_Y", &center_y, NULL, &status) )
        return status;
	else
		sun->center.y = (unsigned int) center_y;
	
	if ( fits_read_key(fptr, TDOUBLE, "CENTER_X", &center_x, NULL, &status) )
        return status;
	else
		sun->center.x = (unsigned int) center_x;
	
	if ( fits_read_key(fptr, TDOUBLE, "R_SUN", &radius, NULL, &status) )
        return status;
	else
		sun->radius = (unsigned int) radius;
	
	
	return status;
	
}

int fitsTool::dateObs(std::string & sdate)
{
	int status = 0;  /* MUST initialize status */
	char date[120];
	
	if (fptr == NULL) /* OPEN file first */ 
	{
		return 999;
	}
	/*
     SOHO CENTER_Y CENTER_X
     SDO  Y0_MP    X0_MP
     */
    
	if ( fits_read_key(fptr, TSTRING, "T_OBS", date, NULL, &status) )
        return status;
	else
		sdate = std::string(date);
		
	
	return status;
	
}

int fitsTool::imScale(double * dImScale)
{
	int status = 0;  /* MUST initialize status */
	double dScale;
	
	if (fptr == NULL) /* OPEN file first */
	{
		return 999;
	}
	/*
     SOHO CENTER_Y CENTER_X
     SDO  Y0_MP    X0_MP
     */
    
	if ( fits_read_key(fptr, TDOUBLE, "IM_SCALE", &dScale, NULL, &status) )
        return status;
	else
		*dImScale = dScale;
	
	return status;
	
}

int fitsTool::obsDist(double * dDist)
{
	int status = 0;  /* MUST initialize status */
	double dTmp;
	
	if (fptr == NULL) /* OPEN file first */
	{
		return 999;
	}
	/*
     SOHO CENTER_Y CENTER_X
     SDO  Y0_MP    X0_MP
     */
    
	if ( fits_read_key(fptr, TDOUBLE, "OBS_DIST", &dTmp, NULL, &status) )
        return status;
	else
		*dDist = dTmp;
	
	return status;
	
}

int fitsTool::readFitPoint(coordinate pt, double *dValue)
{
	int status = 0;
	
	pt.y = 1024 - pt.y;
	
	status = readFitFile(pt);
	
	if (status) 
		return status;
	
	*dValue =  mapFitLines[pt.y][pt.x];
	
	return 0;
}

int fitsTool::readFitFile(coordinate pt)
{
	MAPFITS::iterator iter;
	iter = mapFitLines.find(pt.y);
	
	if (iter == mapFitLines.end()) 
	{
		long firstpix[3] = {1,pt.y,1};
		int status = 0;
		
		mapFitLines.insert(MAPFITS::value_type(pt.y, NULL));
		mapFitLines[pt.y] = (double *) malloc(imgSize[0] * sizeof(double));
		
		fits_read_pix(fptr, TDOUBLE, firstpix, imgSize[0], NULL, mapFitLines[pt.y],NULL, &status);
		return status;
	}
	return 0;
}

int fitsToolMDIContinum::regionGrowing(coordinate point, int *npoints, double *dIntensity, IplImage *src, const double dMin, std::vector<coordinate> &vInPoints)
{	
	if (fptr == NULL) /* OPEN file first */ 
	{
		return 999;
	}
	
	point.y = 1024 - point.y;
	
	std::vector<coordinate> vExploitPoints;
//	std::vector<coordinate> vInPoints;
	std::vector<coordinate> vOutPoints;
	
	int status = 0;
	double dSum = 0;
	coordinate CandidatePoint;
	
	readFitFile(point);
	
	double dValueBase = mapFitLines[point.y][point.x] * DELTAINTENSITY;
    
    if (dValueBase > dMin)
        dValueBase = dMin;
	
	vExploitPoints.push_back(point);
	
	std::map<int, double *>::iterator iter;
	
	while (vExploitPoints.size() > 0) 
	{

		status = readFitFile(vExploitPoints[0]);
		
		if (status)
			return status;
		
		if (mapFitLines[vExploitPoints[0].y][vExploitPoints[0].x] <= dValueBase)
		{
			vInPoints.push_back(vExploitPoints[0]);
			
			dSum += mapFitLines[vExploitPoints[0].y][vExploitPoints[0].x];
			
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

bool fitsTool::validPoint(coordinate point)
{
	if (fptr == NULL) /* OPEN file first */ 
	{
		return false;
	}
	
	return point.x > 0 && point.x < imgSize[0] && point.y > 0 && point.y < imgSize[1];
}

int fitsTool::reposCenterBlob(blob *b)
{
	double dIntensity = 0;
	double dtmp = 0;
	coordinate ptMinInt;
	
	std::vector<coordinate>::iterator inter;
	
	readFitPoint(b->listPoints[0], &dIntensity);
	
    ptMinInt.x = b->listPoints[0].x;
    ptMinInt.y = b->listPoints[0].y;
	
	for (inter = b->listPoints.begin(); inter != b->listPoints.end(); inter++) 
	{
		int status = readFitPoint((*inter), &dtmp);
		if (status)
			return status;
		if (dtmp < dIntensity) 
		{
			dIntensity = dtmp;
			ptMinInt.x = inter->x;
            ptMinInt.y = inter->y;
		}
	}
	b->center.x = ptMinInt.x;
    b->center.y = ptMinInt.y;
	return 0;
}

int fitsToolMDIContinum::Area(const blob *b, double *AreaM2)
{
    double dObsDist;
    double dTmp;
    double dTeta;
    double dP;
    
    int status = 0;
    status = obsDist(&dObsDist);
    
    if (status)
        return status;
    
    status = imScale(&dTeta);
    
    if (status)
        return status;
    
    dTmp =  (dTeta / 3600) * (M_PI / 180) * (150000000.0 * dObsDist);
    
    dP = dTmp / (cos(b->lat) * cos(b->lon));
    
    *AreaM2 = b->listPoints.size() * pow( dP, 2.0 );
    
    return 0;
    
}

int fitsToolMag::magFieldAverage(const std::vector<coordinate> & vPoints, double *dMagAverage, double *dMagAveragePos, double *dMagAverageNeg)
{
    double dtmp = 0;
    double sum = 0;
    double sumpos = 0;
    double sumneg = 0;
    int i = 0;
    int ipos = 0;
    int ineg = 0;
    
    std::vector<coordinate>::const_iterator inter;
    
    for (inter = vPoints.cbegin(); inter != vPoints.cend(); inter++)
	{
        int status = readFitPoint((*inter), &dtmp);
		if (status)
			return status;
        
        sum += dtmp;
        
        if (dtmp >= 0)
        {
            sumpos += dtmp;
            ++ipos;
        }
        else
        {
            sumneg += dtmp;
            ++ineg;
        }
        
        ++i;

    }
    
    *dMagAverage = sum / i;
    *dMagAveragePos = sumpos / ipos;
    *dMagAverageNeg = sumneg / ineg;
    
    return 0;
}

int fitsTool::IntensityAverage(const std::vector<coordinate> & vPoints, double *dIntAverage)
{
    double dtmp = 0;
    double sum = 0;
    int i = 0;

    
    std::vector<coordinate>::const_iterator inter;
    
    for (inter = vPoints.cbegin(); inter != vPoints.cend(); inter++)
	{
        int status = readFitPoint((*inter), &dtmp);
		if (status)
			return status;
        
        sum += dtmp;
        
        ++i;
        
    }
    
    *dIntAverage = sum / i;

    return 0;

}
