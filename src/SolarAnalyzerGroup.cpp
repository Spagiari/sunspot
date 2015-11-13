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
#include "SolarAnalyzerGroup.hpp"
#include "SolarAnalyzerBlob.hpp"
#include "SolarAnalyzerUtils.hpp"

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <cmath>
#include <cv.h>

namespace SolarAnalyzer
{
    using namespace std;
    
    typedef unsigned char pixel_type;
    
    IplImage *circle_sunspots(IplImage *img, blob_collection &blobs)
    {
        IplImage *dst = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 3);
        
        //dst = cvCloneImage( img );
        
        cvCvtColor( img, dst, CV_GRAY2RGB );
        
        blob_collection::const_iterator iter;
        
        for (iter = blobs.begin(); iter != blobs.end(); ++iter) 
        {
            std::cout << "x:" << iter->second.center.x << " y:" <<iter->second.center.y << std::endl;
            cvCircle(dst, cvPoint(iter->second.center.x, iter->second.center.y), SUNRADIUS_6GRAUS , cvScalar(0,255,0), 1);
            
            std::vector<coordinate>::const_iterator iterl;
            for (iterl = iter->second.listPoints.begin(); iterl != iter->second.listPoints.end(); ++iterl) 
            {
                std::cout << "\tx:" << iterl->x << " y:" <<iterl->y << std::endl;
            }
            
        }
        
        return dst;
    }
    
    void put_group (vSunspotGroup &groups, unsigned int blobida, unsigned int blobidb)
    {
        vSunspotGroup::iterator itera;
        vSunspotGroup::iterator iterb;
        vector<unsigned int>::size_type i;
        vector<unsigned int>::iterator p;
        bool finda = false, findb = false;
        
        for (itera = groups.begin(); itera != groups.end(); ++itera) 
        {
            p = find(itera->blobids.begin(), itera->blobids.end(), blobida);
            
            if (p != itera->blobids.end()) 
            {
                finda = true;
                break;
            }
        }
        
        for (iterb = groups.begin(); iterb != groups.end(); ++iterb) 
        {
            p = find(iterb->blobids.begin(), iterb->blobids.end(), blobidb);
            
            if (p != iterb->blobids.end()) 
            {
                findb = true;
                break;
            }
        }
        
        if (!finda && !findb) 
        {
            if (blobida == blobidb)
            {
                i = groups.size();
                SunspotGroup gs = {static_cast<int>(i)};
                groups.push_back(gs);
                groups[i].blobids.push_back(blobida);
            }
            else 
            {
                i = groups.size();
                SunspotGroup gs = {static_cast<int>(i)};
                groups.push_back(gs);
                groups[i].blobids.push_back(blobida);
                groups[i].blobids.push_back(blobidb);		
                
            }
        }
        
        if (finda && !findb) 
        {
            itera->blobids.push_back(blobidb);		
        }
        
        if (!finda && findb) 
        {
            iterb->blobids.push_back(blobida);		
        }
        
        if (finda && findb) 
        {
            if (itera->id_group != iterb->id_group)
            {
                vector<unsigned int>::iterator iterblob;
                
                for (iterblob = iterb->blobids.begin(); iterblob != iterb->blobids.end(); ++iterblob) 
                {
                    itera->blobids.push_back(*iterblob);
                }
                
                groups.erase(iterb);
                
            }
            
        }
        
    }
    
    int count_groups(const SolarStructure &sun, blob_collection &blobs, IplImage *img, vSunspotGroup &vGroups)
    {
        int x=0, y=0;
        
        blob_collection::iterator iter;
        
        for (iter = blobs.begin(); iter != blobs.end(); ++iter) 
        {
            x = iter->second.center.x - sun.center.x;
            y = iter->second.center.y - sun.center.y;
            
            double lat=asin(((double(y)/sun.radius) ) );
            double lon=asin(((double(x)/sqrt(( pow(double(sun.radius),2.0) - pow(y , 2.0)) ))));
            
            if (img != NULL)
                cout << iter->first << " lat:" << lat << "\tLon:" << lon << "Points: "<< iter->second.listPoints.size() << endl;
            
            iter->second.lat = lat;
            iter->second.lon = lon;
            
        }
        
        blob_collection::const_iterator citera;
        for (citera = blobs.begin(); citera != blobs.end(); ++citera) 
        {
            IplImage *dst;
            if (img != NULL) 
            {
                dst = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 3);
                cvCvtColor( img, dst, CV_GRAY2RGB );
                
                cvCircle(dst, cvPoint(citera->second.center.x, citera->second.center.y), SUNRADIUS_6GRAUS , cvScalar(0,255,0), 1);
                
                cvShowImage("Result",dst);
            }
            
            put_group(vGroups, citera->first, citera->first);
            
            blob_collection::const_iterator citerb;
            for (citerb = blobs.begin(); citerb != blobs.end(); ++citerb)
            {
                if (citera->first != citerb->first) 
                {				
                    double coslamb = sin(citera->second.lat) * sin(citerb->second.lat) + 
                    cos(citera->second.lat)*cos(citerb->second.lat) * 
                    cos(citerb->second.lon - citera->second.lon);
                    double lamb = acos(coslamb);
                    
                    if (lamb <= ((6.0/180.0)*PI)) 
                    {
                        put_group(vGroups, citera->first, citerb->first);
                        if (img !=NULL)
                            cvCircle(dst, cvPoint(citerb->second.center.x, citerb->second.center.y), SUNRADIUS_6GRAUS , cvScalar(255,0,0), 1);					
                    }
                    else
                    {
                        if (img !=NULL)
                            cvCircle(dst, cvPoint(citerb->second.center.x, citerb->second.center.y), SUNRADIUS_6GRAUS , cvScalar(255,255,255), 1);
                    }
                    if (img !=NULL)
                    {
                        cvShowImage("Result",dst);
                        cvWaitKey(10);
                    }
                }
            }
            if (img !=NULL)
                cvReleaseImage(&dst);
        }
        return E_OK;
    }
}