#include "Cell.h"

using namespace std;

extern string progname;
extern bool DEB;

Cell::Cell(string idp,float xp,float yp,int trcp,int cpcp,int cccp,int uccp,int dccp,int tcp,float areap,float nuc_areap,bool &err)
{
 id=idp;
 x=xp;
 y=yp;
 trc=trcp;
 cpc=cpcp;
 ccc=cccp;
 ucc=uccp;
 dcc=dccp;
 tc=tcp;
 area=areap;
 nuc_area=nuc_areap;
 
 // Some basic checks;
 if (x<0.0 || y<0.0)
 {
  WARNPROG("A cell coordinate x, or y, or both is/are negative.\n";)
  err=true;
 }
 if (trc<0 || cpc<0 || ccc<0 || ucc<0 || dcc<0 || tc<0)
 {
  WARNPROG("At least one of the counts is negative.\n";)
  err=true; 
 }
 if (trc+cpc+ccc+ucc+dcc!=tc)
 {
  WARNPROG("The sum of all types of counts, which is " <<  trc+cpc+ccc+ucc+dcc << " is not equal to the total number of counts, which is " << tc << ".\n";)
  err=true; 
 }
 if (area<0.0 || nuc_area<0.0)
 {
  WARNPROG("A cell area or cell nucleus area or both (whose values are " << area << " and " << nuc_area << " respectively) is/are negative.\n";)
  err=true;
 }
 if (nuc_area>area)
 {
  WARNPROG("The area of the cell nucleus, which is " << nuc_area << ", is bigger than the area of the cell itself which is " << area << ".\n";)
  err=true;
 }
}

Cell::Cell(InputBufferedFile &f,size_t &nbc,size_t &ncc)
{
 f.readb(id);
 f.readb(x);
 f.readb(y);
 f.readb(trc);
 f.readb(cpc);
 f.readb(ccc);
 f.readb(ucc);
 f.readb(dcc);
 f.readb(tc);
 f.readb(area);
 f.readb(nuc_area);
 f.readb(cluster);
 // These are vectors of a fundamental data type (float, or size_t for the last one) so readb is able to interpret and read them
 f.readb(cboundx);
 f.readb(cboundy);
 if (cboundx.size()!=cboundy.size())
  ERRPROG("Cell constructor from binary file: size of cell contour x and y array is different.\n";)
 if (cboundx.size()>0)
  nbc++;

 f.readb(nuccboundx);
 f.readb(nuccboundy);
 if (nuccboundx.size()!=nuccboundy.size())
  ERRPROG("Cell constructor from binary file: size of nucleus contour x and y array is different.\n";)
 if (nuccboundx.size()>0)
  ncc++;

 f.readb(startpoints);
}

void Cell::SaveAsBinary(OutputBufferedFile &f)
{
 f.writeb(id);
 f.writeb(x);
 f.writeb(y);
 f.writeb(trc);
 f.writeb(cpc);
 f.writeb(ccc);
 f.writeb(ucc);
 f.writeb(dcc);
 f.writeb(tc);
 f.writeb(area);
 f.writeb(nuc_area);
 f.writeb(cluster);
 // These are vectors of a fundamental data type (float, or size_t for the last one) so writeb is able to interpret and write them
 f.writeb(cboundx);
 f.writeb(cboundy);
 f.writeb(nuccboundx);
 f.writeb(nuccboundy);
 f.writeb(startpoints);
}

string Cell::GetAll(char sep,bool getcl)
{
 stringstream st;
 if (getcl)
  st << id << sep << x << sep << y << sep << trc << sep << cpc << sep << ccc << sep << ucc << sep << dcc << sep << tc << sep << area << sep << nuc_area << sep << cluster;
 else
  st << id << sep << x << sep << y << sep << trc << sep << cpc << sep << ccc << sep << ucc << sep << dcc << sep << tc << sep << area << sep << nuc_area;
 return(st.str());
}
