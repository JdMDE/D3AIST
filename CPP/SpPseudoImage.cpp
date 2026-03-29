#include "SpPseudoImage.h"

// @cond DO_NOT_DOCUMENT

#include <iostream>
#include <string>

using namespace std;

extern string progname;
extern bool DEB;

SpPseudoImage::SpPseudoImage(unsigned int width,unsigned int height,Spacing rsp,My3DPoint origin,int cr)
{
  w=width;
  h=height;
  sp[0]=rsp[0];
  sp[1]=rsp[1];
  sp[2]=rsp[2];
  orig[0]=origin[0];
  orig[1]=origin[1];
  orig[2]=0.0;
  corner[0]=orig[0]+w*sp[0];
  corner[1]=orig[1]+h*sp[1];
  coderest=cr;
  if (DEB)
  {
   cout << "Allocating memory for pseudoimage. Please, wait... ";
   cout.flush();
  }
  ima = new PixContentPointer[h];
  for (size_t r=0; r<h; r++)
   ima[r] = new PixContent[w];
  if (DEB)
  {
   cout << "Done.\n";
   cout.flush();
  }
}

SpPseudoImage::~SpPseudoImage()
{
 if (DEB)
 {
  cout << "Deallocating memory for pseudoimage... ";
  cout.flush();
 }
 if (ima!=nullptr)
 {
  for (size_t r=0; r<h; r++)
  {
   if (ima[r]!=nullptr)
    delete[] ima[r];
  }
  delete[] ima;
 }
 if (DEB)
 {
  cout << "Done!\n";
  cout.flush();
 }
}

bool SpPseudoImage::PPhysicalPointToIndex(My3DPoint &p,My3DIndex &idx)
{
 float dx=p[0]-orig[0];
 float dy=p[1]-orig[1];

 int idx0=int(dx/sp[0]);
 int idx1=h-int(dy/sp[1]);
 // This is to take into account the different orientation. Found after a unit test.
 if (idx1>0)
  idx1--;
 if (idx0>=0 && (unsigned int)idx0<w && idx1>=0 && (unsigned int)idx1<h)
 {
  // WARNING: our index convention uses row (vertical) index before column (horizontal) index
  idx[1]=idx0;
  idx[0]=idx1;
  idx[2]=0;
  return(true);
 }
 else
  return(false);
}

bool SpPseudoImage::PIndexToPhysicalPoint(unsigned int r,unsigned int c,MyPoint &p)
{
 p[0]=orig[0]+(float(c)+0.5)*sp[0];
 p[1]=orig[1]+(h-float(r)-0.5)*sp[1];
 return (p[0]>=0.0 && p[0]<=corner[0] && p[1]>=0.0 && p[1]<=corner[1]);
}

void SpPseudoImage::TransfTest()
{
 My3DIndex itest;
 itest[2]=0;
 My3DPoint ptest3;
 MyPoint ptest2;
 for (size_t r=0; r<h; r+=(h/1000))
  for (size_t c=0; c<w; c+=(w/1000))
  {
   cout << "[" << r << "," << c << "] --> ";
   if (PIndexToPhysicalPoint(r,c,ptest2))
   {
    cout << "(" << ptest2[0] << "," << ptest2[1] << ") --> ";
    ptest3[0]=ptest2[0];
    ptest3[1]=ptest2[1];
    ptest3[2]=0.0;
    if (PPhysicalPointToIndex(ptest3,itest))
     cout << "[" << itest[0] << "," << itest[1] << "]\n";
    else
     cout << "INVALID\n";
   }
   else
    cout << "INVALID\n";
  }
 cout << "Limit cases:\n";

 cout << "[" << 0 << "," << w-1 << "] --> ";
 if (PIndexToPhysicalPoint(0,w-1,ptest2))
 {
  cout << "(" << ptest2[0] << "," << ptest2[1] << ") --> ";
  ptest3[0]=ptest2[0];
  ptest3[1]=ptest2[1];
  ptest3[2]=0.0;
  if (PPhysicalPointToIndex(ptest3,itest))
   cout << "[" << itest[0] << "," << itest[1] << "]\n";
  else
  cout << "INVALID\n";
 }
 else
  cout << "INVALID\n";

 cout << "[" << h-1 << "," << 0 << "] --> ";
 if (PIndexToPhysicalPoint(h-1,0,ptest2))
 {
  cout << "(" << ptest2[0] << "," << ptest2[1] << ") --> ";
  ptest3[0]=ptest2[0];
  ptest3[1]=ptest2[1];
  ptest3[2]=0.0;
  if (PPhysicalPointToIndex(ptest3,itest))
   cout << "[" << itest[0] << "," << itest[1] << "]\n";
  else
  cout << "INVALID\n";
 }
 else
  cout << "INVALID\n";

 cout << "[" << h-1 << "," << w-1 << "] --> ";
 if (PIndexToPhysicalPoint(h-1,w-1,ptest2))
 {
  cout << "(" << ptest2[0] << "," << ptest2[1] << ") --> ";
  ptest3[0]=ptest2[0];
  ptest3[1]=ptest2[1];
  ptest3[2]=0.0;
  if (PPhysicalPointToIndex(ptest3,itest))
   cout << "[" << itest[0] << "," << itest[1] << "]\n";
  else
  cout << "INVALID\n";
 }
 else
  cout << "INVALID\n";
}

void SpPseudoImage::SimpleInsertOrIncrementPixelValue(My3DIndex &idx)
{
 PixContent::iterator itb=ima[idx[0]][idx[1]].begin();
 PixContent::iterator ite=ima[idx[0]][idx[1]].end();

 // If this pixel/location has not yet any event, just create a new list with this event, store it and return. Value is 1 since now there is just one item.
 // The Realcoord list is not altered (always empty) since in this case, locations will be in all cases the center of the corresponding pseudopixel.
 if (itb==ite)
 {
  entry e;
  e.code=idx[2];
  e.v=1;
  ima[idx[0]][idx[1]].push_back(e);
  return;
 }

 // This pixel has already something at it. Let's find if any of the events present come from the same gene (have the same code)
 PixContent::iterator itl=itb;
 while (itl!=ite && itl->code!=int(idx[2]))
  ++itl;

 // If not (we are at the end of the list and codes were different) add the new event at the end of those already present.
 if (itl==ite)
 {
  entry e;
  e.code=idx[2];
  e.v=1;
  ima[idx[0]][idx[1]].push_back(e);
 }
 // Otherwise, we have the same type of event at the same pseudopixel. Just increment the count.
 else
  itl->v++;
}

void SpPseudoImage::InsertOrIncrementPixelValue(My3DIndex &idx,float dnuc,Cell &cell,size_t &rep_locations,size_t &rep_locations_and_codes,float x,float y)
{
 PixContent::iterator itb=ima[idx[0]][idx[1]].begin();
 PixContent::iterator ite=ima[idx[0]][idx[1]].end();

 // If this pixel/location has not yet any event, just create a new list with this event, store it and return. Value is 1 since now there is just one item.
 if (itb==ite)
 {
  entry e;
  e.code=idx[2];
  e.Realcoord.push_back(pair<float,float>(x,y));
  e.DistNuc.push_back(dnuc);
  e.vcells.push_back(cell);
  //e.AreaRatios.push_back(arearatio);
  //e.AttClust.push_back(attclus);
  e.v=1;
  ima[idx[0]][idx[1]].push_back(e);
  return;
 }

 // This pixel has already something at it. Let's find if any of the events present come from the same gene (have the same code)
 PixContent::iterator itl=itb;
 while (itl!=ite && itl->code!=int(idx[2]))
   ++itl;

 // If not (we are at the end of the list and codes were different) add the new event at the end of those already present.
 // Notice that the location COULD be repeated (two _different_ genes expressing just at the very same point)
 if (itl==ite)
 {
  entry e;
  e.code=idx[2];
  e.Realcoord.push_back(pair<float,float>(x,y));
  e.DistNuc.push_back(dnuc);
  e.vcells.push_back(cell);
  e.v=1;
  ima[idx[0]][idx[1]].push_back(e);
  /*
  ListCoord::iterator itlcb=itl->Realcoord.begin();
  ListCoord::iterator itlce=ite->Realcoord.end();
  bool found_placed=false;
  for (ListCoord::iterator itlc=itlcb; itlc!=itlce; ++itlc)
   if ((itlc->first==x) && (itlc->second==y))
    found_placed=true;

  if (found_placed)
   rep_locations++;
  */
 }
 else
 {
  /*
  ListCoord::iterator itlcb=itl->Realcoord.begin();
  ListCoord::iterator itlce=itl->Realcoord.end();

  bool found_identical=false;
  for (ListCoord::iterator itlc=itlcb; itlc!=itlce; ++itlc)
   if ((itlc->first==x) && (itlc->second==y))
    found_identical=true;

  if (!found_identical)
  {
  */
   itl->Realcoord.push_back(pair<float,float>(x,y));
   itl->DistNuc.push_back(dnuc);
   itl->vcells.push_back(cell);
   itl->v++;
  /*

  }
  else
   rep_locations_and_codes++;
  */
 }
}

void SpPseudoImage::GetCounts(size_t *sums,size_t &totsum,size_t *maxvalues,size_t &nv)
{
 int code;
 MyValPixelType v;
 for (size_t r=0; r<h; r++)
  for (size_t c=0; c<w; c++)
  {
   nv++;
   for (size_t pos=0; pos<ima[r][c].size(); pos++)
   {
    code=(ima[r][c])[pos].code;
    v=(ima[r][c])[pos].v;
    if (v>maxvalues[code])
     maxvalues[code]=v;
    totsum+=v;
    sums[code]+=v;
   }
  }
}
// @endcond
