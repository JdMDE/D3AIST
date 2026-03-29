// @cond DO_NOT_DOCUMENT
#ifndef _SPPSEUDOIMAGE_H
#define _SPPSEUDOIMAGE_H

#define ITK_LEGACY_FUTURE_REMOVE

#include <vector>
#include "Cell.h"

class SpPseudoImage
{
 public:
  typedef float MyPoint[2];
  typedef float My3DPoint[3];
  typedef unsigned int My3DIndex[3];
  typedef float Spacing[3];
  typedef unsigned short MyValPixelType;
   typedef std::vector<std::pair<float,float>> ListCoord;
   typedef struct
   {
    int code;
    ListCoord Realcoord;
    std::vector<float> DistNuc;
    std::vector<std::reference_wrapper<Cell>> vcells;
    MyValPixelType v;
   } entry;
   typedef std::vector<entry> PixContent;
   typedef PixContent *PixContentPointer;
   SpPseudoImage(unsigned int width,unsigned int height,Spacing rsp,My3DPoint origin,int cr);
   ~SpPseudoImage();
   unsigned int GetWidth() { return(w); };
   unsigned int GetHeight() { return(h); };
   PixContent GetPixList(unsigned int r,unsigned int c) { return(ima[r][c]); };
   int GetCoderest() { return(coderest); };
   bool PPhysicalPointToIndex(My3DPoint &p,My3DIndex &idx);
   bool PIndexToPhysicalPoint(unsigned int r,unsigned int c,MyPoint &p);
   void SimpleInsertOrIncrementPixelValue(My3DIndex &idx);
   void InsertOrIncrementPixelValue(My3DIndex &idx,float dnuc,Cell &cell,size_t &rep_locations,size_t &rep_locations_and_codes,float x=0.0,float y=0.0);
   void GetCounts(size_t *sums,size_t &totsum,size_t *maxvalues,size_t &nv);
   float GetXOrig() { return(orig[0]); };
   float GetYOrig() { return(orig[1]); };
   float GetXCorner() { return(corner[0]); };
   float GetYCorner() { return(corner[1]); };
   float GetXSpacing() { return(sp[0]); };
   float GetYSpacing() { return(sp[1]); };
   void TransfTest();
 private:
   unsigned int w,h;
   Spacing sp;
   My3DPoint orig;
   MyPoint corner;
   PixContent **ima;
   int coderest;
};

#endif
// @endcond

