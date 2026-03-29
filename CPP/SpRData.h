#ifndef _SPRDATA_H
#define _SPRDATA_H

#define ITK_LEGACY_FUTURE_REMOVE

#include <functional>

#include "Gene.h"
#include "Genepanel.h"
#include "Cell.h"
#include "Cellset.h"
#include "Transcript.h"
#include "Transcriptset.h"

#include "imlimits.h"
#include "imparams.h"
#include "imtypes.h"

#include "SpPseudoImage.h"

//! SpRData class.
/*!
 SpRData: a class to generate the R files and/or output images with transcription event locations
*/
class SpRData
{
 public:
  /*!
   Constructor that integrates all the input data inside the SpRdata information
   \param ipar ImParams& reference to an ImParam structure with all the needed characteristics and requests to build the output
   \param panelref GenePanel& reference to a GenePanel with the genes that can be used (all of them, even those not actually used)
   \param csref CellSet& reference to a CellSet with information on the cells present in the input image
   \param tsref TranscriptSet& reference to a TranscriptSet with all the transcription events appearing in the input image (all of them, even those outside the area of interest)
  */
  SpRData(ImParams &ipar,GenePanel &panelref,CellSet &csref,TranscriptSet &tsref);
  /*!
    Default destructor. Calls internally to the destructor of the subordinated class, SpPseudoImage
  */
  ~SpRData() { if (psima!=nullptr) psima->~SpPseudoImage(); };
  /*!
   Procedure to build the R files and/or images, as requested, to be called after the constructor.
   \param histname string Name of the histogram file (csv file with number of events of each gene inside the interest area) or none (default value) for not generating such file
  */
  void Build(std::string histname="none");
  /*!
   Procedure to save the R files and/or images generated.
   \param outfname string Base name of the generated files (see direction to fill the configuration file, parameters OutRName and OutImName)
   \param sep char Separation character to be used in the tables that R will read. Default: a comma
  */
  void Save(std::string outfname,char sep=',');
  /*!
   Procedure to save the transcription event locations as points in binary images
   \param outfname string Base name of the generated files (see direction to fill the configuration file, parameters OutRName and OutImName)
   \param idescname string Name of the parameter file from which the task was described, without extension. Used to build the .fuse file to be used by imfuse and the groovy script to extract the requested area of the original image.
  */
  void SaveAsImages(std::string outfname,std::string idescname);
 private:
  std::string SkipComments(std::ifstream &f);
  bool IsInside(float x,float y) { return ((ulx<=x) && (x<=lrx) && (uly<=y) && (y<=lry)); };
  bool FindOrder(int ecode,size_t &order);
  //void FillBackgroundImage(size_t w,size_t h,std::string bgn);
  void ScriptBackgroundImage(std::string scr,std::string bgn);
  GenePanel &panel;
  CellSet &cs;
  TranscriptSet &ts;
  bool pseudoimage_built;
  SpPseudoImage::Spacing orig_sp,reduced_sp;
  int pix_per_bin;
  unsigned long origwidth,origheight;
  float ulx,uly,lrx,lry;
  std::string fov;
  bool extract_by_fov;
  bool only_targets;
  int qvth;
  std::vector<int> ecodes;
  bool genmarks;
  SpPseudoImage *psima;
};

#endif

