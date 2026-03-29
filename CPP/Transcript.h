#ifndef _TRANSCRIPT_H
#define _TRANSCRIPT_H

#define ITK_LEGACY_FUTURE_REMOVE

#include <sstream>
#include <cmath>

#include "debugging.h"
#include "auxbin.h"

// Forward declaration of TranscriptSet, since it needs to be a friend class of Transcript
class TranscriptSet;

//! Transcript class.
/*!
 Transcript: a class to hold all information about a transcription event from a spatial transcriptomics description file in .csv or binary format
*/
class Transcript
{
 public:
  /*!
   Constructor that reads all relevant information from the .csv file and store it inside the class after doing some basic sanity checks
   \param  idp unsigned long long: the transcription event identifier
   \param  cell_idp string: unique identifier of the cell inside which the event transcription happens, or UNASSIGNED if happens outside any cell
   \param  overlap_nucp bool: true if the event happens inside the cell nucleus, false otherwise (even if outside cell)
   \param  feat_namep string: the identifier of the gene or clonotype that has generated this transcription event
   \param  xp float: the x-coordinate of the event, in micrometers
   \param  yp float: the y-coordinate of the event, in micrometers
   \param  zp float: the z-coordinate of the event, in micrometers
   \param  qvp int: the quality value that measures the confidence on correct assignment of this event to a gene or clonotype. Values in 0..40
   \param  fov_namep string: identifier of the region in histopathological image where the event has happened
   \param  nuc_distp float: distance between the event and the closest cell nucleus, or 0 if event happend inside cell nucleus
   \param  codewordp int: number asiigned to the gene or clonotype that has generated this trancription event, i.e., feat_namep.
   \param err &bool: Test variable that returns with true if all values pass the checks and transcript has been constructed and false if not.
  */
  Transcript(unsigned long long idp,std::string cell_idp,bool overlap_nucp,std::string feat_namep,float xp,float yp,float zp,int qvp,std::string fov_namep,float nuc_distp,int codewordp,bool &err);
  /*!
   Constructor that read all relevant information from a binary opened InputBuferedFile
   \param    f   InputBuferedFile Reference to a previously opened binary file
  */
  Transcript(InputBufferedFile &f);
  /*!
   @brief Function to return the unique event identifier
   @return long long The number assigned to the event
  */
  unsigned long long GetId() { return(id); };
  /*!
   @brief Function to return the unique cell identifier inside which this even happened, which can be UNASSIGNED if outside any cell
   @return  string The cell id or UNASSIGNED
  */
  std::string GetCellId() { return(cell_id); };
  /*!
   @brief Function to check if the event happened inside the cell nucleus or not
   @return bool true if event happened inside nucleus, false otherwise (including outside any cell)
  */
  bool OverlapNucleus() { return(overlap_nuc); };
  /*!
   @brief Function to return the name of the gene or clonotype that originates this event
   @return string Gene or clonotype identifier
  */
  std::string GetFeatureName() { return(feat_name); };
  /*!
   @brief Function to return the x coordinate of the event in micrometers
   @return  float x-value
  */
  float GetX() { return(x); };
  /*!
   @brief Function to return the y coordinate of the event in micrometers
   @return  float y-value
  */
  float GetY() { return(y); };
  /*!
   @brief Function to return the z coordinate of the event in micrometers
   @return  float z-value
  */
  float GetZ() { return(z); };
  /*!
   @brief Function to return the quality index assigned to the event
   @return  int Value in 0 (totally unconfident) to 40 (totally confident)
  */
  float GetQv() { return(qv); };
  /*!
   @brief Function to return the identifier of the area in which the event happened
   @return string An area indentifier
  */
  std::string GetFOV() { return(fov_name); };
  /*!
   @brief Function to return the distance from the event to the closest cell nucleus
   @return float Distance event-nucleus boundary, or 0 if event happens inside nucleus
  */
  float GetNucleusDist() { return(nuc_dist); };
  /*!
   @brief Function to return the codeword of the gene or clonotype that generated this event
   @return int The gene or clonotype codeword
  */
  int GetCodeword() { return(codeword); };
  /**
  @brief Function to return all the information about a transcription event as a text line with fields separated by the requested character
  @param csep char: The separation character
  @return string: The line with the information (without any final carriage return/new line character)
  */
  std::string GetAll(char csep);
  friend TranscriptSet;
 private:
  bool ValidRegionName( std::string rname );
  unsigned long long id;
  std::string cell_id;
  bool overlap_nuc;
  std::string feat_name;
  float x;
  float y;
  float z;
  float qv;
  std::string fov_name;
  float nuc_dist;
  int codeword;
  void SaveAsBinary(OutputBufferedFile &f);
};

#endif

