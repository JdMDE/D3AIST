#ifndef _CELL_H
#define _CELL_H

#define ITK_LEGACY_FUTURE_REMOVE

#include <sstream>
#include <cmath>
#include <unordered_map>

#include "debugging.h"
#include "auxbin.h"

// Forward declaration, since we need CellSet to be a friend class of Cell
class CellSet;

//! Cell class.
/*!
 Class Cell: a class to hold all information about a cell read from a spatial transcriptomics cell description file in .csv or binary format
*/
class Cell
{
 public:
  /*!
   Constructor that reads all relevant information from the .csv file and store it inside the class after doing some basic sanity checks
   \param   id string: the cell identifier
   \param    x float: x-position of cell centroid in micrometers
   \param    y float: y-position of cell centroid in micrometers
   \param  trc int: Number of transcript counts inside cell
   \param  cpc int: Number of control probe counts inside cell
   \param  ccc int: Number of control codeword counts inside cell
   \param  ucc int: Number of unassigned codeword counts inside cell
   \param  dcc int: Number of deprecated codeword counts inside cell
   \param   tc int: Total number of counts inside cell. It must be trc+cpc+ucc+dcc
   \param area float: Cell area in sq. micrometers
   \param nuc_area float: Nucleus cell area in sq. micrometers
   \param err &bool: Test variable that returns with true if all values pass the checks and cell has been constructed and false if not.
  */
  Cell(std::string id,float x,float y,int trc,int cpc,int ccc,int ucc,int dcc,int tc,float area,float nuc_area,bool &err);
  /*!
   Constructor that read all relevant information from a binary opened InputBuferedFile
   \param    f   InputBuferedFile Reference to a previously opened binary file
   \param    nbc &size_t Variable to count the number of cell boundaries read from the binary file. Passed as reference to be updated
   \param    nnc &size_t Variable to count the number of nucleus boundaries read from the binary file. Passed as reference to be updated
  */
  Cell(InputBufferedFile &f,size_t &nbc,size_t &nnc);
  /*!
   @brief Function to return the unique cell identifier
   @return  string The cell id
  */
  std::string GetId() { return(id); };
  /*!
   @brief Function to return the x coordinate of the cell centroid in micrometers
   @return  float x-value
  */
  float GetX() { return(x); };
  /*!
   @brief Function to return the y coordinate of the cell centroid in micrometers
   @return  float y-value
  */
  float GetY() { return(y); };
  /*!
   @brief Function to return the number of transcription events happened inside the cell
   @return int Number of transcription events
  */
  int GetTranscriptC() { return(trc); };
  /*!
   @brief Function to return the number of control probe events happened inside the cell
   @return int Number of control probe events
  */
  int GetControlProbeC() { return(cpc); };
  /*!
   @brief Function to return the number of events corresponding to control codewords happened inside the cell
   @return int Number of events corresponding to codewords
  */
  int GetControlCodewordC() { return(ccc); };
  /*!
   @brief Function to return the number of events happened inside the cell not assigned to any known codeword
   @return int Number of unassigned events
  */
  int GetUnassignedCodewordC() { return(ucc); };
  /*!
   @brief Function to return the number of events happened inside the cell assigned to deprecated codeword
   @return int Number of deprecated events
  */
  int GetDeprecatedCodewordC() { return(dcc); };
  /*!
   @brief Function to return the total number of events happened inside the cell
   @return int Total number of events
  */
  int GetTotalC() { return(tc); };
  /*!
   @brief Function to return the cell area in sq. micrometers
   @return float Area of the cell
  */
  float GetArea() { return(area); };
  /*!
   @brief Function to return the cell nucleus area in sq. micrometers
   @return float Area of the cell nucleus
  */
  float GetNucArea() { return(nuc_area); };
  /*!
   @brief Function to return the number of the cluster the cell is assigned to
   @return int Number of the attributed cluster
  */
  int GetCluster() { return(cluster); };
  /*!
   @brief Procedure to set the number of the cluster the cell is assigned to
  */
  void SetCluster(int clc) { cluster=clc; };
  /**
  @brief Function to return the ratio between cell nucleus area and cell area, a real number in ]0..1[. Used as suggestive indicator of tumoral character.
  @return float: Quotient nucleus cell area / cell area
  */
  float GetAreaRatio() { return(nuc_area/area); };
  /**
  @brief Function to return all the information about a cell as a text line with fields separated by the requested character
  @param csep char: The separation character
  @param getclus bool: true to indicate that the attributed cluster number should be written, too. False if not.
  @return string: The line with the information (without any final carriage return/new line character)
  */
  std::string GetAll(char csep,bool getclus);
  friend CellSet;
 private:
  std::string id;                 /*!< Cell identifier */
  float x;                          /*!< x-position of cell centroid */
  float y;                          /*!< y-position of cell centroid */
  int trc;                          /*!< Transcript counts inside cell */
  int cpc;                          /*!< Control probe counts inside cell */
  int ccc;                          /*!< Control codeword counts inside cell */
  int ucc;                          /*!< Unassigned codeword counts inside cell */
  int dcc;                          /*!< Deprecated codeword counts inside cell */
  int tc;                           /*!< Total counts inside cell */
  float area;                       /*!< Cell area in sq. micrometers */
  float nuc_area;                   /*!< Cell nucleus area in sq. micrometers */
  int cluster;                      /*!< The number of the cluster attributed to this cell */
  std::vector<float> cboundx;
  std::vector<float> cboundy;     /*!< The boundary of the cell, as a closed polygon (vector of 2D points) */
  std::vector<float> nuccboundx;
  std::vector<float> nuccboundy;  /*!< The boundary of the cell nucleus, as a succession of closed polygons */
  std::vector<size_t> startpoints;  /*!< The indices of the nuccboundx,nuccboundy arrays where a new chain starts. It has always at least a first element with value 0. */
  void SaveAsBinary(OutputBufferedFile &f);
};



#endif
