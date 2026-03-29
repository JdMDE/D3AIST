#ifndef _TRANSCRIPTSET_H
#define _TRANSCRIPTSET_H

#include "Transcript.h"
#include "ProcessCSV.h"
#include "ndplaces.h"

//! TranscriptSet class.
/*!
 Class TranscriptSet: a class to hold all information about a set of transcription events read from the event description file in .csv or binary format
*/
class TranscriptSet
{
 public:
  const unsigned short default_chunksize_for_TranscriptSet=4096; /*!< Size of memory chunck to save the transcript set. Even 4 GB is not enough, but it depends on the memory of your machine */
  const std::string UnassignedCellCode="UNASSIGNED";   /*!< This is the string used in the csv files to indicate that the transcript is outside any cell. */
  /*!
   Constructor that takes the requested binary or .csv cell description file, reads from it all relevant information and stores it inside the class.
   File type (binary or .csv) is determined from its extension.
   \param tr_file_name string: The name (if accessible) or complete path of the .csv or binary file
   \param check bool: true to check that there are no any transcritps with duplicate identifiers. This check is slow, so the parameter is by default false.
  */
  TranscriptSet(std::string tr_file_name,bool check=false);
  /**
   @brief Procedure to show all the data from the transcript set in .csv format
   \param out &ostream Reference to an output stream to show the information, like cout, cerr or an opened ofstream object
   \param sep char: character to be used as .csv separator, by default a comma
  */
  void Show(std::ostream &out,char sep=',');
  /**
  @brief Function to return the number of transcription events in the transcript set
  \return size_t: the number of transcription events currently in the transcript set
  */
  size_t GetNumTr() { return(tset.size()); };
  /*!
   @brief Function to return the x coordinate of event with given identifier, in micrometers
   \param t size_t The event identifier
   @return  float x-value
  */
  float GetX(size_t t) { return(tset[t].GetX()); };
  /*!
   @brief Function to return the y coordinate of event with given identifier, in micrometers
   \param t size_t The event identifier
   @return  float y-value
  */
  float GetY(size_t t) { return(tset[t].GetY()); };
  /*!
   @brief Function to return the identifier of the area in which an event with given identifier happened
   \param t size_t The event identifier
   @return string An area indentifier
  */
  std::string GetFOV(size_t t) { return(tset[t].GetFOV()); };
  /*!
   @brief Function to return the identifier of the cell inside which an event with given identifier happened, or UNASSIGNED if outside any cell
   \param t size_t The event identifier
   @return string A cell identifier
  */
  std::string GetCellId(size_t t) { return(tset[t].GetCellId()); };
  /*!
   @brief Function to check if an event with given identifier happened inside any cell or not
   \param t size_t The event identifier
   @return bool True if the event happened inside any cell (including inside its nucleus) or false if not
  */
  bool InsideCell(size_t t) { return(tset[t].GetCellId()!=UnassignedCellCode); };
  /*!
   @brief Function to check if an event with given identifier happened inside any cell nucleus
   \param t size_t The event identifier
   @return bool True if the event happened inside any cell nucleus or false if not (even if outside any cell)
  */
  bool OverlapNucleus(size_t t) { return(tset[t].OverlapNucleus()); };
  /*!
   @brief Function to return the quality index of the vent with given identifier
   \param t size_t The event identifier
   @return int The quality index, an integer in 0..40
  */
  int GetQv(size_t t) { return(tset[t].GetQv()); };
  /*!
   @brief Function to return the codeword of the gene or clonotype the originated the event with given identifier
   \param t size_t The event identifier
   @return int The codeword of the gene or clonotype
  */
  int GetCW(size_t t) { return(tset[t].GetCodeword()); };
  /*!
   @brief Function to return the distance to the closest nucleus of the event with given identifier
   \param t size_t The event identifier
   @return float The distance to the closest nucleus in micrometers
  */
  float GetDNuc(size_t t) { return(tset[t].GetNucleusDist()); };
  /**
  @brief Procedure to save the cell set in binary format in a file debugging the writing bytes process
  \param fname string The name of the file to be created
  */
  void SaveAsBinaryWithDebug(std::string fname);
   /**
  @brief Procedure to save the cell set in binary format in a file
  \param fname string The name of the file to be created
  */
  void SaveAsBinary(std::string fname);
 private:
  std::vector<Transcript> tset;                     /*!< A vector of transcription events to hold all of them */
  bool ValidRegName( std::string rname ) { return(tset[0].ValidRegionName(rname)); };
  void TranscriptSetFromCSV(std::string fname);
  void TranscriptSetFromBin(std::string fname);
  bool CheckTranscriptSet();
};

#endif

