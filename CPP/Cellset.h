#ifndef _CELLSET_H
#define _CELLSET_H

#define ITK_LEGACY_FUTURE_REMOVE

#include "Cell.h"
#include "ProcessCSV.h"
#include "ndplaces.h"


//! Cellset class.
/*!
 Class CellSet: a class to hold all information about a set of cells read from the cell description file in .csv or binary format
*/
class CellSet
{
 public:
  const unsigned short default_chunksize_for_CellSet=128;   /*!< Size of memory chunk to save the cell set. 32 MB seems more than enough, even if we would have 900.000 cells... */
  const int UnknownCluster=-1;    /*!< Constant to indicate that a cell has not been attributed (at least, not yet) to any cluster */

  //! CellSet constructor
  /*!
   Constructor that takes the requested binary or .csv cell description file, reads from it all relevant information and stores it inside the class.
   File type (binary or .csv) is determined from its extension.
   \param cells_file_name string: The name (if accessible) or complete path of the .csv or binary file
   \param check bool: true to check that there are no any cells with duplicate identifiers. This check is slow, so the parameter is by default false.
  */
  CellSet(std::string cells_file_name,bool check=false);
  /**
  @brief Procedure to show all the data from the cell set in .csv format
  \param out &ostream Reference to an output stream to show the information, like cout, cerr or an opened ofstream object
  \param sep char: character to be used as .csv separator, by default a comma
  \param getcl bool: true if the attributed cluster should be shown, too. False if not. Default: false
   */
  void Show(std::ostream &out,char sep=',',bool getcl=false);
  /**
  @brief Function to return the number of cells in the cell set
  \return size_t: the number of cells currently in the cell set
   */
  size_t NumCells() { return(cset.size()); };
  /**
  @brief Function to return a reference to the cell from its cell identifier
  \param idc string the cell identifier
  \return &Cell: the reference to the cell
  */
  Cell& GetCellRef(std::string idc);
  /**
  @brief Function to return the ratio between cell nucleus area and cell area of a cell from its cell identifier.
  \param idc string the cell unique identifier
  \return float: the nucleus area/cell area ratio, in [0..1]
  */
  float GetAreaRatio(std::string idc) { return( id2place.contains(idc) ? cset[id2place[idc]].GetAreaRatio() : -1.0 ); };
  /**
  @brief Function to add the attibuted cluster to each cell reading it from an external .csv file
  \param clfile string with the name of the .csv clusters file
  */
  void AddClusterInfo(std::string clfile);
  /**
  @brief Function to load the cell or nucleus boundaries from an external .csv file
  \param bfile string with the name of the .csv cell or nucleus boundaries file
  \param cellboundary bool, true for setting the cell boundaries and false to set nucleus boundaries
  */
  void AddBoundaryInfo(std::string bfile,bool cellboundary);
  /**
  @brief Function to return the attributed cluster of a cell from its cell identifier.
  \param idc string the cell unique identifier
  \return int: the number of the cluster attributed to this cell
  */
  int GetAttClust(std::string idc) { return( id2place.contains(idc) ? cset[id2place[idc]].GetCluster() : UnknownCluster ); };

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
  std::vector<Cell> cset;                              /*!< A vector of cells to hold all of them */
  void CellSetFromCSV(std::string fname);
  void CellSetFromBin(std::string fname);
  void CheckCellSet();
  std::unordered_map<std::string,size_t>  id2place; /*!< A map to recover fast a cell from its identifier */
  void SetCellCluster(std::string idc,int clc,bool &err); /*!< Sets the cluster of a particular cell with given cell identifier */
  void SetBoundary(std::string idc,std::vector<float> const &bx,std::vector<float> const &by,std::vector<size_t> const &startp,bool cellboundary,bool &err);
  void FillId2place();
  size_t num_cells_with_boundary;
  size_t num_nuc_with_boundary;
};

#endif
