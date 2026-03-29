#ifndef _GENEPANEL_H
#define _GENEPANEL_H

#include <limits>
#include <algorithm>
#include <map>

#include "Gene.h"

//! GenePanel class.
/*!
 Class GenePanel: a class to hold all information about a gene panel read from the description file in .csv or binary format
*/
class GenePanel
{
 public:
  const unsigned short default_chunksize_for_GenePanel=16; /*!< Size of memory chunk to save the gene panel. 16 MB seems more than enough, even if we would have 20.0000 genes... */
  const int Target=0;   /*!< This is because genes in category current are the targets, and we have stored current at place 0 */
  //! GenePanel constructor
  /*!
   Constructor that takes the requested binary or .csv gene panel file, reads from it all relevant information and stores it inside the class.
   File type (binary or .csv) is determined from its extension.
   \param panel_file_name string: The name (if accessible) or complete path of the .csv or binary file
   \param check bool: true to check that there are no any genes with duplicate identifiers. With many genes this check could be slow, so the parameter is by default false.
  */
  GenePanel(std::string panel_file_name,bool check=false);
  /**
   @brief Procedure to show all the data from the gene panel in .csv format
   \param out &ostream Reference to an output stream to show the information, like cout, cerr or an opened ofstream object
   \param sep char: character to be used as .csv separator, by default a comma
  */
  void Show(std::ostream &out,char sep=',');

  /**
   @brief Function to get the minimum value of the codewords of all genes/clonotypes/control probes
   @return int The smallest of the codewords
  */
  int GetMinCodeword() { return(mincod); };
  /**
   @brief Function to get the maximum value of the codewords of all genes/clonotypes/control probes
   @return int The bigest of the codewords
  */
  int GetMaxCodeword() { return(maxcod); };
  /**
   @brief Function to get the minimum value of the coverage of all genes/clonotypes/control probes
   @return int The smallest of the coverages
  */
  int GetMinCoverage() { return(mincov); };
  /**
   @brief Function to get the maximum value of the coverage of all genes/clonotypes/control probes
   @return int The bigest of the coverages
  */
  int GetMaxCoverage() { return(maxcov); };

  /**
   @brief Function to get a vector with the codewords of all the genes considered as targets, sorted from lower to higher
   @return vector<int> The sorted codewords
  */
  std::vector<int> GetSortedTargets();
  /**
   @brief Function to get a vector with the codewords of all the genes, either targets or not, sorted from lower to higher
   @return vector<int> The sorted codewords
  */
  std::vector<int> GetSortedGenes();

  /**
   @brief Function to return the gene ENSEMS identifier from its codeword
   \param cw int The gene codeword
   @return string The gene identifier
  */
  std::string GetIdFromCodeword(int cw);
  /**
   @brief Function to return the gene short name from its codeword
   \param cw int The gene codeword
   @return string The gene name
  */
  std::string GetNameFromCodeword(int cw);
  /**
   @brief Function to return the gene codeword from its short name
   \param genename string The gene name
   @return int The gene codeword
  */
  int GetCodewordFromName(std::string genename);

  /**
   @brief Function to check if the entity with a given codeword is a gene or not
   \param cw int The entity codeword
   @return bool True if it is a gene, false otherwise
  */
  bool IsGene(int cw) { return( isgene.find(cw)!=isgene.end() ? isgene[cw] : false); };

  /**
   @brief Procedure to fill the internal table of colors associated to each gene (for representation purposes) from a vector of colors
   \param v vector<ColorTabEntry> A vector of entries to the color table
  */
  void FillColors(std::vector<ColortabEntry> &v);
  /**
   @brief A function to return the color associated to the gene with a given codeword, as stored in the internal table
   \param cw int The codeword
   @return A RGB_color structure with the color
  */
  RGB_color GetColorFromCodeword(int cw);

  /**
  @brief Procedure to save the gene panel in binary format in a file debugging the writing bytes process
  \param fname string The name of the file to be created
  */
  void SaveAsBinaryWithDebug(std::string fname);

  /**
  @brief Procedure to save the cell set in binary format in a file
  \param fname string The name of the file to be created
  */
  void SaveAsBinary(std::string fname);

 private:
  std::vector<Gene> panel;       /*!< A vector of genes to hold all of them */
  std::map<int,bool> isgene;     /*!< A map that relates the gene/entity codeword with a boolean value: true for a gene, false otherwise */
  int mincod;                     /*!< The minimum codeword of the genes in the list */
  int maxcod;                     /*!< The maximum codeword of the genes in the list */
  int mincov;                     /*!< The minimum coverage of the genes in the list */
  int maxcov;                     /*!< The maximum coverage of the genes in the list */

  void CheckPanel();
  void GenePanelFromCSV(std::string fname);
  void GenePanelFromBin(std::string fname);
  void FillIsGeneMap();
  void SetColorOfGene(int cw,RGB_color &c);
};
#endif
