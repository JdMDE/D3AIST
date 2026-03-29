#ifndef _GENE_H
#define _GENE_H

#include <sstream>

#include "debugging.h"
#include "auxbin.h"
#include "colors.h"

// @defgroup ConstantsForGenes
/**
 @addtogroup ConstantsForGenes
 @{
*/
const int max_codeword=500;                                               /*!< The maximum number of codewords allowed */

const int NumCat=2;                                                        /*!< The number of entity (gene/clonotype) categories, currently only two */
const std::string catnames[NumCat]={"current","base"};                   /*!< The names of entity categories, current or base */
const int NumDesc=2;                                                       /*!< The number of possible classes an entity can belong to, currently only two */
const std::string descnames[NumDesc]={"gene","negative_control"};         /*!< The names of the classes for an entity: gene or negative_control */
/**
 @}
*/
const std::string NoName="Name_UNKNOWN";
const std::string NoId="Id_UNKNOWN";
const int NoCode=-1;

// Forward declaration, since we need GenePanel to be a friend class of Gene
class GenePanel;

//! Gene class.
/*!
 Gene: a class to hold all information about a gene, clonotype or control from a spatial transcriptomics description file in a .csv or binary format gene panel
*/
class Gene
{
 public:
  /*!
   Constructor that reads all relevant information from the .csv file and store it inside the class after doing some basic sanity checks
   \param  cw int: the codeword assigned to the gene or clonotype
   \param  gc int: the gene coverage
   \param  cat string: the gene category (current, if selected for the analysis, or base otherwise)
   \param  gid string: the ENSEMS gene identifier
   \param  gname string: the gene short name
   \param  gd string: the descriptor of the entity (gene or negative_control)
   \param err &bool: Test variable that returns with true if all values pass the checks and transcript has been constructed and false if not.
  */
  Gene(int cw,int gc,std::string cat,std::string gid,std::string gname,std::string gd,bool &err);
  /*!
   Constructor that read all relevant information from a binary opened InputBuferedFile
   \param    f   InputBuferedFile Reference to a previously opened binary file
  */
  Gene(InputBufferedFile &f);
  /*!
   @brief Function to return the gene codeword
   @return int The codeword
  */
  int GetCodeword() { return(codeword); };
  /*!
   @brief Function to return the gene coverage
   @return int The coverage
  */
  int GetCoverage() { return(gene_coverage); };
  /*!
   @brief Function to return the gene category (current or base)
   @return string The category
  */
  std::string GetCat() { return(catnames[gene_cat]); };
  /*!
   @brief Function to return the gene ENSEMS identifier
   @return string The category
  */
  std::string GetId() { return(id); };
  /*!
   @brief Function to return the gene short name
   @return string The name
  */
  std::string GetName() { return(name); };
  /*!
   @brief Function to return the gene descriptor (gene or negative_control)
   @return string The descriptor
  */
  std::string GetDesc() { return(descnames[gene_desc]); };
  /*!
   @brief Function to check if the entity is a gene or not
   @return bool True for a gene, false otherwise
  */
  bool IsGene() { return(descnames[gene_desc]=="gene"); };
  /**
  @brief Function to return all the information about a gene as a text line with fields separated by the requested character
  @param sep char: The separation character
  @return string: The line with the information (without any final carriage return/new line character)
  */
  std::string GetAll(char sep=',');
  friend GenePanel;
 private:
  int GetCatNum() { return(gene_cat); };
  int GetDescNum() { return(gene_desc); };
  RGB_color GetColor() { return(color); };
  void SetColor(const RGB_color &c) { color=c; };
  int codeword;                   /*!< Gene codeword */
  int gene_coverage;             /*!< Gene coverage */
  int gene_cat;                   /*!< Gene category */
  std::string id;                /*!< Gene ENSEMS identfier */
  std::string name;              /*!< Gene short name */
  int gene_desc;                  /*!< Gene descriptor */
  RGB_color color;                 /*!< Gene color. Color will not be stored in the binary or .csv files. It is not intrinsic to the gene, it is something assigned when used. */
  void SaveAsBinary(OutputBufferedFile &f);
};

#endif
