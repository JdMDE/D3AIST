#ifndef _BUILDIM_H
#define _BUILDIM_H

#define ITK_LEGACY_FUTURE_REMOVE

#include <array>

#include "debugging.h"
#include "SpRData.h"

// @defgroup ParametersForDescriptionFile
/**
 @addtogroup ParametersForDescriptionFile
 @{
*/
/**
 @brief ParamTypes documentation
         Enumeration with the types of values we can expect for our parameters. Int2 means 'integer integer',
         Float2 means 'float float', Float4 means 'float float float float'. All of these are separated by one or more blank spaces.
         ListInt means 'int,int,....,int', i.e., a list of comma-separated integers.
         Similarly for ListString. In these cases, spaces are not allowed and there must be no comma at the end.
*/
enum class ParamTypes { Boolean,String,Int,Int2,Float2,Float4,ListInt,ListString,NoType };

/**
 @brief Structure to parse each parameter of the configuration file
*/
typedef struct
{
 std::string ParamName;          /*!< Name of the parameter, as given in the configuration file */
 std::string AltParamName;       /*!< Alternative parameter name when two forms of the parameter are mutually exclusive */
 ParamTypes   ParType;             /*!< Type of parameter. One of the values of the enum ParamTypes */
 ParamTypes   AltParType;          /*!< Type of alternative parameter when two forms of the parameter are mutually exclusive */
} Parameter;

const unsigned int CurrentNumberParameters=20;     /*!< This constant must be changed if parameters of the configuration file are added or removed. */

const std::string dirslash="/";         /*!< The directory slash in Unix-like systems. */

const Parameter Params[CurrentNumberParameters]=             /*!< Array of Parameter structures with all the currently allowed parameters */
{
 {"PrependPath",    "",            ParamTypes::String, ParamTypes::NoType},       // 0
 {"Genes",          "",            ParamTypes::String, ParamTypes::NoType},       // 1
 {"Cells",          "",            ParamTypes::String, ParamTypes::NoType},       // 2
 {"Transcripts",    "",            ParamTypes::String, ParamTypes::NoType},       // 3
 {"Clusters",       "",            ParamTypes::String, ParamTypes::NoType},       // 4
 {"CellBoundaries", "",            ParamTypes::String, ParamTypes::NoType},       // 5
 {"NucBoundaries",  "",            ParamTypes::String, ParamTypes::NoType},       // 6
 {"Spacing",        "",            ParamTypes::Float2, ParamTypes::NoType},       // 7
 {"Dimensions",     "",            ParamTypes::Int2,   ParamTypes::NoType},       // 8
 {"PixPerBin",      "",            ParamTypes::Int,    ParamTypes::NoType},       // 9
 {"AoIWin",         "AoIID",       ParamTypes::Float4, ParamTypes::String},       // 10
 {"QvThres",        "",            ParamTypes::Int,    ParamTypes::NoType},       // 11
 {"OnlyTargets",    "",            ParamTypes::Boolean,ParamTypes::NoType},       // 12
 {"TargetIds",      "TargetNames", ParamTypes::ListInt,ParamTypes::ListString},   // 13
 {"GenMarks",       "",            ParamTypes::Boolean,ParamTypes::NoType},       // 14
 {"OutRfile",       "",            ParamTypes::String, ParamTypes::NoType},       // 15
 {"OutImfile",      "",            ParamTypes::String, ParamTypes::NoType},       // 16
 {"ColorFile",      "",            ParamTypes::String, ParamTypes::NoType},       // 17
 {"OutHistfile",    "",            ParamTypes::String, ParamTypes::NoType},       // 18
 {"CheckThisFile",  "",            ParamTypes::Boolean,ParamTypes::NoType}        // 19
};
/**
 @}
*/

#endif

