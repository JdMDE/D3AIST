#ifndef _COLORS_H
#define _COLORS_H

/**
 @brief Structure to account for a color
*/
typedef struct
{
 unsigned char r; /*!< The red component, in [0.255] */
 unsigned char g; /*!< The green component, in [0.255] */
 unsigned char b; /*!< The blue component, in [0.255] */
} RGB_color;

const RGB_color black={0,0,0};
const RGB_color white={255,255,255};

/**
 @brief Structure to account for each entry in a color table
*/
typedef struct
{
 int line_in_file;         /*!< The line in the color assignment .csv file where this color appears */
 int gene_codeword;        /*!< The codeword of the gene to which the color is to be assigned */
 std::string gene_name;   /*!< The short name of the gene to which the color is to be assigned */
 RGB_color col;             /*!< The color to assign, as a RGB_color structure */
} ColortabEntry;

#endif
