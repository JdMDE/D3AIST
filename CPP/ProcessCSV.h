#ifndef _PROCESSCSV_H
#define _PROCESSCSV_H

#include <string_view>
#include <vector>
#include <iostream>
#include <fstream>
#include <string_view>
std::vector<std::string_view> parseCSVRow(std::string_view row,std::string sep);
std::vector<std::vector<std::string_view>> readCSV(const std::string& filename);

#endif

