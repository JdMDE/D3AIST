#include "ProcessCSV.h"

using namespace std;

vector<string_view> parseCSVRow(string_view row,string sep)
{
    vector<string_view> fields;
    size_t start = 0;

    size_t i = 0;
    while (i < row.length())
    {
        if (sep.find(row[i])!=string::npos)
        {
            while ((sep.find(row[i])!=string::npos) && (i<row.length()))
             i++;
            i--;
            if (row[start]=='"')
                fields.push_back(row.substr(start+1, i - start-2));
            else
                fields.push_back(row.substr(start, i - start));
            start = i + 1;
        } 
        i++;
    }
    fields.push_back(row.substr(start));
    
    return fields;
}

vector<vector<string_view>> readCSV(const string& filename,string sep)
{
    vector<vector<string_view>> data;
    ifstream file(filename);
    
    if (!file.is_open())
    {
        cerr << "Failed to open file: " << filename << endl;
        return data;
    }

    string line;
    while (getline(file, line))
      data.push_back(parseCSVRow(line,sep));
   

    file.close();
    return data;
}
