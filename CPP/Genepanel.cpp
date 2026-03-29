#include "Genepanel.h"

using namespace std;

extern string progname;
extern bool DEB;

GenePanel::GenePanel(string fname,bool check)
{
 if (fname.find(".csv")!=string::npos || fname.find(".CSV")!=string::npos)
 {
  if (DEB)
   cout << "Reading gene panel from .csv file " << fname << ".\n";
  ifstream f(fname.c_str());
  if (!f.is_open())
   ERRPROG("Cannot open panel file " << fname << " to read. Check for existence and permissions.\n";)
  f.close();
  GenePanelFromCSV(fname);
 }
 else
 {
  if (DEB)
   cout << "Reading gene panel from binary file " << fname << ".\n";
  GenePanelFromBin(fname);
 }

 FillIsGeneMap();

 if (check)
  CheckPanel();
}

void GenePanel::GenePanelFromCSV(std::string fname)
{
 mincod=mincov=numeric_limits<int>::max();
 maxcod=maxcov=numeric_limits<int>::min();

 string gene_entry;
 unsigned long line=0;
 int codeword,genecov;
 string cat,idg,name,desc;
 bool err;

 ifstream f(fname.c_str());
 // First, read the header line
 getline(f,gene_entry);
 if (DEB)
 {
  cout << "Header from file " << fname << ":\n";
  cout << gene_entry << endl;
 }
 line++;

 // Then, start reading content lines
 getline(f,gene_entry);
 line++;
 while (!f.eof())
 {
  stringstream st(gene_entry);

  codeword=genecov=0;
  cat=idg=name=desc="";

  st >> codeword >> genecov >> cat >> idg >> name >> desc;

  // Some lines (negative controls) have no name
  if (desc=="")
  {
   desc=name;
   name=NoName;
  }

  err=false;
  Gene g(codeword,genecov,cat,idg,name,desc,err);
  if (err)
   ERRPROG("Error at line " << line << " of input file " << fname << ". Line is:\n" << gene_entry << endl;)

  if (codeword<mincod)
   mincod=codeword;
  if (codeword>maxcod)
   maxcod=codeword;
  if (genecov<mincov)
   mincov=genecov;
  if (genecov>maxcov)
   maxcov=genecov;

  panel.push_back(g);

  getline(f,gene_entry);
  if (!f.eof())
   line++;
 }
 f.close();

 if (DEB)
  cout << panel.size() << " items read.\n";
}

void GenePanel::GenePanelFromBin(std::string fname)
{
 InputBufferedFile f(fname);

 f.readb(mincod);
 f.readb(maxcod);
 f.readb(mincov);
 f.readb(maxcov);
 size_t ps=0;
 f.readb(ps);
 for (size_t i=0; i<ps; i++)
 {
  Gene g(f);
  panel.push_back(g);
 }

 if (DEB)
  cout << panel.size() << " items read.\n";
}

void GenePanel::FillIsGeneMap()
{
 for (vector<Gene>::iterator it=panel.begin(); it!=panel.end(); ++it)
  isgene[it->GetCodeword()]=it->IsGene();
}

void GenePanel::CheckPanel()
{
 if (DEB)
 {
  cout << "Checking " << panel.size() << " items...";
  cout.flush();
 }

 for (size_t i=0; i<panel.size()-1; i++)
  for (size_t j=i+1; j<panel.size(); j++)
  {
   if (panel[i].GetCodeword()==panel[j].GetCodeword())
    ERRPROG("Duplicated codeword (" << panel[i].GetCodeword() << ") between items " << i+1 << " and " << j+1 << endl;)
   if (panel[i].GetId()==panel[j].GetId())
    ERRPROG("Duplicated gene identifier (" << panel[i].GetId() << ") between items " << i+1 << " and " << j+1 << endl;)
   if (panel[i].GetName()!=NoName && panel[i].GetName()==panel[j].GetName())
    ERRPROG("Duplicated name identifier (" << panel[i].GetName() << ") between items " << i+1 << " and " << j+1 << endl;)
  }

 if (DEB)
  cout << "Done. No duplicated items found.\n";
}

vector<int> GenePanel::GetSortedTargets()
{
 vector<int> ret;
 for (size_t g=0; g<panel.size(); g++)
  if (panel[g].GetCatNum()==Target)
   ret.push_back(panel[g].GetCodeword());
 sort(ret.begin(),ret.end());
 return(ret);
}

vector<int> GenePanel::GetSortedGenes()
{
 vector<int> ret;
 for (size_t g=0; g<panel.size(); g++)
  if (panel[g].IsGene())
   ret.push_back(panel[g].GetCodeword());
 sort(ret.begin(),ret.end());
 return(ret);
}

void GenePanel::Show(ostream &out,char sep)
{
 for (size_t g=0; g<panel.size(); g++)
  out << panel[g].GetAll(sep) << endl;
 out << "Gene codewords in [" << mincod << "," << maxcod << "]\n";
 out << "Gene coverages in [" << mincov << "," << maxcov << "]\n";
}

string GenePanel::GetIdFromCodeword(int cw)
{
 size_t i=0;
 while (i<panel.size() && panel[i].GetCodeword()!=cw)
  i++;
 if (i==panel.size())
 {
  WARNPROG("GetIdFromCodeword: Target with codeword " << cw << " not found in gene panel. Returning " << NoId << " as Id\n";)
  return(NoId);
 }
 return(panel[i].GetId());
}

RGB_color GenePanel::GetColorFromCodeword(int cw)
{
 size_t i=0;
 while (i<panel.size() && panel[i].GetCodeword()!=cw)
  i++;
 if (i==panel.size())
 {
  WARNPROG("GetColorFromCodeword: Target with codeword " << cw << " not found in gene panel. Returning " << NoId << " as Id\n";)
  return(black);
 }
 RGB_color col=panel[i].GetColor();
 if ((col.r==black.r) && (col.g==black.g) && (col.b==black.b))
 {
  WARNPROG("GetColorFromCodeword: Target with codeword " << cw << " not found in color table. Returning white.\n";)
  return(white);
 }
 return(panel[i].GetColor());
}


string GenePanel::GetNameFromCodeword(int cw)
{
 size_t i=0;
 while (i<panel.size() && panel[i].GetCodeword()!=cw)
  i++;
 if (i==panel.size())
 {
  cerr << "Target with codeword " << cw << " not found in gene panel. Returning " << NoName << " as Id\n";
  return(NoName);
 }
 return(panel[i].GetName());
}

int GenePanel::GetCodewordFromName(string genename)
{
 size_t i=0;
 while (i<panel.size() && panel[i].GetName()!=genename)
  i++;
 if (i==panel.size())
 {
  cerr << "Target with name " << genename << " not found in gene panel. Returning " << NoCode << " as gene code\n";
  return(NoCode);
 }
 return(panel[i].GetCodeword());
}

void GenePanel::SetColorOfGene(int cw,RGB_color &col)
{
 size_t i=0;
 while (i<panel.size() && panel[i].codeword!=cw)
  i++;
 if (i>=panel.size())
  ERRPROG("GenePanel: gene with codeword " << cw << " not found when trying to set its color.\n";)
 panel[i].color=col;
}

void GenePanel::FillColors(vector<ColortabEntry> &v)
{
 string gname;
 for (vector<ColortabEntry>::iterator it=v.begin(); it!=v.end(); ++it)
 {
  gname=GetNameFromCodeword(it->gene_codeword);
  if (gname==NoName)
   ERRPROG("FillColors: error at line " << it->line_in_file << " of the color file. " << it->gene_codeword << " is not a valid gene codeword.\n";)
  if (gname!=it->gene_name)
   ERRPROG("FillColors: error at line " << it->line_in_file << " of the color file. " << it->gene_name << " is not the gene name given in the gene table, which is " << gname << ".\n";)
  SetColorOfGene(it->gene_codeword,it->col);
 }
}

void GenePanel::SaveAsBinaryWithDebug(string fname)
{
 if (DEB)
  cout << "Writing panel as binary file " << fname << ".\n";

 OutputBufferedFile f(fname,default_chunksize_for_GenePanel,true);
 f.writeb(mincod);
 f.writeb(maxcod);
 f.writeb(mincov);
 f.writeb(maxcov);

 size_t ps=panel.size();
 f.writeb(ps);

 for (size_t i=0; i<ps; i++)
  panel[i].SaveAsBinary(f);
}

void GenePanel::SaveAsBinary(string fname)
{
 if (DEB)
  cout << "Writing panel as binary file " << fname << ".\n";

 OutputBufferedFile f(fname,default_chunksize_for_GenePanel,false);
 f.writeb(mincod);
 f.writeb(maxcod);
 f.writeb(mincov);
 f.writeb(maxcov);

 size_t ps=panel.size();
 f.writeb(ps);

 for (size_t i=0; i<ps; i++)
  panel[i].SaveAsBinary(f);
}
