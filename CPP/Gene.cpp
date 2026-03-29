#include "Gene.h"

using namespace std;

extern string progname;
extern bool DEB;

Gene::Gene(int cw,int gc,std::string cat,std::string gid,std::string gname,std::string gd,bool &err)
{
 if ((cw<0) || (cw>max_codeword))
 {
  WARNPROG("Incorrect codeword for a gene. It must be an integer in [0.." << max_codeword << "]\n";)
  err=true;
  return;
 }
 codeword=cw;
 
 if (gc<0)
 {
  WARNPROG("Incorrect coverage for a gene. It must be an integer >=0.\n";)
  err=true;
  return;
 }
 gene_coverage=gc;
 
 int i=0;
 while ((i<NumCat) && (catnames[i]!=cat))
  i++;
 if (i>NumCat)
 {
  WARNPROG("Incorrect category for a gene.\n";)
  err=true;
  return;
 }
 gene_cat=i;
 
 if ((gid.substr(0,3)!="ENS") && (gid.substr(0,3)!="Clo") && (gid.substr(0,3)!="Neg"))
 {
  WARNPROG("Incorrect identifier for a gene. It must start by 'ENS', 'Clo' or 'Neg'\n";)
  err=true;
  return;
 }
 
 id=gid;
 
 name=gname; 
 
 i=0;
 while ((i<NumDesc) && (descnames[i]!=gd))
  i++;
 if (i>NumDesc)
 {
  WARNPROG("Incorrect descriptor for a gene.\n";)
  err=true;
  return;
 }
 gene_desc=i; 
 
 if ((id.substr(0,3)=="Neg") && (descnames[gene_desc]!="negative_control"))
 {
  WARNPROG("A gene id starting by Neg is not described as 'negative_control'\n";)
  err=true;
  return;
 }
 
 color=black;

 err=false;
}

Gene::Gene(InputBufferedFile &f)
{
 f.readb(codeword);
 f.readb(gene_coverage);
 f.readb(gene_cat);
 f.readb(id);
 f.readb(name);
 f.readb(gene_desc);
}

void Gene::SaveAsBinary(OutputBufferedFile &f)
{
 f.writeb(codeword);
 f.writeb(gene_coverage);
 f.writeb(gene_cat);
 f.writeb(id);
 f.writeb(name);
 f.writeb(gene_desc);
}

string Gene::GetAll(char sep)
{
 stringstream st;
 st << codeword << sep << gene_coverage << sep << catnames[gene_cat] << sep << id << sep << name << sep << descnames[gene_desc];
 return(st.str());
}
