#include "Transcriptset.h"

using namespace std;

extern string progname;
extern bool DEB;

TranscriptSet::TranscriptSet(std::string fname,bool check)
{
 if (fname.find(".csv")!=string::npos || fname.find(".CSV")!=string::npos)
 {
  if (DEB)
   cout << "Reading transcript set from .csv file " << fname << ".\n";
  ifstream f(fname.c_str());
  if (!f.is_open())
   ERRPROG("Cannot open transcript file " << fname << " to read. Check for existence and permissions.\n";)
  TranscriptSetFromCSV(fname);
 }
 else
 {
  if (DEB)
   cout << "Reading transcript set from binary file " << fname << ".\n";
  TranscriptSetFromBin(fname);
 }

 if (check)
  CheckTranscriptSet();
}

void TranscriptSet::TranscriptSetFromCSV(std::string fname)
{
 unsigned long line=0;
 unsigned long long idp;
 string cell_idp;
 bool overlap_nucp;
 string feat_namep;
 float xp,yp,zp;
 float qvp;
 string fov_namep;
 float nuc_distp;
 int codewordp;
 bool err;

 float cround=1;
 for (int i=0; i<num_dec_places; i++)
  cround *= 10;


 ifstream f(fname.c_str());
 // First, read the header line
 string trans_entry;
 getline(f,trans_entry);
 if (DEB)
 {
  cout << "Header from file " << fname << ":\n";
  cout << trans_entry << endl;
  cout << "Reading transcripts. Progress in millions. Please, wait...\n";
 }
 line++;

 // Then, start reading content lines
 getline(f,trans_entry);
 line++;
 while (!f.eof())
 {
  vector<string_view> vs=parseCSVRow(trans_entry,",");

  //for (size_t i=0; i<vs.size(); i++)
  // cout << "vs[" << i << "]=" << vs[i] << endl;

  if (vs.size()!=11)
   ERRPROG("Error at line " << line << " of input file " << fname << ". It hasn't 11 fields but " << vs.size() << ". Line is:\n" << trans_entry << endl;)

  idp=atol(vs[0].data());
  cell_idp=vs[1];
  overlap_nucp=(atoi(vs[2].data())==1);
  feat_namep=vs[3];
  xp=atof(vs[4].data());
  yp=atof(vs[5].data());
  zp=atof(vs[6].data());
  qvp=atof(vs[7].data());
  fov_namep=vs[8];
  nuc_distp=atof(vs[9].data());
  codewordp=atoi(vs[10].data());

  // Round to num_dec_places the position and distance
  xp=floor(cround*xp+0.5)/cround;
  yp=floor(cround*yp+0.5)/cround;
  zp=floor(cround*zp+0.5)/cround;
  qvp=floor(cround*qvp+0.5)/cround;
  nuc_distp=floor(cround*nuc_distp+0.5)/cround;

  err=false;
  Transcript t(idp,cell_idp,overlap_nucp,feat_namep,xp,yp,zp,qvp,fov_namep,nuc_distp,codewordp,err);
  if (err)
   WARNPROG("Error at line " << line << " of input file " << fname << ". Line is:\n" << trans_entry << endl;)
  else
   tset.push_back(t);

  if (DEB && !(line%1000000))
  {
   cout << line/1000000 << "M ";
   cout.flush();
  }
  getline(f,trans_entry);
  if (!f.eof())
   line++;
 }
 if (DEB)
  cout << endl;

 f.close();

 if (DEB)
  cout << tset.size() << " transcripts read.\n";
}

void TranscriptSet::TranscriptSetFromBin(std::string fname)
{
 InputBufferedFile f(fname);

 if (DEB)
  cout << "Reading transcripts. Progress in millions. Please, wait...\n";
 size_t ps=0;
 f.readb(ps);
 for (size_t i=0; i<ps; i++)
 {
  if (DEB && (i>0) && !(i%1000000))
  {
   cout << i/1000000 << "M ";
   cout.flush();
  }
  Transcript t(f);
  tset.push_back(t);
 }

 if (DEB)
  cout << endl << tset.size() << " items read.\n";
}

bool TranscriptSet::CheckTranscriptSet()
{
 if (DEB)
 {
  cout << "Checking " << tset.size() << " transcripts...";
  cout.flush();
 }
 for (size_t i=0; i<tset.size()-1; i++)
 {
  if (DEB && !(i%1000))
  {
   cout << i << " ";
   cout.flush();
  }
  for (size_t j=i+1; j<tset.size(); j++)
   if (tset[i].GetId()==tset[j].GetId())
    ERRPROG("Duplicated transcript identifier (" << tset[i].GetId() << ") between items " << i+1 << " and " << j+1 << endl;)
 }
 if (DEB)
  cout << " Done. No duplicated identifiers found.\n";

 return(true);
}

void TranscriptSet::Show(ostream &out,char sep)
{
 for (size_t g=0; g<tset.size(); g++)
  out << tset[g].GetAll(sep) << endl;
}

void TranscriptSet::SaveAsBinaryWithDebug(string fname)
{
 OutputBufferedFile f(fname,default_chunksize_for_TranscriptSet,true);
 if (DEB)
  cout << "Writing transcript set as binary file " << fname << ".\n";

 size_t ps=tset.size();
 f.writeb(ps);

 for (size_t i=0; i<ps; i++)
  tset[i].SaveAsBinary(f);
}

void TranscriptSet::SaveAsBinary(string fname)
{
 OutputBufferedFile f(fname,default_chunksize_for_TranscriptSet,false);
 if (DEB)
  cout << "Writing transcript set as binary file " << fname << ".\n";

 size_t ps=tset.size();
 f.writeb(ps);

 for (size_t i=0; i<ps; i++)
  tset[i].SaveAsBinary(f);
}

