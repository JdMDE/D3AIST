#include "Cellset.h"

using namespace std;

extern string progname;
extern bool DEB;

CellSet::CellSet(string fname,bool check)
{
 if (fname.find(".csv")!=string::npos || fname.find(".CSV")!=string::npos)
 {
  ifstream f(fname.c_str());
  if (!f.is_open())
   ERRPROG("Cannot open cell file " << fname << " to read. Check for existence and permissions.\n";)
  f.close();
  CellSetFromCSV(fname);
  if (DEB)
   cout << cset.size() << " cells read.\n";
 }
 else
 {
  CellSetFromBin(fname);
  if (DEB)
  {
   cout << cset.size() << " cells read.\n";
   if (num_cells_with_boundary>0)
    cout << num_cells_with_boundary << " cells whose boundary has been read.\n";
   if (num_nuc_with_boundary>0)
    cout << num_nuc_with_boundary << " cells whose nucleus boundary has been read.\n";
  }
 }

 FillId2place();

 if (check)
  CheckCellSet();
}

void CellSet::CellSetFromCSV(string fname)
{
 unsigned long line=0;
 string idp;
 float xp,yp,areap,nuc_areap;
 int trcp,cpcp,cccp,uccp,dccp,tcp;
 bool err;

 float cround=1;
 for (int i=0; i<num_dec_places; i++)
  cround *= 10;

 ifstream f(fname.c_str());
  // First, read the header line
 string cell_entry;
 getline(f,cell_entry);
 if (DEB)
 {
  cout << "Header from file " << fname << ":\n";
  cout << cell_entry << endl;
 }
 line++;

 // Then, start reading content lines
 getline(f,cell_entry);
 line++;
 size_t pos=0;
 while (!f.eof())
 {
  vector<string_view> vs=parseCSVRow(cell_entry,",");

  if (vs.size()!=11)
   ERRPROG("Error at line " << line << " of input file " << fname << ". It hasn't 11 fields but " << vs.size() << ". Line is:\n" << cell_entry << endl;)

  idp=vs[0];
  xp=atof(vs[1].data());
  yp=atof(vs[2].data());
  trcp=atoi(vs[3].data());
  cpcp=atoi(vs[4].data());
  cccp=atoi(vs[5].data());
  uccp=atoi(vs[6].data());
  dccp=atoi(vs[7].data());
  tcp=atoi(vs[8].data());
  areap=atof(vs[9].data());
  nuc_areap=atof(vs[10].data());

  //st >> idp >> xp >> yp >> trcp >> cpcp >> cccp >> uccp >> dccp >> tcp >> areap >> nuc_areap;

  // Round to num_dec_places the position and area values.
  xp=floor(cround*xp+0.5)/cround;
  yp=floor(cround*yp+0.5)/cround;
  areap=floor(cround*areap+0.5)/cround;
  nuc_areap=floor(cround*nuc_areap+0.5)/cround;

  err=false;
  Cell c(idp,xp,yp,trcp,cpcp,cccp,uccp,dccp,tcp,areap,nuc_areap,err);
  if (err)
   ERRPROG("Error at line " << line << " of input file " << fname << ". Line is:\n" << cell_entry << endl;)

  id2place[idp]=pos;
  cset.push_back(c);
  pos++;

  getline(f,cell_entry);
  if (!f.eof())
   line++;
 }

 f.close();
}

void CellSet::CellSetFromBin(string fname)
{
 InputBufferedFile f(fname);

 size_t ps=0;
 f.readb(ps);
 num_cells_with_boundary=num_nuc_with_boundary=0;
 for (size_t i=0; i<ps; i++)
 {
  Cell c(f,num_cells_with_boundary,num_nuc_with_boundary);
  cset.push_back(c);
 }

}

Cell& CellSet::GetCellRef(std::string idc)
{
 if (id2place.contains(idc))
  return(cset[id2place[idc]]);
 else
  ERRPROG("GetCellRef: cannot get cell with identifier " << idc << ".\n";)
}

void CellSet::FillId2place()
{
 size_t pos=0;
 for (vector<Cell>::iterator it=cset.begin(); it!=cset.end(); it++)
 {
  id2place[it->GetId()]=pos;
  pos++;
 }
}

void CellSet::CheckCellSet()
{
 if (DEB)
 {
  cout << "Cheking " << cset.size() << " cells...";
  cout.flush();
 }
 for (size_t i=0; i<cset.size()-1; i++)
 {
  if (DEB && !(i%1000))
  {
   cout << i << " ";
   cout.flush();
  }
  for (size_t j=i+1; j<cset.size(); j++)
   if (cset[i].GetId()==cset[j].GetId())
    ERRPROG("Duplicated cell identifier (" << cset[i].GetId() << ") between items " << i+1 << " and " << j+1 << endl;)
 }
 if (DEB)
  cout << " Done. No duplicated identifiers found.\n";
}

void CellSet::Show(ostream &out,char sep,bool getcl)
{
 for (size_t c=0; c<cset.size(); c++)
  out << cset[c].GetAll(sep,getcl) << endl;
}

void CellSet::SetCellCluster(string idc,int clc,bool &err)
{
 unordered_map<string,size_t>::iterator it = id2place.find(idc);
 if (it==id2place.end())
 {
  WARNPROG("No cell with cell identifier " << idc << " currently stored in cell set.";)
  err=true;
  return;
 }
 cset[id2place[idc]].SetCluster(clc);
}

void CellSet::AddClusterInfo(string clfile)
{
 ifstream f(clfile.c_str());
 if (!f.is_open())
  ERRPROG("Cannot open cluster file " << clfile << " to read. Check for existence and permissions.\n";)

 unsigned long line=0;
 string cellid;
 int cellclust;

 string clust_entry;
 // First, read the header line
 getline(f,clust_entry);
 if (DEB)
 {
  cout << "Header from file " << clfile << ":\n";
  cout << clust_entry << endl;
 }
 line++;

 bool err;
 // Then, start reading content lines
 getline(f,clust_entry);
 line++;
 while (!f.eof())
 {
  vector<string_view> vs=parseCSVRow(clust_entry,",");

  if (vs.size()!=2)
   ERRPROG("Error at line " << line << " of input file " << clfile << ". It hasn't 2 fields but " << vs.size() << ". Line is:\n" << clust_entry << endl;)

  cellid=vs[0];
  cellclust=atoi(vs[1].data());

  err=false;
  SetCellCluster(cellid,cellclust,err);
  if (err)
   ERRPROG("Error at line " << line << " of input file " << clfile << ". Line is:\n" << clust_entry << endl;)

  getline(f,clust_entry);
  if (!f.eof())
   line++;
 }
 f.close();
}

void CellSet::SetBoundary(string idc, const vector<float>& bx, const vector<float>& by,const vector<size_t>& startp,bool cellboundary, bool& err)
{
 unordered_map<string,size_t>::iterator it = id2place.find(idc);
 if (it==id2place.end())
 {
  WARNPROG("No cell with cell identifier " << idc << " currently stored in cell set.";)
  err=true;
  return;
 }
 if (cellboundary)
 {
  cset[id2place[idc]].cboundx=bx;
  cset[id2place[idc]].cboundy=by;
 }
 else
 {
  cset[id2place[idc]].nuccboundx=bx;
  cset[id2place[idc]].nuccboundy=by;
  cset[id2place[idc]].startpoints=startp;
 }
}

void CellSet::AddBoundaryInfo(std::string bfile, bool cellboundary)
{
 ifstream f(bfile.c_str());
 if (!f.is_open())
  ERRPROG("Cannot open " << (cellboundary ? "cell " : "nucleus ") << "boundary file " << bfile << " to read. Check for existence and permissions.\n";)

 unsigned long line=0;
 string cellid,thiscellid;
 float x,y;
 long cellnum,thiscellnum;
 vector<float> boundaryx,boundaryy;
 vector<size_t> startpoints;

 string b_entry;
 // First, read the header line
 getline(f,b_entry);
 if (DEB)
 {
  cout << "Header from file " << bfile << ":\n";
  cout << b_entry << endl;
 }
 line++;

 bool err;
 // Then, start reading content lines
 cellnum=-1;         // This is to mark the abolute beginning
 boundaryx.clear();
 boundaryy.clear();

 getline(f,b_entry);
 line++;
 size_t num_cells=0;
 size_t num_segments=0;
 size_t thiscell_segments;
 size_t num_cells_multiple_bound=0;
 while (!f.eof())
 {
  vector<string_view> vs=parseCSVRow(b_entry,",");

  if (vs.size()!=4)
   ERRPROG("Error at line " << line << " of input file " << bfile << ". It hasn't 4 fields but " << vs.size() << ". Line is:\n" << b_entry << endl;)

  thiscellid=vs[0];
  x=atof(vs[1].data());
  y=atof(vs[2].data());
  thiscellnum=atoi(vs[3].data());

  // This happens only at the very beginning...
  if (cellnum==-1)
  {
   cellnum=thiscellnum;
   cellid=thiscellid;
   thiscell_segments=1;
  }

  // This happens when a new cell boundary, or a new segment of the nucleus cell boundary, starts.
  if (cellnum!=thiscellnum)
  {
   if (cellboundary)
   {
    if (cellid==thiscellid)  // In the case of cell boundaries, change of id and number must coincide so the new cell id must differ from the former one
     err=true;
    else
    {
     err=false;
     SetBoundary(cellid,boundaryx,boundaryy,startpoints,cellboundary,err);  // err will go back being false, unless the cell id does not identify any known cell

     // Reset the name and id to the new values....
     cellnum=thiscellnum;
     cellid=thiscellid;
     // and start a new boundary
     boundaryx.clear();
     boundaryy.clear();
     num_cells++;
     num_segments++;
    }
   }
   else  // Here we may be changing of cell, too, or not.
   {
    if (cellid!=thiscellid)   // cell has finished...
    {
     err=false;
     SetBoundary(cellid,boundaryx,boundaryy,startpoints,cellboundary,err);  // err will go back being false, unless the cell id does not identify any known cell

     // Reset the name and id to the new values....
     cellnum=thiscellnum;
     cellid=thiscellid;
     // and the number of segments o this cell also to 0
     num_segments += thiscell_segments;
     if (thiscell_segments>1)
      num_cells_multiple_bound++;
     thiscell_segments=1;
     // and start a new boundary
     boundaryx.clear();
     boundaryy.clear();
     startpoints.clear();
     startpoints.push_back(0);   // The first new point (position 0) will be the start of a new boundary segment, too
     num_cells++;
    }
    else                       // more segments of this cell remain...
    {
     cellnum=thiscellnum;  // Reset the number, but dont' reset the id, and don't start a new boundary...
     thiscell_segments++;
     startpoints.push_back(boundaryx.size());  // Take note of the point where the new boundary segment starts
    }
   }
  }
  // The point will be added either to the current boundary (if no change of cell has happened) or as the first point of the new cell
  boundaryx.push_back(x);
  boundaryy.push_back(y);

  if (err)
   ERRPROG("Error at line " << line << " of input file " << bfile << ". Line is:\n" << b_entry << endl;)

  getline(f,b_entry);
  if (!f.eof())
   line++;
 }
 f.close();

 // We went out of the loop because the file ended, but this means the last cell/nucleus had finished and has not been added to the set, so let's do it...
 err=false;
 SetBoundary(cellid,boundaryx,boundaryy,startpoints,cellboundary,err);
 // ... and let's increment the cell count, too.
 num_cells++;

 if (DEB)
 {
  cout << "Found " << (cellboundary ? "cell " : "nucleus ") << "boundaries of " << num_cells << " cells.\n";
  if (num_cells<cset.size())
    WARNPROG("AddBoundaryInfo: there is no " << (cellboundary ? "cell " : "nucleus ") << "boundary for " << cset.size()-num_cells << " cells.\n";)
  if (num_cells_multiple_bound>0)
    cout << num_cells_multiple_bound << " cells have nucleus with multiple segment boundaries.\n";
 }
}

void CellSet::SaveAsBinary(string fname)
{
 OutputBufferedFile f(fname,default_chunksize_for_CellSet,false);
 if (DEB)
  cout << "Writing cell set as binary file " << fname << ".\n";

 size_t ps=cset.size();
 f.writeb(ps);

 for (size_t i=0; i<ps; i++)
  cset[i].SaveAsBinary(f);
}

void CellSet::SaveAsBinaryWithDebug(string fname)
{
 OutputBufferedFile f(fname,default_chunksize_for_CellSet,true);
 if (DEB)
  cout << "Writing cell set as binary file " << fname << ".\n";

 size_t ps=cset.size();
 f.writeb(ps);

 for (size_t i=0; i<ps; i++)
  cset[i].SaveAsBinary(f);
}
