#include "Transcript.h"

using namespace std;

extern string progname;
extern bool DEB;

bool Transcript::ValidRegionName(string rname)
{
 size_t l=rname.length();
 bool OKname;
 switch (l)
 {
  case 2: OKname  = ((rname[0]>='A') && (rname[0]<='Z'));
          OKname &= ((rname[1]>='0') && (rname[1]<='9'));
          break;
  case 3: OKname  = ((rname[0]>='A') && (rname[0]<='Z'));
          OKname &= ( ((rname[1]>='A') && (rname[1]<='Z')) || ((rname[1]>='0') && (rname[1]<='9')) );
          OKname &= ((rname[2]>='0') && (rname[2]<='9'));
          break;
  case 4: OKname  = ((rname[0]>='A') && (rname[0]<='Z'));
          OKname &= ((rname[1]>='A') && (rname[1]<='Z'));
          OKname &= ((rname[2]>='0') && (rname[2]<='9'));
          OKname &= ((rname[3]>='0') && (rname[3]<='9'));
          break;
  default: OKname=false; break;                  
 }
 return(OKname);
}

Transcript::Transcript(unsigned long long idp,std::string cell_idp,bool overlap_nucp,std::string feat_namep,float xp,float yp,float zp,int qvp,std::string fov_namep,float nuc_distp,int codewordp,bool &err)
{
 id=idp;
 cell_id=cell_idp;
 overlap_nuc=overlap_nucp;
 feat_name=feat_namep;
 x=xp;
 y=yp;
 z=zp;
 qv=qvp;
 fov_name=fov_namep;
 nuc_dist=nuc_distp;
 codeword=codewordp;
 
 // Some basic checks;
 if (x<0.0 || y<0.0 || z<0.0)
 {
  WARNPROG("A transcript coordinate x, y, z (here, " << x << "," << y << "," << z << " respectively) or several of them is/are negative.\n";)
  err=true;
 }
 
 if (qv<0.0 || qv>40.0)
 {
  WARNPROG("A transcript quality index qv is " << qv << ", which is not in [0.0 .. 40.0].\n";)
  err=true;
 }
 
 if (!ValidRegionName(fov_name))
 {
  WARNPROG("The transcript FOV name " << fov_name << " is invalid (it should be LD or LLD or LDD or LLDD where L is an uppercase letter and D a digit).\n";)
  err=true;
 }
 
}

Transcript::Transcript(InputBufferedFile &f)
{
 f.readb(id);
 f.readb(cell_id);
 f.readb(overlap_nuc);
 f.readb(feat_name);
 f.readb(x);
 f.readb(y);
 f.readb(z);
 f.readb(qv);
 f.readb(fov_name);
 f.readb(nuc_dist);
 f.readb(codeword);
}

void Transcript::SaveAsBinary(OutputBufferedFile &f)
{
 f.writeb(id);
 f.writeb(cell_id);
 f.writeb(overlap_nuc);
 f.writeb(feat_name);
 f.writeb(x);
 f.writeb(y);
 f.writeb(z);
 f.writeb(qv);
 f.writeb(fov_name);
 f.writeb(nuc_dist);
 f.writeb(codeword);
}

string Transcript::GetAll(char sep)
{
 stringstream st;
 st << id << sep << cell_id << sep << (overlap_nuc ? "1" : "0") << sep << feat_name << sep << x << sep << y << sep << z << sep << qv << sep << fov_name << sep << nuc_dist << sep << codeword;
 return(st.str());
}
