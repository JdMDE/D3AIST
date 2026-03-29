#include "buildim.h"

using namespace std;

string progname="buildim";
bool DEB=true;

void Usage()
{
 cerr << "Usage:\n   " << progname << " <task_desc_file>\n";
 exit(1);
}

string SkipComments(ifstream &f,int &lnum)
{
 string line;
 do
 {
  if (!f.eof())
  {
      getline(f,line);
      lnum++;
  }
 }
 while (!f.eof() && (line=="" || line[0]=='#'));
 return line;
}

void Accessible(string fname)
{
 ifstream f(fname.c_str());
 if (!f.is_open())
  ERRPROG("File " << fname << " mentioned in the description file cannot be opened. Check for existance and permissions.\n";)
}

void ParseColorTable(string fcolname,vector<ColortabEntry> &v)
{
 ifstream f(fcolname.c_str());
 if (!f.is_open())
  ERRPROG("ParseColorTable: cannot open file " << fcolname << " mentioned in the configuration file to read the color table from it.\n";)
 int lnum=1;
 string line;
 string cwn,name,rn,gn,bn;
 do
 {
  if (!f.eof())
  {
   getline(f,line);

   if (line!="")
   {
    std::replace(line.begin(),line.end(),',',' ');
    cwn=name=rn=gn=bn="";
    istringstream ss(line);
    ss >> cwn >> name >> rn >> gn >> bn;

    if ((cwn=="") || (name=="") || (rn=="") || (gn=="") || (bn==""))
     ERRPROG("ParseColorTable: error reading line " << lnum << " of color file " << fcolname << ". Incorrect format (it should be 'codeword,name,r,g,b')\n";)

    ColortabEntry cte;
    cte.line_in_file=lnum;
    cte.gene_codeword=atoi(cwn.c_str());
    cte.gene_name=name;
    cte.col.r=atoi(rn.c_str());
    cte.col.g=atoi(gn.c_str());
    cte.col.b=atoi(bn.c_str());
    v.push_back(cte);

    lnum++;
   }
  }
 }
 while (!f.eof());
}

bool FindParName(string pname,int &npar,ParamTypes &what_expected,bool &alternative)
{
 size_t i=0;
 bool found=false;
 while (i<CurrentNumberParameters && !found)
 {
  if (pname==Params[i].ParamName)
  {
   npar=i;
   what_expected=Params[i].ParType;
   found=true;
   alternative=false;
  }
  if (pname==Params[i].AltParamName)
  {
   npar=i;
   what_expected=Params[i].AltParType;
   found=true;
   alternative=true;
  }
  i++;
 }

 return(found);
}

string DivideLine(string line,int lnum,string &values)
{
 size_t i=0;
 while ( i<line.size() && ((line[i]==' ') || (line[i]=='\t')) )
  i++;
 if (i==line.size())
  ERRPROG("DivideLine: line " << lnum << " of configuration file seems to be empty...\n";)

 string parname="";
 while ( i<line.size() && line[i]!=':' )
 {
  parname.push_back(line[i]);
  i++;
 }
 if (i==line.size())
  ERRPROG("DivideLine: line " << lnum << " of configuration file does not contain a semmicolon after the parameter name.\n";)
 i++;

 while ( i<line.size() && ((line[i]==' ') || (line[i]=='\t')) )
  i++;
 if (i==line.size())
  ERRPROG("DivideLine: line " << lnum << " of configuration file gives no value to parameter " << parname << ".\n";)

 values=line.substr(i);
 return(parname);
}

vector<string> ListString(string args)
{
 vector<string> ret;
 size_t i=0;
 string sn;
 while (i<args.size())
 {
  size_t j=i;
  sn="";
  while (j<args.size() && args[j]!=',' && args[j]!=' ')
  {
   sn.push_back(args[j]);
   j++;
  }
  ret.push_back(sn);
  while (j<args.size() && (args[j]==',' || args[j]==' '))
   j++;
  i=j;
 }
 return(ret);
}

void LexicalAnalysis(string dfname,array<pair<string,bool>,CurrentNumberParameters> &values)
{
 ifstream dfile(dfname.c_str());
 if (!dfile.is_open())
  ERRPROG("ParseDescFile: Description file " << dfname << " cannot be opened. Check existence and permissions.\n";)

 // Line number to emit errors
 int lnum=0;

 // Arrays to know if a parameter has been used before
 array<bool,CurrentNumberParameters> already_found;
 fill(already_found.begin(),already_found.end(),false);

 // Skip initial comment lines, if any
 string line=SkipComments(dfile,lnum);

 do
 {
  // The parameter name in a line, and the arguments to it
  string args;

  // This returns the parameter name (without the final semmicolon) and the arguments (without initial blank spaces, if any)
  string par_name=DivideLine(line,lnum,args);

  // Variables to be filled by FindParName
  int pnum;
  ParamTypes ptype;
  bool alt;

  // Error processing.
  // Error type I: parameter is not in the list of known parameters. FindParName returns false.
  if (!FindParName(par_name,pnum,ptype,alt))
   ERRPROG("ParseDescFile: error in line " << lnum << " of configuration file.\nIncorrect parameter " << par_name << ".\n";)
  // Error type II: parameter had been given before....
  if (already_found[pnum]==true)
  {
   // Type IIa: ... becuase it has been repeated...
   if (Params[pnum].ParType==ParamTypes::NoType)
    ERRPROG("ParseDescFile: error in line " << lnum << " of configuration file.\nRepeated parameter " << par_name << ".\n";)
   else
   // Type IIb: ... or because the alternative form has been used before.
    ERRPROG("ParseDescFile: error in line " << lnum << " of configuration file.\nUse of parameter " << par_name << " after having used " << ((alt) ? Params[pnum].ParamName : Params[pnum].AltParamName) << ".\nBoth cannot be set simultaneously.\n";)
  }

  // Parameter is correct. Fill the information arrays with the appropriate values
  values[pnum]=pair(args,alt);
  already_found[pnum]=true;

  // Go to next not-comment line
  line=SkipComments(dfile,lnum);
 }
 while (!dfile.eof());
 dfile.close();
}

string PrependPath(string path,string prpath)
{
 if (prpath=="" || path[0]==dirslash[0])
  return(path);     // We don't prepend anything either if we don't want to (empty prpath) or if the path is absolute
 else
  return(prpath+path);
}

void ParseDescFile(string dfname,string &genepanelname,string &cname,string &cbname,string &ncbname,string &tname,string &clname,string &outrname,string &outimname,string &outhistname,ImParams &ipar)
{
 array<pair<string,bool>,CurrentNumberParameters> values;
 values.fill(pair("",false));

 LexicalAnalysis(dfname,values);

 // The prepend path for all file names. It might be empty if file names are given as absolute paths or paths relative to the current working directory
 string prpath="";

 bool cells_read_as_bin=false;

 for (size_t par=0; par<CurrentNumberParameters; par++)
 {
  string v=values[par].first;
  ParamTypes ptype = (values[par].second ? Params[par].AltParType : Params[par].ParType);
  istringstream stline(v);
  string vs="";
  float vf[4];
  vf[0]=numeric_limits<float>::quiet_NaN();
  int vi[2];
  vi[0]=numeric_limits<int>::quiet_NaN();
  vector<int> vil;
  vector<string> vsl;

  if (v=="")
    vs="";
  else
  {
   switch (ptype)
   {
    case ParamTypes::Boolean:
      vs=v;
      transform(v.begin(),v.end(),vs.begin(),::tolower);
      if (vs.find("true")!=string::npos)
       vs="true";
      else
       if (vs.find("false")!=string::npos)
        vs="false";
       else
        ERRPROG("ParseDescFile: error in configuration file.\nPossible values for " << Params[par].ParamName << " are only true or false.\n";)
      break;
    case ParamTypes::Float2:
      stline >> vf[0] >> vf[1];
      break;
    case ParamTypes::Float4:
      stline >> vf[0] >> vf[1] >> vf[2] >> vf[3];
      break;
    case ParamTypes::Int:
      stline >> vi[0];
      break;
    case ParamTypes::Int2:
      stline >> vi[0] >> vi[1];
      break;
    case ParamTypes::ListInt:
      vsl=ListString(v);
      for (vector<string>::iterator it=vsl.begin(); it!=vsl.end(); ++it)
       vil.push_back(atoi(it->c_str()));
      break;
    case ParamTypes::ListString:
      vsl=ListString(v);
      break;
    case ParamTypes::String:
      stline >> vs;
      break;
    case ParamTypes::NoType:
      ERRPROG("ParseDescFile: error in configuration file.\nParameter " << Params[par].ParamName << " has no type.\n";)
      break;
   }
  }

  if (par==0)  // Prepend path. This parameter is analyzed previously to any other. That's why it is in the 0-place of parameter array
  {
   if (vs!="none" && vs!="")
   {
    if (vs[0]!=dirslash[0])
     ERRPROG("ParseDescFile: parameter PrependPath must be an absolute path (i.e.: start with '" << dirslash << "'). You have given " << vs << ".\n";)
    if (vs.back() != dirslash[0])
     vs.push_back(dirslash[0]);
    prpath=vs;
   }
   else
    prpath="";
  }

  if (par==1)  // Genes
  {
   if (vs=="")
    ERRPROG("ParseDescFile: parameter Genes is compulsory.\n";)
   genepanelname = PrependPath(vs,prpath);
   Accessible(genepanelname);
  }

  if (par==2) // Cells
  {
   if (vs=="")
    ERRPROG("ParseDescFile: parameter Cells is compulsory.\n";)
   cname=PrependPath(vs,prpath);
   Accessible(cname);
   if (cname.find(".bin")!=string::npos || cname.find("BIN")!=string::npos)
    cells_read_as_bin=true;
  }

  if (par==3) // Transcripts
  {
   if (vs=="")
    ERRPROG("ParseDescFile: parameter Transcripts is compulsory.\n";)
   tname=PrependPath(vs,prpath);
   Accessible(tname);
  }

  if (par==4)  // Clusters. It could be none or not mentioned.
  {
   if (vs=="none" || vs=="")
    clname="none";
   else
   {
    clname = PrependPath(vs,prpath);
    Accessible(clname);
   }
  }

  if (par==5) // Cell boundaries. It could be none or not mentioned.
  {
   if (vs=="none" || vs=="")
    cbname="none";
   else
   {
    if (cells_read_as_bin)
    {
     WARNPROG("ParseDescFile: cell set seems to be a binary file. Therefore, it should contain the cell boundaries so you can't read them from a .csv file.\n";)
     cbname="none";
    }
    else
    {
     cbname = PrependPath(vs,prpath);
     Accessible(cbname);
    }
   }
  }

  if (par==6) // Nucleus boundaries. It could be none or not mentioned.
  {
   if (vs=="none" || vs=="")
    ncbname="none";
   else
   {
    if (cells_read_as_bin)
    {
     WARNPROG("ParseDescFile: cell set seems to be a binary file. Therefore, it should contain the nucleus boundaries so you can't read them from a .csv file.\n";)
     ncbname="none";
    }
    else
    {
     ncbname = PrependPath(vs,prpath);
     Accessible(ncbname);
    }
   }
  }

  if (par==7) // Spacing.
  {
   if (isnan(vf[0]))
    ERRPROG("ParseDescFile: parameter Spacing is compulsory.\n";)
   if (vf[0]<=0.0 || vf[1]<=0.0)
    ERRPROG("Spacing given in description file does not seem to be a pair of possitive real numbers.\n";)
   ipar.osp[0]=vf[0];
   ipar.osp[1]=vf[1];
  }

  if (par==8) // Dimensions
  {
   if (isnan(vi[0]))
    ERRPROG("ParseDescFile: parameter Dimensions is compulsory.\n";)
   if (vi[0]<=0 || size_t(vi[0])>MaxImWidth || vi[1]<0 || size_t(vi[1])>MaxImHeight)
    ERRPROG("Image dimensions given in description file does not seem to be a pair of possitive integer numbers (or they are too big).\n";)
   ipar.w=vi[0];
   ipar.h=vi[1];
  }

  if (par==9) // PixPerBin
  {
   if (isnan(vi[0]))
    ERRPROG("ParseDescFile: parameter PixPerBin is compulsory.\n";)
   if (vi[0]<0 || vi[0]>MaxPPB)
    ERRPROG("Pixels per bin (macropixel) given in description file does not seem to be a possitive integer number or is bigger than " << MaxPPB << ".\n";)
   ipar.ppbin=vi[0];
  }

  if (par==10)
  {
   if (vs=="" && isnan(vf[0]))
     ERRPROG("ParseDescFile: one of the parameters AoIWin or AoIID is compulsory.\n";)

   if (values[par].second==true)  // AoIID
   {
    ipar.idfov=vs;
    ipar.exfov=true;
   }
   else                             // AoIWin
   {
    ipar.exfov=false;
    ipar.idfov="none";
    for (int i=0; i<4; i++)
    {
     if (vf[i]<0.0)
      ERRPROG("Value for at least one coordinate of the Area of interest window given in description file is negative.\n";)
     if ( (i%2==0 && vf[i]>float(ipar.w)*ipar.osp[0]) || (i%2 && vf[i]>float(ipar.h)*ipar.osp[1]) )
      ERRPROG("Value for at least one coordinate of the Area of interest window given in description file seem to be outside the image limits.\n";)
     ipar.fovwin[i]=vf[i];
    }
   }
  }

  if (par==11) // QvThres
  {
   if (isnan(vi[0]))
    ERRPROG("ParseDescFile: parameter QvThres is compulsory.\n";)
   if (vi[0]<0 || vi[0]>40)
    ERRPROG("Parameter QvThres given in description file must have an integer value in 1..40, not " << vi[0] << "\n";)
   ipar.qvthres=vi[0];
  }

  if (par==12) // OnlyTarget
   ipar.onlytargets = (v=="true" && vs!="");  // If not given, i.e. v is empty, the default value is false

  if (par==13)
  {
   if (vsl.empty() && vil.empty())
    ERRPROG("ParseDescFile: one of the parameters TargetNames or TargetIds is compulsory.\n";)

   if (values[par].second==true) // TargetNames
   {
    ipar.targetnames=vsl;
    ipar.numcodes.clear();
   }
   else                          // TargetIds
   {
    ipar.numcodes=vil;
    ipar.targetnames.clear();
   }
  }

  if (par==14) // GenMarks
   ipar.gmarks = (vs=="true" && v!="");  // If not given, i.e. v is empty, the default value is false

  if (par==15) // OutRfile. It could be none or not mentioned.
  {
   if (vs=="none" || vs=="")
    outrname="none";
   else
    outrname = PrependPath(vs,prpath);
  }

  if (par==16) // OutImfile. It could be none or not mentioned.
  {
   if (vs=="none" || vs=="")
    outimname="none";
   else
    outimname = PrependPath(vs,prpath);
  }

  if (par==17) // ColorFile. It could be none or not mentioned.
  {
   string colfname;
   if (vs=="none" || vs=="")
    colfname="none";
   else
   {
    colfname = PrependPath(vs,prpath);
    Accessible(colfname);
    ParseColorTable(colfname,ipar.colortab);
   }
  }

  if (par==18) // OutHistfile. It could be none or not mentioned.
  {
   if (vs=="none" || vs=="")
    outhistname="none";
   else
    outhistname = PrependPath(vs,prpath);
  }

  if (par==19)  // CheckThisFile.
   ipar.checkandexit = (vs=="true" && vs!="");  // If not given, i.e. v is empty, the default value is false
 } // end for

}

void CheckConfigFile(string idescname,string genepanelname,string cname,string cbname,string ncbname,string tname,string clname,string outrname,string outimname,string outhistname,ImParams const &ipar)
{
 cout << "=======================================================================\n";
 cout << "Check of the parameters set in description file " << idescname << ":\n";

 cout << "genepanelname: " << genepanelname << endl;
 cout << "cname        : " << cname << endl;
 cout << "tname        : " << tname << endl;
 cout << "clname       : " << clname << endl;
 cout << "cbname       : " << cbname << endl;
 cout << "ncbname      : " << ncbname << endl;
 cout << "outrname     : " << outrname << endl;
 cout << "outimname    : " << outimname << endl;
 cout << "outhistname  : " << outhistname << endl;
 cout << "ipar.osp        : {" << ipar.osp[0] << "," << ipar.osp[1] << "}\n";
 cout << "ipar.ppbin      : " << ipar.ppbin << endl;
 cout << "ipar.w          : " << ipar.w << endl;
 cout << "ipar.h          : " << ipar.h << endl;
 cout << "ipar.exfov      : " << (ipar.exfov ? "true" : "false") << endl;
 cout << "ipar.idfov      : " << ipar.idfov << endl;
 cout << "ipar.fovwin     : {[" << ipar.fovwin[0] << "," << ipar.fovwin[1] << "][" << ipar.fovwin[2] << "," << ipar.fovwin[3] << "]}\n";
 cout << "ipar.onlytargets: " << (ipar.onlytargets ? "true" : "false") << endl;
 cout << "ipar.qvthres    : " << ipar.qvthres << endl;;
 cout << "ipar.numcodes   : ( ";
 if (ipar.numcodes.empty())
  cout << "empty list ";
 else
  for (size_t i=0; i<ipar.numcodes.size(); i++)
   cout << ipar.numcodes[i] << " ";
 cout << ")\n";
 cout << "ipar.targetnames: ( ";
 if (ipar.targetnames.empty())
  cout << "empty list ";
 for (size_t i=0; i<ipar.targetnames.size(); i++)
  cout << ipar.targetnames[i] << " ";
 cout << ")\n";
 cout << "ipar.gmarks     : " << (ipar.gmarks ? "true" : "false") << endl;
 cout << "ipar.colortab   :";
 if (!ipar.colortab.empty())
 {
  cout << " entries set for the following " << ipar.colortab.size() << " items in the form [line in file, GeneID, GeneName, (R,G,B)]\n";
  for (size_t i=0; i<ipar.colortab.size(); i++)
  {
   cout << " [" << ipar.colortab[i].line_in_file << "," << ipar.colortab[i].gene_codeword << "," << ipar.colortab[i].gene_name << ",";
   cout << "(" << int(ipar.colortab[i].col.r) << " " << int(ipar.colortab[i].col.g) << " " << int(ipar.colortab[i].col.b) << ")]\n";
  }
 }
 else
  cout << " ( empty tab )\n";
 cout << "=======================================================================\n";
 exit(0);
}

int main(int argc,char *argv[])
{
 progname=string(argv[0]);

 if (argc!=2)
  Usage();
  
 string idescname=string(argv[1]);
 string genepanelname,cname,tname,clname,cbname,ncbname,outrname,outimname,outhistname;
 ImParams ipar;

 ParseDescFile(idescname,genepanelname,cname,cbname,ncbname,tname,clname,outrname,outimname,outhistname,ipar);

 if (ipar.checkandexit)
  CheckConfigFile(idescname,genepanelname,cname,cbname,ncbname,tname,clname,outrname,outimname,outhistname,ipar);

 if (DEB)
  cout << "=====================================================================\nReading gene panel from file " << genepanelname << endl;
 GenePanel panel(genepanelname);
 if (!ipar.colortab.empty())
  panel.FillColors(ipar.colortab);

 /*
 for (size_t i=panel.GetMinCodeword(); i<=panel.GetMaxCodeword(); i++)
 {
  string name=panel.GetNameFromCodeword(i);
  cout << i << " " << panel.GetIdFromCodeword(i) << " " << name << " " << panel.GetCodewordFromName(name) << endl;
 }
 exit(1);
 */

 if (DEB)
  cout << "=====================================================================\nReading cells from file " << cname << endl;
 CellSet cs(cname);
 
 if (DEB)
  cout << "=====================================================================\nReading transcripts from file " << tname << endl;
 TranscriptSet ts(tname);

 if (clname!="none")
 {
  if (DEB)
   cout << "=====================================================================\nReading clusters from file " << clname << endl;
  cs.AddClusterInfo(clname);
 }

 if (cbname!="none")
 {
  if (DEB)
   cout << "=====================================================================\nReading cell boundaries from file " << cbname << endl;
  cs.AddBoundaryInfo(cbname,true);
 }

 if (ncbname!="none")
 {
  if (DEB)
   cout << "=====================================================================\nReading nucleus boundaries from file " << ncbname << endl;
  cs.AddBoundaryInfo(ncbname,false);
 }

 if (outrname != "none" || outimname != "none" || outhistname != "none")
 {
  if (DEB)
   cout << "======================================================================\nCreating pseudoimage\n";
  SpRData Spr(ipar,panel,cs,ts);

  if (DEB)
   cout << "======================================================================\nBuilding pseudoimage\n";

  if (outhistname == "none")
   Spr.Build();
  else
   Spr.Build(outhistname);

  if (outrname != "none")
  {
   if (DEB)
   {
    string bname=outrname.substr(0,outrname.find_last_of("."))+"xxx.csv";
    cout << "======================================================================\nSaving data in files " << bname << endl;
   }
   Spr.Save(outrname);
  }

  if (outimname != "none")
  {
   string binimname=outimname.substr(0,outimname.find_last_of("."))+".tif";
   if (DEB)
    cout << "=======================================================================\nSaving binary image in file " << binimname << endl;
   Spr.SaveAsImages(binimname,idescname);
  }

 }

 if (DEB)
  cout << "******** Whole process successfully finished *********\n";

 return(0);
}

