#include "SpRData.h"

using namespace std;
using namespace itk;

extern string progname;
extern bool DEB;

SpRData::SpRData(ImParams &ipar,GenePanel &panelref,CellSet &csref,TranscriptSet &tsref) : panel(panelref), cs(csref), ts(tsref)
{
 pseudoimage_built=false;

 psima=nullptr;

 orig_sp[0]=ipar.osp[0];
 orig_sp[1]=ipar.osp[1];
 origwidth=ipar.w;
 origheight=ipar.h;
 pix_per_bin=ipar.ppbin;
 reduced_sp[0]=pix_per_bin*orig_sp[0];
 reduced_sp[1]=pix_per_bin*orig_sp[1];

 extract_by_fov=ipar.exfov;
 if (extract_by_fov)
 {
  fov=ipar.idfov;
  ulx=uly=lrx=lry=0.0;
 }
 else
 {
  fov="";
  ulx=ipar.fovwin[0];
  uly=ipar.fovwin[1];
  lrx=ipar.fovwin[2];
  lry=ipar.fovwin[3];
  if (ulx<0.0 || uly<0.0 || lrx<0.0 || lry<0.0)
   ERRPROG("Error in area of interest window: at least one of the values is not a possitive real number.\n";)
 }

 only_targets=ipar.onlytargets;

 qvth=ipar.qvthres;

 genmarks=ipar.gmarks;

 if (DEB)
 {
  cout << "Output image to be extracted from image with (heigth,width)=(" << origheight << " x " << origwidth << ") and pixel spacing (vert,hor)=(" << orig_sp[1] << " x " << orig_sp[0] << ") micrometers.\n";
  cout << "Output image will have each pixel equivalent to ";
  if (pix_per_bin==1)
   cout << "one original pixel.";
  else
   cout << "(" << pix_per_bin << " x " << pix_per_bin << ") original pixels.";
  cout << " Therefore, pixel spacing will be (" << reduced_sp[1] << " x " << reduced_sp[0] << ") micrometers.\n";

  cout << "Frame to analize: ";
  if (extract_by_fov)
   cout << fov << endl;
  else
   cout << "(" << ulx << "," << uly << ") to (" << lrx << "," << lry << ")\n";
 }

 vector<int> requestedcodes = (only_targets ? panel.GetSortedTargets() : panel.GetSortedGenes());

 int code,min_code=panel.GetMinCodeword(),max_code=panel.GetMaxCodeword();

 if ( ( (ipar.numcodes.size()==0 && ipar.targetnames.size()==0) ) ||
      ( (ipar.numcodes.size()!=0 && ipar.targetnames.size()!=0) ) )
  ERRPROG("SpImage constructor: exactly one of the lists of targets by numbers or by names must be not empty.\n";)

 if (ipar.targetnames.size()==0)
 {
  for (size_t i=0; i<ipar.numcodes.size(); i++)
  {
   code=ipar.numcodes[i];
   if (code<min_code || code>max_code)
    ERRPROG("Error in target code list: " << code << " is not in the code list of the gene panel.\n";)
   if (!binary_search(requestedcodes.begin(),requestedcodes.end(),code))
   {
    if (only_targets)
     ERRPROG("Error in target code list: " << code << " is a valid code, but not a target.\n";)
    else
     ERRPROG("Error in gene code list: " << code << " is a valid code, but not a gene.\n";)
   }
   ecodes.push_back(code);
  }
 }
 else
 {
  for (size_t i=0; i<ipar.targetnames.size(); i++)
  {
   code = panel.GetCodewordFromName(ipar.targetnames[i]);
   if (code<min_code || code>max_code)
    ERRPROG("Error in target gene list: " << ipar.targetnames[i] << " is not in the list of names of the gene panel.\n";)
   if (!binary_search(requestedcodes.begin(),requestedcodes.end(),code))
   {
    if (only_targets)
     ERRPROG("Error in target code list: " << code << " is a valid code, but not a target.\n";)
    else
     ERRPROG("Error in gene code list: " << code << " is a valid code, but not a gene.\n";)
   }
   ecodes.push_back(code);
  }
 }

 if (DEB)
 {
  cout << " Scan will be performed in the target list for these " << ecodes.size() << " targets:\n[";
  for (size_t i=0; i<ecodes.size()-1; i++)
   cout << ecodes[i] << ",";
  cout << ecodes[ecodes.size()-1] << "]\n";
 }
}

bool SpRData::FindOrder(int ecode,size_t &order)
{
 size_t i=0;
 while (i<ecodes.size() && ecodes[i]!=ecode)
  i++;
 if (i<ecodes.size())
 {
  order=i;
  return(true);
 }
 return(false);
}

void SpRData::Build(string histname)
{
 vector<unsigned long> selected;
 float x,y;
 float minx=numeric_limits<float>::max();
 float miny=numeric_limits<float>::max();
 float maxx=numeric_limits<float>::min();
 float maxy=numeric_limits<float>::min();
 
 size_t lowqv=0;
 size_t notgene=0;
 if (extract_by_fov)
 {
  for (size_t i=0; i<ts.GetNumTr(); i++)
   if (fov=="ALL" || ts.GetFOV(i)==fov)
   {
    if (panel.IsGene(ts.GetCW(i)))
    {
     if (ts.GetQv(i)>=qvth)
     {
      x=ts.GetX(i);
      y=ts.GetY(i);
      if (x<minx)
       minx=x;
      if (x>maxx)
       maxx=x;
      if (y<miny)
       miny=y;
      if (y>maxy)
       maxy=y;
      selected.push_back(i);
     }
     else
      lowqv++;
    }
    else
     notgene++;
   }
   // There is no else here (if not inside FOV, nothing happens)
 }
 else
 {
  for (size_t i=0; i<ts.GetNumTr(); i++)
  {
   x=ts.GetX(i);
   y=ts.GetY(i);
   if (IsInside(x,y))
   {
    if (panel.IsGene(ts.GetCW(i)))
    {
     if (ts.GetQv(i)>=qvth)
      selected.push_back(i);
     else
      lowqv++;
    }
    else
     notgene++;
   }
  }
 }
 
 if (selected.size()==0)
  ERRPROG("No events found inside the requested region. No image can be written.\n";)
 
 SpPseudoImage::MyPoint ulp,lrp;
 if (extract_by_fov)
 {
  if (fov=="ALL")
  {
   ulp[0]=ulp[1]=0;
   // OJO: Así, o al revés???
   lrp[0]=origwidth*orig_sp[0];
   lrp[1]=origheight*orig_sp[1];
  }
  else
  {
   ulp[0]=minx-0.005*(maxx-minx);
   ulp[1]=miny-0.005*(maxy-miny);
   lrp[0]=maxx+0.005*(maxx-minx);
   lrp[1]=maxy+0.005*(maxy-miny);
  }
 }
 else
 {
  ulp[0]=ulx;
  ulp[1]=uly;
  lrp[0]=lrx;
  lrp[1]=lry;
 }
 
 // OJO: repasar. ¿Cómo solucionar esto en el caso general?
 unsigned long origidx,origidy,w,h;
 float dx=lrp[0]-ulp[0];
 float dy=lrp[1]-ulp[1];

 if (extract_by_fov && fov=="ALL")
 {
  origidx=origidy=0;
  // Integer division, without +0.5
  w=origwidth/pix_per_bin;
  h=origheight/pix_per_bin;
 }
 else
 {
  // Quitamos el 0.5 para ser consistentes con la imagen reducida real
  origidx=(unsigned long)((ulp[0]/orig_sp[0]));
  origidy=(unsigned long)((ulp[1]/orig_sp[1]));
 
  // OJO: Aquí también. Pero hay que añadir 1 por los puntos en el límite
  // O no, ya ni lo se...
  w=(unsigned long)((dx/reduced_sp[0]))+1;
  h=(unsigned long)((dy/reduced_sp[1]))+1;
 }

 if (DEB)
 {
  cout << selected.size() << " transcription events of genes with sufficient quality inside requested region.\n";
  if (lowqv>0)
   cout << lowqv << " transcripts discarded since their qv was lower than " << qvth << ".\n";
  if (notgene>0)
   cout << notgene << " transcripts discarded since they did not come from genes (probes, or others).\n";
  cout << "Region goes from (" << ulp[0] << "," << ulp[1] << ") to (" << lrp[0] << "," << lrp[1] << ")\n";
  cout << "Region size is (" << dx << " x " << dy << ") micrometers, (" << w << " x " << h << ") " << ((pix_per_bin>1) ? "macro" : "") << "pixels.\n";
 }

 size_t nsel=ecodes.size();
 float xev,yev;

 psima = new SpPseudoImage(w,h,reduced_sp,ulp,nsel);

 if (DEB)
 {
  cout << "Output image with (startrow,startcol)-(height,width,depth)=(" << origidx << "," << origidy << ")-(" << h << ", " << w << ", " << nsel << ")\n";
  cout << "Events not included in the selected codes will be ignored.\n";
  cout.flush();
 }

 size_t ord;
 size_t sumnsel=0;
 SpPseudoImage::My3DPoint p;
 SpPseudoImage::My3DIndex index;
 size_t reploc=0,replocplace=0;
 float dnuc;
 //float arearatio;
 //int attclust;
 for (size_t t=0; t<selected.size(); t++)
 {
  if (FindOrder(ts.GetCW(selected[t]),ord))
  {
   xev=ts.GetX(selected[t]);
   yev=ts.GetY(selected[t]);
   if (genmarks && ts.InsideCell(selected[t]))
   {
    dnuc=ts.GetDNuc(selected[t]);
    //arearatio=cs.GetAreaRatio(ts.GetCellId(selected[t]));
    //attclust=cs.GetAttClust(ts.GetCellId(selected[t]));
    //cout << selected[t] << " " << ts.GetCellId(selected[t]) << " " << arearatio << endl;
   }
   else
    dnuc=-1.0;
   /*
   {
    dnuc=arearatio=-1.0;
    attclust=UnknownCluster;
   }
   */
   p[0]=xev;
   p[1]=yev;
   p[2]=0.0;
   if (!psima->PPhysicalPointToIndex(p,index))
    ERRPROG("Point (" << xev << "," << yev << ") seems to be outside requested region. Indices are (" << index[0] << "," << index[1] <<") and origin index is at (" << origidx << "," << origidy << ")\n";)
   index[2]=ord;
   if (pix_per_bin>1)
    psima->SimpleInsertOrIncrementPixelValue(index);
   else
    psima->InsertOrIncrementPixelValue(index,dnuc,cs.GetCellRef(ts.GetCellId(selected[t])),reploc,replocplace,xev,yev);
  }
  else   // This means the event is not in the set of targets
   sumnsel++;
 }

 if (DEB)
 {
  if (pix_per_bin==1)
  {
   if (reploc>0)
    WARNPROG(reploc << " precise locations (up to " << num_dec_places << " decimals) have more than one event of expression (even they come from differente genes).";)
   if (replocplace>0)
    WARNPROG(reploc << " precise locations (up to " << num_dec_places << " decimals) have more than one event of expression FROM THE SAME GENE.\nThis is anomalous. Such events have not been stored.\n";)
  }
  cout << "Checking...";
  cout.flush();
 }

 if (DEB || histname != "none")
 {
  size_t sums[nsel],vmax[nsel],tsums=0;
  for (size_t k=0; k<nsel; k++)
  {
    sums[k]=0;
    vmax[k]=numeric_limits<size_t>::min();
  }
  size_t numv=0;
  psima->GetCounts(sums,tsums,vmax,numv);

  if (histname != "none")
  {
   ofstream f(histname.c_str());
   if (!f.is_open())
    WARNPROG("SpRData::Build: cannot open file " << histname << " to write.\n";)
   else
   {
    if (DEB)
       cout << "Saving histogram of genes present in file " << histname << endl;
    for (size_t k=0; k<ecodes.size(); k++)
     f << panel.GetNameFromCodeword(ecodes[k]) << "\t" << panel.GetIdFromCodeword(ecodes[k]) << "\t" << sums[k] << endl;
    f.close();
   }
  }

  if (DEB)
  {
   cout << " Found a total of " << tsums+sumnsel << " events along the " << numv << " voxels,";
   cout << "From the " << tsums+sumnsel << " events, " << tsums << " correspond to selected targets and " << sumnsel << " to others.\n";

   //cout << "Distribution by targets:\n";
   //cout << "Num_target  Target_code  Totals_for_target  Max_value_at_a_point\n";
   size_t minacc=sums[0],maxacc=0,minp=vmax[0],maxp=0;
   int code_of_min=0,code_of_max=0;
   for (size_t k=0; k<nsel; k++)
   {
    //cout << k << " " << ecodes[k] << " " << sums[k] << " " << vmax[k] << endl;
    if (sums[k]<minacc)
    {
     minacc=sums[k];
     code_of_min=ecodes[k];
    }
    if (sums[k]>maxacc)
    {
     maxacc=sums[k];
     code_of_max=ecodes[k];
    }
    if (vmax[k]>maxp)
     maxp=vmax[k];
    if (vmax[k]<minp)
     minp=vmax[k];
   }
   //if (accumulate_rest)
   // cout << nsel << " " << "REST" << " " << sums[nsel] << " " << vmax[nsel] << endl << endl;

   cout << "Target with identifier " << code_of_min << " has the smallest global expression with " << minacc << " events.\n";
   cout << "Target with identifier " << code_of_max << " has the largest global expression with " << maxacc << " events.\n";
   cout << "The smallest accumulated value (apart from 0) at a single point is " << minp << " events.\n";
   cout << "The largest accumulated value at a single point is " << maxp << " events.\n";
  }
 }
   
 pseudoimage_built=true;
}

void SpRData::ScriptBackgroundImage(string scr,string bgn)
{
 size_t x0=size_t(((pix_per_bin*ulx)/psima->GetXSpacing())+0.5);
 size_t y0=size_t(((pix_per_bin*uly)/psima->GetYSpacing())+0.5);
 size_t w=pix_per_bin*psima->GetWidth();
 size_t h=pix_per_bin*psima->GetHeight();

 ofstream f(scr);
 if (!f.is_open())
  ERRPROG("ScriptBackgroundImage: cannot open file " << scr << " to write.\n";)
 f << "import qupath.lib.gui.images.servers.*\n";
 f << "import qupath.lib.gui.viewer.overlays.*\n";;
 f << "def viewer=getCurrentViewer()\n";
 f << "def imageData=viewer.getImageData()\n";
 f << "def region = RegionRequest.createInstance(getCurrentServer().getPath()," << pix_per_bin << "," << x0 << "," << y0 << "," << w << "," << h << ")\n\n";
 f << "def newserver = new RenderedImageServer.Builder(imageData)\n";
 f << "                                       .downsamples(1)\n";
 f << "                                       .layers(new HierarchyOverlay(viewer.getImageRegionStore(), viewer.getOverlayOptions(), imageData))\n";
 f << "                                       .build()\n\n";
 f << "writeImageRegion(newserver, region,\"" << bgn << "\")\n";
 f.close();
}

void SpRData::SaveAsImages(string outname,string idescname)
{
 if (!pseudoimage_built)
  ERRPROG("SaveAsImages: Error in save: cannot save an image that has not been built. Call function BuildImage() first.\n";)
 if (psima==nullptr)
  ERRPROG("SaveAsImages: Error in save: image is supposed to have been built but pseudoimage pointer is null.\n";)

 string fusefilename=idescname.substr(0,idescname.find_last_of("."))+".fuse";
 ofstream ff(fusefilename.c_str());
 if (!ff.is_open())
  ERRPROG("SaveAsImages: cannot open fuse file " << fusefilename << " to write.\n";)

 int nc=psima->GetCoderest();
 vector<string> outnames;
 string basename=outname.substr(0,outname.find_last_of("."));
 string background_iname=basename+"backg.tif";
 ff << background_iname << endl;
 char num[16],name[64];
 RGB_color col;
 for (int i=0; i<nc; i++)
 {
  strcpy(name,panel.GetNameFromCodeword(ecodes[i]).c_str());
  sprintf(num,"%03d_%s",ecodes[i],name);
  outnames.push_back(basename+string(num)+".tif");
  col=panel.GetColorFromCodeword(ecodes[i]);
  ff << basename+string(num)+".tif" << " " << int(col.r) << " " << int(col.g) << " " << int(col.b) << endl;
 }
 ff.close();

 ImageBinType::Pointer imres = ImageBinType::New();
 ImageBinType::SizeType s;

 s[0]=psima->GetWidth();
 s[1]=psima->GetHeight();

 ImageBinType::RegionType r(s);
 imres->SetRegions(r);
 ImageSpType ressp;

 ressp[0]=psima->GetXSpacing();
 ressp[1]=psima->GetYSpacing();
 imres->SetSpacing(ressp);

 imres->Allocate();
 SpPseudoImage::PixContent lp;
 SpPseudoImage::MyPoint p;
 ImageBinType::IndexType idx;
 size_t numout=0;

 for (size_t outs=0; outs<outnames.size(); outs++)
 {
  if (DEB)
   cout << "Generating image " << outnames[outs] << " corresponding to gene code " << ecodes[outs] << endl;
  imres->FillBuffer(0);
  size_t nmarkpoints=0;
  for (size_t r=0; r<s[1]; r++)
   for (size_t c=0; c<s[0]; c++)
   {
    if (psima->PIndexToPhysicalPoint(r,c,p))
    {
     lp=psima->GetPixList(r,c);

     if (lp.size()>0)
      for (size_t pos=0; pos<lp.size(); pos++)
       if (size_t(lp[pos].code)==outs)
       {
        idx[0]=c;
        idx[1]=s[1]-r-1;
        imres->SetPixel(idx,255);
        nmarkpoints++;
      }
    }
    else
     numout++;
   }

  if (DEB)
  {
   cout << "   " << nmarkpoints << " points marked.\n";
   if (numout>0)
    WARNPROG("   In this image "<< numout << " points are out of image\n";)
   else
    cout << "  No points out of image. Good.\n";
  }

  if (outnames[outs].find(".tif")==string::npos && outnames[outs].find(".TIF")==string::npos)
   ERRPROG("Output image format must be tiff.";)

  ImageBinWriterType::Pointer imw=ImageBinWriterType::New();
  imw->SetFileName(outnames[outs]);

  TIFFIOType::Pointer imtif=TIFFIOType::New();
  imw->SetImageIO(imtif);
  imw->SetInput(imres);
  try
  {
   imw->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
   cerr << e.what() << endl;
   ERRPROG("SaveAsImages: cannot write file " << outname[outs] << ".\n";)
  }
 }
 if (DEB)
  cout << "Generating and saving groovy script to extract background image " << background_iname << endl;

 string script_name=idescname.substr(0,idescname.find_last_of("."))+".groovy";
 ScriptBackgroundImage(script_name,background_iname);
}

void SpRData::Save(string outname,char sep)
{
 if (!pseudoimage_built)
  ERRPROG("Save: Error in save: cannot save an image that has not been built. Call function BuildImage() first.\n";)
 if (psima==nullptr)
  ERRPROG("Save: Error in save: image is supposed to have been built but pseudoimage pointer is null.\n";)

 int nc=psima->GetCoderest();
 vector<string> outnames;
 string basename=outname.substr(0,outname.find_last_of("."));
 char num[16];
 for (int i=0; i<nc; i++)
 {
  sprintf(num,"%03d",i);
  outnames.push_back(basename+string(num)+".csv");
 }

 ofstream f[nc];
 for (int i=0; i<nc; i++)
 {
   f[i].open(outnames[i].c_str());
   if (!f[i].is_open())
    ERRPROG("Error in save: cannot open file " << outnames[i] << " to write in it.\n";)
   if (genmarks)
    f[i] << "x" << sep << "y" << sep << "NucDist" << sep << "AreaRatio" << sep << "ClustNum\n";
   else
    f[i] << "x" << sep << "y\n";

 }

 SpPseudoImage::PixContent lp;
 SpPseudoImage::MyPoint p;
 //MyValPixelType *marks = new MyValPixelType[nc];
 size_t nev=0;

 size_t nlocs=0;
 char sp[128];
 for (int i=0;i<128; i++)
  sp[i]=0;
 char fmt[16]="%-P.Qf%c%-P.Qf";
 fmt[2]=fmt[10]='0'+num_int_places;
 fmt[4]=fmt[12]='0'+num_dec_places;
 size_t ndups=0;

 float dn,ar;
 int attcl;
 float xmin,ymin;
 xmin=ymin=numeric_limits<float>::max();
 float xmax,ymax;
 xmax=ymax=numeric_limits<float>::min();

 for (size_t r=0; r<psima->GetHeight(); r++)
  for (size_t c=0; c<psima->GetWidth(); c++)
  {
   if (psima->PIndexToPhysicalPoint(r,c,p))
   {
    lp=psima->GetPixList(r,c);

    if (lp.size()>0)
    {
     nlocs++;
     if (pix_per_bin>1)
     {
      //for (int i=0; i<nc; i++)
      // marks[i]=0;
      sprintf(sp,fmt,p[0],sep,p[1]);
      if (p[0]<xmin)
       xmin=p[1];
      if (p[0]>xmax)
       xmax=p[0];
      if (p[1]<ymin)
       ymin=p[1];
      if (p[1]>ymax)
       ymax=p[1];
      for (size_t pos=0; pos<lp.size(); pos++)
      {
      // marks[lp[pos].code]=lp[pos].v;
       nev+=lp[pos].v;
       f[lp[pos].code] << sp << sep << lp[pos].v << endl;
      }
      //for (int i=0; i<nc-1; i++)
      // f[lp[0].code] << marks[i] << "\t";
      //f[lp[0].code] << marks[nc-1] << endl;
     }
     else
     {
      for (size_t pos=0; pos<lp.size(); pos++)
      {
       if (lp[pos].Realcoord.size()==2 &&
           lp[pos].Realcoord[0]==lp[pos].Realcoord[1])
       {
        WARNPROG("Two events are at the same location, at least at the given precision, " << num_dec_places << " decimals. We'll write only one.\nIf this is a problem, you might consider change the constant num_dec_places at file ndplaces.h and recompile.\nThe common location is (" << lp[pos].Realcoord[0].first << "," << lp[pos].Realcoord[0].second << ")\n";)
        //for (int i=0; i<nc; i++)
        //  marks[i]=0;
        // marks[lp[pos].code]=1;
        sprintf(sp,fmt,lp[pos].Realcoord[0].first,lp[pos].Realcoord[0].second);
        if (lp[pos].Realcoord[0].first<xmin)
         xmin=lp[pos].Realcoord[0].first;
        if (lp[pos].Realcoord[0].first>xmax)
         xmax=lp[pos].Realcoord[0].first;
        if (lp[pos].Realcoord[0].second<ymin)
         ymin=lp[pos].Realcoord[0].second;
        if (lp[pos].Realcoord[0].second>ymax)
         ymax=lp[pos].Realcoord[0].second;
        if (genmarks)
        {
         dn = lp[pos].DistNuc[0];
         //ar = lp[pos].AreaRatios[0];
         //attcl = lp[pos].AttClust[0];
         Cell &c=lp[pos].vcells[0];
         ar = c.GetAreaRatio();
         attcl = c.GetCluster();
         if (dn<0)
          f[lp[pos].code] << sp << sep << "NA" << sep << "NA" << sep << attcl << endl;
         else
         {
          f[lp[pos].code] << sp << sep << dn << sep;
          if (isnan(ar))
           f[lp[pos].code] << "NA" << sep << attcl << endl;
          else
           f[lp[pos].code] << ar << sep << attcl << endl;
         }
        }
        else
         f[lp[pos].code] << sp << endl;
        //for (int i=0; i<nc-1; i++)
        // f << marks[i] << "\t";
        //f << marks[nc-1] << endl;
        nev++;
        ndups++;
       }
       else
       {
        if (lp[pos].Realcoord.size()>2)
        {
         WARNPROG("There are " << lp[pos].Realcoord.size() << " (v=" << lp[pos].v << ") events of gene " << lp[pos].code << " at image location (" << r << "," << c << ").\nThis is strange. No event will be written.";)
         ndups++;
        }
        else
        {
         //for (int i=0; i<nc; i++)
         // marks[i]=0;
         //marks[lp[pos].code]=1;
         for (size_t kp=0; kp<lp[pos].Realcoord.size(); kp++)
         {
          sprintf(sp,fmt,lp[pos].Realcoord[kp].first,sep,lp[pos].Realcoord[kp].second);
          if (lp[pos].Realcoord[kp].first<xmin)
           xmin=lp[pos].Realcoord[kp].first;
          if (lp[pos].Realcoord[kp].first>xmax)
           xmax=lp[pos].Realcoord[kp].first;
          if (lp[pos].Realcoord[kp].second<ymin)
           ymin=lp[pos].Realcoord[kp].second;
          if (lp[pos].Realcoord[kp].second>ymax)
           ymax=lp[pos].Realcoord[kp].second;
          if (genmarks)
          {
           dn = lp[pos].DistNuc[kp];
           //ar = lp[pos].AreaRatios[0];
           //attcl = lp[pos].AttClust[0];
           Cell &c=lp[pos].vcells[0];
           ar = c.GetAreaRatio();
           attcl = c.GetCluster();
           if (dn<0)
            f[lp[pos].code] << sp << sep << "NA" << sep << "NA" << sep << attcl << endl;
           else
           {
            f[lp[pos].code] << sp << sep << dn << sep;
            if (isnan(ar))
             f[lp[pos].code] << "NA" << sep << attcl << endl;
            else
             f[lp[pos].code] << ar << sep << attcl << endl;
           }
          }
          else
           f[lp[pos].code] << sp << endl;

          //for (int i=0; i<nc-1; i++)
          // f << marks[i] << "\t";
          //f << marks[nc-1] << endl;
          nev++;
         }
        }
       }
      }
     }
    }
   }
  }
 //delete[] marks;
 if (DEB)
 {
  cout << nev << " events saved at " << nlocs << " locations. " << 100.0*float(nev-nlocs)/float(nlocs) << "% of the locations shared more than 1 expression event\n";
  if (pix_per_bin==1 && ndups>0)
  cout << ndups << " events discarded since they had common location and gene expression.\n";
 }
 for (int i=0; i<nc; i++)
  f[i].close();

 string limname=basename+"limits.R";
 float dx=xmax-xmin;
 float dy=ymax-ymin;
 xmin -= 0.1*dx;
 if (xmin<0.0)
  xmin=0.0;
 ymin -= 0.1*dy;
 if (ymin<0.0)
  ymin=0.0;
 xmax += 0.1*dx;
 ymax += 0.1*dy;
 ofstream g(limname);
 if (!g.is_open())
  ERRPROG("Cannot open file " << limname << " to write limits in it.\n";)
 g << "xlimits <- c(" << xmin << "," << xmax << ")\n";
 g << "ylimits <- c(" << ymin << "," << ymax << ")\n";
 g.close();

 string genetabname=basename+"genes.csv";
 ofstream h(genetabname);
 if (!h.is_open())
  ERRPROG("Cannot open file " << genetabname << " to write gene table in it.\n";)
 h << "GeneNum Codeword Id\n";
 for (size_t i=0; i<ecodes.size(); i++)
  h << i << " " << ecodes[i] << " " << panel.GetIdFromCodeword(ecodes[i]) << endl;
 h.close();
}

