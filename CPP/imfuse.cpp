#include "imfuse.h"

using namespace std;
using namespace itk;

string progname="imfuse";
bool DEB=true;

void Usage()
{
 cerr << "Usage:\n   " << progname << " <imagelist_fusion_file> <fusion_image> [-adjcup|-adjcdown]\n";
 exit(1);
}

ImageBinType::Pointer AdjustContrast(ImageColType::Pointer ima,AdjustModes adjc)
{
 ImageBinType::Pointer grayima=ImageBinType::New();
 grayima->SetRegions(ima->GetLargestPossibleRegion());
 grayima->SetSpacing(ima->GetSpacing());
 grayima->Allocate();

 ImageRegionIteratorWithIndex<ImageColType> inpit(ima,ima->GetLargestPossibleRegion());
 ImageRegionIteratorWithIndex<ImageBinType> outit(grayima,grayima->GetLargestPossibleRegion());
 inpit.GoToBegin();
 outit.GoToBegin();
 ColPixelType c;
 BinPixelType b;
 while (!inpit.IsAtEnd())
 {
  c=inpit.Get();
  b=(c[0]+c[1]+c[2])/3;
  outit.Set(b);
  ++inpit;
  ++outit;
 }

 ImageCalculatorFilterType::Pointer calc = ImageCalculatorFilterType::New();
 calc->SetImage(grayima);
 calc->Compute();
 BinPixelType minv=calc->GetMinimum();
 BinPixelType maxv=calc->GetMaximum();

 if (maxv==minv)
  ERRPROG("AdjustContrast: the extracted image has all the same gray value (" << int(minv) << ")\n";)

 if (DEB)
  cout << "Background image with values in [" << int(minv) << "," << int(maxv) << "]\n";

 if (adjc!=NoAdj)
 {
  BinPixelType abmax=NumericTraits<BinPixelType>::max();
  double adjust_gamma = ((adjc == Up) ? adjust_gamma_up : adjust_gamma_down );
  double gamma;
  if (maxv==abmax)
   gamma=log(adjust_gamma)/log(0.9);
  else
   gamma=log(adjust_gamma)/log(double(maxv)/double(abmax));

  if (DEB)
  {
   double nv=double(maxv)/double(abmax);
   double tv=exp(gamma*log(nv));
   BinPixelType rv=BinPixelType(tv*double(abmax));
   cout << "gamma value of " << gamma << " will transform value " << int(maxv) << " (normalized: " << nv << ") into value " << int(rv) << " (normalized: " << tv << ").\n";
  }

  GammaFilterType::Pointer gfilter = GammaFilterType::New();
  gfilter->SetInput(grayima);
  gfilter->GetFunctor().SetGammaValue(gamma);
  gfilter->Update();
  ImageBinType::Pointer adjusted=gfilter->GetOutput();

  // It was saved just as a test...
  /*
  TIFFIOType::Pointer tf=TIFFIOType::New();
  ImageBinWriterType::Pointer wrb=ImageBinWriterType::New();
  wrb->SetImageIO(tf);
  wrb->SetFileName("adjusted_uc.tif");
  wrb->SetInput(adjusted);
  wrb->Update();
  */

  return(adjusted);
 }

 // If we are here, we don't want to adjust the contrast
 return(grayima);
}

ImageColType::Pointer CreateImageLike(pair<string,ColPixelType> ct,AdjustModes adjc,string backname)
{
 {
  ifstream f(backname.c_str());
  if (!f.is_open())
   ERRPROG("CreateImageLike: background image " << backname << " cannot be read.\n";)
  f.close();
 }
 // First, read the background image that will dictate the obligued dimensions for the others

 ImageColReaderType::Pointer imcol_reader = ImageColReaderType::New();
 imcol_reader->SetFileName(backname);
 TIFFIOType::Pointer imtif = TIFFIOType::New();
 imcol_reader->SetImageIO(imtif);
 try
 {
  imcol_reader->Update();
 }
 catch (const itk::ExceptionObject & e)
 {
  cerr << e.what() << endl;
  ERRPROG("CreateImageLike: cannot read file " << backname << ". Is it a valid TIFF image?\n";)
 }

 ImageColType::Pointer imcol_pre = imcol_reader->GetOutput();

 ImageBinType::Pointer imgray_uc = AdjustContrast(imcol_pre,adjc);

 ImageRegionType r=imgray_uc->GetLargestPossibleRegion();
 ImageSizeType s=r.GetSize();
 ImageSpType sp=imgray_uc->GetSpacing();

 // Now, read the first image and check if it has the correct dimensions
 ImageBinReaderType::Pointer imbin_reader = ImageBinReaderType::New();
 imbin_reader->SetFileName(ct.first);
 imbin_reader->Update();
 ImageBinType::Pointer imbin = imbin_reader->GetOutput();

 ImageRegionType rbin = imbin->GetLargestPossibleRegion();
 ImageSizeType sbin = rbin.GetSize();
 ImageSpType spbin = imbin->GetSpacing();

 if (s[0]!=sbin[0] || s[1]!=sbin[1])
  ERRPROG("CreateImageLike: background image " << backname << " has different size than first image, " << ct.first << ".\nBackground image is " << s[0] << " x " << s[1] << " whereas first image is " << sbin[0] << " x " << sbin[1] << ")\n";)

 if (sp[0]!=spbin[0] || sp[1]!=spbin[1])
  WARNPROG("CreateImageLike: background image " << backname << " has different spacing than first image, " << ct.first << ".\nBackground image spacing is " << sp[0] << " x " << sp[1] << " whereas first image is " << spbin[0] << " x " << spbin[1] << ")\n";)

 // Now, create the empty color image with the same dimensions and spacing than the background image
 ImageColType::Pointer imcol = ImageColType::New();
 ImageRegionType::SizeType scol;
 scol[0] = s[0];
 scol[1] = s[1];
 ImageRegionType rcol(scol);
 imcol->SetRegions(rcol);

 ImageSpType spcol;
 spcol[0] = sp[0];
 spcol[1] = sp[1];
 imcol->SetSpacing(spcol);

 imcol->Allocate();

 if (DEB)
 {
  cout << "Copying background from image " << backname << " converting to color with all components equal...\n";
 }
 ImageRegionIteratorWithIndex<ImageBinType> itback(imgray_uc,r);
 itback.GoToBegin();
 ImageIndexType idx;
 ColPixelType v;
 BinPixelType vg;
 while (!itback.IsAtEnd())
 {
  idx=itback.GetIndex();
  vg=imgray_uc->GetPixel(idx);
  v[0]=v[1]=v[2]=vg;
  imcol->SetPixel(idx,v);
  ++itback;
 }
 if (DEB)
  cout << "Done.!\n";

 if (DEB)
 {
  cout << "Running through image " << ct.first << endl;
  cout.flush();
 }

 ImageRegionIteratorWithIndex<ImageBinType> itfirst(imbin,rbin);
 itfirst.GoToBegin();
 size_t np=0;
 ImageIndexType idxn;
 while (!itfirst.IsAtEnd())
 {
  idx=itfirst.GetIndex();
  if (itfirst.Get()!=0)
  {
   unsigned inf0=(idx[0]>=neisize) ? (idx[0]-neisize) : 0;
   unsigned inf1=(idx[1]>=neisize) ? (idx[1]-neisize) : 0;
   unsigned sup0=(idx[0]+neisize<long(sbin[0])) ? idx[0]+neisize : sbin[0]-1;
   unsigned sup1=(idx[1]+neisize<long(sbin[1])) ? idx[1]+neisize : sbin[1]-1;
   for (idxn[0]=inf0; idxn[0]<=sup0; idxn[0]++)
    for (idxn[1]=inf1; idxn[1]<=sup1; idxn[1]++)
     imcol->SetPixel(idxn,ct.second);
   np++;
  }
  ++itfirst;
 }

 if (DEB)
 {
  cout << "Done. " << np << " pixels set to color [" << int(ct.second[0]) << "," << int(ct.second[1]) << "," << int(ct.second[2]) << "].\n";
  cout.flush();
 }

 return(imcol);
}

void OverlayOnImage(ImageColType::Pointer ima,pair<string,ColPixelType> ct)
{
 ImageBinReaderType::Pointer imbin_reader = ImageBinReaderType::New();
 imbin_reader->SetFileName(ct.first);
 imbin_reader->Update();

 ImageBinType::Pointer imbin = imbin_reader->GetOutput();
 ImageRegionType rbin = imbin->GetLargestPossibleRegion();
 ImageSizeType sbin = rbin.GetSize();
 ImageSpType spbin = imbin->GetSpacing();

 ImageRegionType rcol = ima->GetLargestPossibleRegion();
 ImageSizeType scol = rcol.GetSize();
 ImageSpType spcol = ima->GetSpacing();

 if (sbin[0]!=scol[0] || sbin[1]!=scol[1])
  ERRPROG("OverlayOnImage: image " << ct.first << " has not the correct size. It is " << sbin[0] << " x " << sbin[1] << " instead of " << scol[0] << " x " << scol[1] << ".\n";)

 if (spbin[0]!=spcol[0] || spbin[1]!=spcol[1])
  WARNPROG("OverlayOnImage: image " << ct.first << " has not the correct spacing. It is (" << spbin[0] << "," << spbin[1] << ") instead of (" << spcol[0] << "," << spcol[1] << ").\n";)

 if (DEB)
  cout << "Read image " << ct.first << " which has at least the correct size.\n";

 if (DEB)
 {
  cout << "Running through image " << ct.first << endl;
  cout.flush();
 }
 ImageRegionIteratorWithIndex<ImageBinType> itcurrent(imbin,rbin);
 ImageIndexType idx;
 itcurrent.GoToBegin();
 ColPixelType vc;
 size_t num_conflicts=0;
 size_t np=0;
 ImageIndexType idxn;
 while (!itcurrent.IsAtEnd())
 {
  idx=itcurrent.GetIndex();
  if (itcurrent.Get()!=0)
  {
   vc = ima->GetPixel(idx);
   if (! ( vc[0]==vc[1] && vc[1]==vc[2] ) )
    num_conflicts++;
   unsigned inf0=(idx[0]>=neisize) ? (idx[0]-neisize) : 0;
   unsigned inf1=(idx[1]>=neisize) ? (idx[1]-neisize) : 0;
   unsigned sup0=(idx[0]+neisize<long(sbin[0])) ? idx[0]+neisize : sbin[0]-1;
   unsigned sup1=(idx[1]+neisize<long(sbin[1])) ? idx[1]+neisize : sbin[1]-1;
   for (idxn[0]=inf0; idxn[0]<=sup0; idxn[0]++)
    for (idxn[1]=inf1; idxn[1]<=sup1; idxn[1]++)
     ima->SetPixel(idxn,ct.second);
   np++;
  }
  ++itcurrent;
 }
 if (DEB)
 {
  cout << "Done. " << np << " pixels set to color [" << int(ct.second[0]) << "," << int(ct.second[1]) << "," << int(ct.second[2]) << "].\n";
  cout.flush();
 }
 if (DEB && num_conflicts>0)
  cout << "Found " << num_conflicts << " conflicts when overlaying image " << ct.first << ".\n";
}

int main(int argc,char *argv[])
{
 progname=string(argv[0]);

 if ((argc!=3) && (argc!=4))
  Usage();

 string fname=string(argv[1]);
 string resname=string(argv[2]);

 AdjustModes adjust=NoAdj;
 if (argc==4)
 {
  if (string(argv[3])=="-adjcup")
   adjust=Up;
  else
  {
   if (string(argv[3])=="-adjcdown")
    adjust=Down;
   else
    Usage();
  }
 }

 ifstream f(fname.c_str());
 if (!f.is_open())
  ERRPROG("main: cannot open file " << fname << " to read.\n";)

 // First, read the background image name
 string backname;
 f >> backname;

 // Now, read each line with the name of one overlay image and the RGB color to overlay it.
 vector<pair<string,ColPixelType>> coltab;
 string cname;
 int r,g,b;
 int line=1;
 f >> cname >> r >> g >> b;
 do
 {
  pair<string,ColPixelType> ct;
  ct.first=cname;
  if (r<0 || r>255 || g<0 || g>255 || b<0 || b>255)
   ERRPROG("main: incorrect color number in line " << line << " of file " << fname << ". It is not an integer number in [0..255]\n";)
  ct.second[0]=(unsigned char)r;
  ct.second[1]=(unsigned char)g;
  ct.second[2]=(unsigned char)b;
  coltab.push_back(ct);
  line++;
  if (!f.eof())
   f >> cname >> r >> g >> b;
 }
 while (!f.eof());
 f.close();

 // A color image with the same dimensions and spacing of the gray background image is created and filled with the gray image, possibly gamma-transformed
 ImageColType::Pointer ima=CreateImageLike(coltab[0],adjust,backname);
 // ... and each gene with its corresponding color is overlaid on it.

 for (size_t i=1; i<coltab.size(); i++)
  OverlayOnImage(ima,coltab[i]);

 if (DEB)
  cout << "\nWriting image " << resname << "...";

 ImageColWriterType::Pointer imw = ImageColWriterType::New();
 imw->SetFileName(resname);
 TIFFIOType::Pointer imtif = TIFFIOType::New();
 imw->SetImageIO(imtif);
 imw->SetInput(ima);
 imw->Update();

 if (DEB)
  cout << "Done!\n";

 return(0);
}

