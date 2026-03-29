#include "csv2bin.h"

using namespace std;

bool DEB=true;
string progname;

void Usage()
{
 cerr << "Usage:\n" << progname << " <object_type> <csv_file> [<csv_cell_boundary_file>] [<csv_nucleus_boundary_file>] <bin_file>\n";
 cerr << " where object_type is:\n";
 cerr << "  'p' to convert a gene panel,\n";
 cerr << "  'c' to convert a cell set,\n";
 cerr << "  'b' to convert a cell set plus a cell boundary,\n";
 cerr << "  'n' to convert a cell set plus a nucleus boundary,\n";
 cerr << "  'a' to convert a cell set, a cell boundary and a nucleus boundary, or\n";
 cerr << "  't' to convert a transcript set.\n";
 cerr << " Arguments between [ ] will be required for the appropriate types (b,n or a) according to the type.\n";
 cerr << " The resulting binary file must be always the last argument.\n";
 exit(1);
}

int main(int argc,char *argv[])
{
 progname=string(argv[0]);

 if (argc<4)
  Usage();

 string ncsv,nbcsv,nncsv,nbin;

 char otype=tolower(argv[1][0]);

 if (otype=='p' || otype=='c' || otype=='t')
 {
  if (argc!=4)
   Usage();
  ncsv=string(argv[2]);
  nbin=string(argv[3]);
 }

 if (otype=='b')
 {
  if (argc!=5)
   Usage();
  ncsv=string(argv[2]);
  nbcsv=string(argv[3]);
  nbin=string(argv[4]);
 }

 if (otype=='n')
 {
  if (argc!=5)
   Usage();
  ncsv=string(argv[2]);
  nncsv=string(argv[3]);
  nbin=string(argv[4]);
 }

 if (otype=='a')
 {
  if (argc!=6)
   Usage();
  ncsv=string(argv[2]);
  nbcsv=string(argv[3]);
  nncsv=string(argv[4]);
  nbin=string(argv[5]);
 }

 switch (otype)
 {
     case 'p':
     {
      GenePanel p(ncsv);
      p.SaveAsBinary(nbin);
     }
     break;
     case 'c':
     {
      CellSet cs(ncsv);
      cs.SaveAsBinary(nbin);
     }
     break;
     case 't':
     {
      TranscriptSet tset(ncsv);
      tset.SaveAsBinary(nbin);
     }
     break;
     case 'b':
     {
      CellSet cs(ncsv);
      cs.AddBoundaryInfo(nbcsv,true);
      cs.SaveAsBinary(nbin);
     }
     break;
     case 'n':
     {
      CellSet cs(ncsv);
      cs.AddBoundaryInfo(nncsv,false);
      cs.SaveAsBinary(nbin);
     }
     break;
     case 'a':
     {
      CellSet cs(ncsv);
      cs.AddBoundaryInfo(nbcsv,true);
      cs.AddBoundaryInfo(nncsv,false);
      cs.SaveAsBinary(nbin);
     }
     break;
     default:
         cerr << "Unknown object type.\n";
         exit(1);
         break;
 }

 return(0);
}
