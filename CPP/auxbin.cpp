#include "auxbin.h"
// @cond DO_NOT_DOCUMENT

using namespace std;

InputBufferedFile::InputBufferedFile(std::string fname,bool debug) : BufferedFile(debug)
{
 filename=fname;
 mem=nullptr;
 memsize=0;
 off=0;

 f.open(filename.c_str(),ios::binary);
 if (!f.is_open())
 {
  cerr << "Low level error from InputBufferedFile constructor: file " << fname << " cannot be opened to read. Check existence and permissions.\n";
  exit(1);
 }
 // Get the file size:
 f.seekg(0,f.end);
 memsize=f.tellg();
 // Book a chunk of memory to read it
 mem = new (nothrow) char[memsize];
 if (mem==nullptr)
 {
  cerr << "Low level error from InputBufferedFile constructor: Cannot book " << memsize << " bytes of memory to read file " << fname << ".\n";
  f.close();
  exit(1);
 }
 // If everything has gone well, go the beginning of the file and read it to memory as a whole
 f.seekg(0,f.beg);
 f.read(mem,memsize);
 if (!f)
 {
  cerr << "Low level error from InputBufferedFile constructor: only " << f.gcount() << " of the " << memsize << " bytes of the file " << fname << " could be read to memory.\n";
  f.close();
  if (mem!=nullptr)
   delete[] mem;
  exit(1);
 }

 if (deb)
  cout << "Opening ifstream associated to file " << filename << ". InputBufferedFile constructed with memory chunk of " << memsize << " bytes.\n";
}

InputBufferedFile::InputBufferedFile(std::string fname) : InputBufferedFile(fname,false)
{}

InputBufferedFile::~InputBufferedFile()
{
 if (f.is_open())
  f.close();
}

OutputBufferedFile::OutputBufferedFile(std::string fname, unsigned short chunk_size, bool debug) : BufferedFile(debug)
{
 filename=fname;
 mem=nullptr;
 memsize=0;
 off=0;

 f.open(filename.c_str(),ios::binary);
 if (!f.is_open())
 {
  cerr << "Low level error from OutputBufferedFile constructor: file " << fname << " cannot be opened to write. Check existence and permissions.\n";
  exit(1);
 }
 // Book a chunk of memory to store temporarily the data before physically writing them
 memsize=1024*1024*chunk_size;
 mem = new (nothrow) char[memsize];
 if (mem==nullptr)
 {
  cerr << "Low level error from OutputBufferedFile constructor: Cannot book " << memsize << " bytes of memory to write file " << fname << ".\n";
  if (f.is_open())
   f.close();
  exit(1);
 }

 if (deb)
  cout << "OutputBufferedFile constructed with memory chunk of " << memsize << " bytes associated to file " << fname << endl;
}

OutputBufferedFile::OutputBufferedFile(std::string fname) : OutputBufferedFile(fname,default_out_membuf,false)
{}

OutputBufferedFile::OutputBufferedFile(std::string fname,unsigned short chsize) : OutputBufferedFile(fname,chsize,false)
{}

OutputBufferedFile::OutputBufferedFile(std::string fname,bool debug) : OutputBufferedFile::OutputBufferedFile(fname, default_out_membuf, debug)
{}

OutputBufferedFile::~OutputBufferedFile()
{
 if (mem!=nullptr)
 {
  // Write to the file the part of the chunk still pending...
  f.write(mem,off);
  f.flush();
  if (deb)
   cout << "Last " << off << " bytes of OutputBufferedFile associated to file " << filename << " have been written. Closing ofstream.\n";
 }
 if (f.is_open())
  f.close();
}

BufferedFile::BufferedFile()
{
 deb=false;
}

BufferedFile::BufferedFile(bool debug)
{
 deb=debug;
}

BufferedFile::~BufferedFile()
{
 if (mem!=nullptr)
  delete[] mem;
}

// @endcond


