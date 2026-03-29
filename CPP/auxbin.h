// @cond DO_NOT_DOCUMENT
#ifndef _AUXBIN_H
#define _AUXBIN_H

#include <iostream>
#include <fstream>
#include <string.h> // To use memcpy
#include <vector>
#include <type_traits>
#include <typeinfo>

// This is the maximum length allowed for strings
const unsigned int max_string_len=254;
// This is the length of the temporary array used to hold the strings before being written.
// Its first byte will contain the string length, so this amount must be the previous one, plus 1
constexpr unsigned int tmp_max_string_len=max_string_len+1;

// The buffer size for output buffered streams. Let's say 128 Mb
constexpr size_t default_out_membuf=128;

class BufferedFile
{
 public:
  BufferedFile();
  BufferedFile(bool debug);
  ~BufferedFile();
 protected:
  std::string filename;
  char *mem;
  size_t memsize;
  size_t off;
  bool deb;
};

class InputBufferedFile: public BufferedFile
{
 public:
  InputBufferedFile(std::string fname);
  InputBufferedFile(std::string fname,bool debug);
  ~InputBufferedFile();

  // Special function to read strings, which are stored in memory as a first byte, indicating its length (from 0 to 254) followed by the string bytes. No ending \0, added when read.
  void inline readb(std::string &s)
  {
   unsigned char sl;
   if (off+1 < memsize)
   {
    sl=mem[off];
    off++;
   }
   else
   {
    std::cerr << "Low level error from InputBufferedFile::readb: trying to read beyond the end of the booked memory.\n";
    exit(1);
   }
   if (off+sl <= memsize)
   {
    char dummy[tmp_max_string_len];
    for (size_t i=0; i<sl; i++)
     dummy[i]=mem[off+i];
    off+=sl;
    dummy[sl]='\0';
    s=std::string(dummy);
    if (deb)
     std::cout << "file " << filename << ": " << sl+1 << " bytes assigned from buffer for a string. Current offset: " << off << std::endl;
   }
  };

  // Function to read any other fundamental data type
  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value, void>::type
  inline readb(T &v)
  {
   // If we are inside the booked memory...
   if (off+sizeof(T) <= memsize)
   {
    // ... copy as many bytes as required in the variable to be filled
    memcpy(reinterpret_cast<char*>(&v),mem+off,sizeof(T));
    // ... and advance the memory offset.
    off += sizeof(T);
    if (deb)
     std::cout << "file " << filename << ": " << sizeof(T) << " bytes assigned from buffer for a " << typeid(T).name() << ". Current offset: " << off << std::endl;
   }
   else
   {
    std::cerr << "Low level error from InputBufferedFile::readb: trying to read beyond the end of the booked memory.\n";
    exit(1);
   }
  };

  // Function to read vectors of any fundamental data type
  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value, void>::type
  inline readb(std::vector<T> &v)
  {
   // First, let's read the length of the vector (number of elements)
   size_t sv;
   if (off+sizeof(size_t) <= memsize)
   {
    memcpy(reinterpret_cast<char*>(&sv),mem+off,sizeof(size_t));
    v.resize(sv);
    // Once read, advance the memory offset
    off += sizeof(size_t);
    // Now, read one by one the elements of the vector
    if (off+sv*sizeof(T) <= memsize)
    {
     for (size_t i=0; i<sv; i++,off+=sizeof(T))
      memcpy(reinterpret_cast<char*>(&(v[i])),mem+off,sizeof(T));
     // Advance of the memory offset has been done by the loop update part
     if (deb)
      std::cout << "file " << filename << ": vector of " << sv << " " << typeid(T).name() << " so " << sv*sizeof(T)+sizeof(size_t) << " bytes assigned from buffer. Current offset: " << off << std::endl;
    }
    else
    {
     std::cerr << "Low level error from InputBufferedFile::readb: trying to read beyond the end of the booked memory.\n";
     exit(1);
    }
   }
   else
   {
    std::cerr << "Low level error from InputBufferedFile::readb: trying to read beyond the end of the booked memory.\n";
    exit(1);
   }
  };

 private:
  std::ifstream f;
};

class OutputBufferedFile: public BufferedFile
{
 public:
  OutputBufferedFile(std::string fname);
  OutputBufferedFile(std::string fname,unsigned short chunk_size);
  OutputBufferedFile(std::string fname,bool debug);
  OutputBufferedFile(std::string fname,unsigned short chunk_size,bool debug);
  ~OutputBufferedFile();
  // Special function to write strings, which are stored in memory as a first byte, indicating its length (from 0 to 254) followed by the string bytes. No ending \0 is stored.
  void inline writeb(std::string &s)
  {
   size_t sl=s.length();
   if (sl>max_string_len)
   {
    std::cerr << "Low level error from OutputBufferedFile::writeb: trying to write a string with length of " << sl << ", bigger than the maximum allowed, which is " << max_string_len<< std::endl;
    std::cerr << "The offending string is '" << s << "'\n";
    exit(1);
   }
   // If we would be beyond the end of the chunk...
   if (off+sl > memsize)
   {
    // Write the chunk to the disk up to the current position:
    f.write(mem,off);
    // Copy the pending object at the beginning of the chunk.
    // First, a byte with the sring length...
    *mem=(unsigned char)sl;
    // Then, the characters themselves
    memcpy(mem+1,s.c_str(),sl);
    if (deb)
     std::cout << "file " << filename << ": " <<  off << " bytes of string written to reuse memory chunk. Then, " << sl << " bytes stored.\n";
    // and update the offset to the string length, plus 1 to account for the first byte with the length
    off = sl+1;
   }
   else
   {
    // We have still empty room in the memory. Just copy the characters...
    *(mem+off)=(unsigned char)sl;
    memcpy(mem+off+1,s.c_str(),sl);
    // and update the offset
    off += (sl+1);
    if (deb)
     std::cout << "file " << filename << ": " << sl+1 << " bytes of string stored in buffer. Current offset: " << off << std::endl;
   }
  };

  // Function to write any other fundamental data type
  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value, void>::type
  inline writeb(T &v)
  {
   // If we would be beyond the end of the chunk...
   if (off+sizeof(T) > memsize)
   {
    // Write the chunk to the disk up to the current position:
    f.write(mem,off);
    f.flush();
    // Copy the pending object at the beginning of the chunk...
    memcpy(mem,reinterpret_cast<char *>(&v),sizeof(T));
    if (deb)
     std::cout << "file " << filename << ": " <<  off << " bytes written to reuse memory chunk. Then, " << sizeof(T) << " bytes stored.\n";
    // and update the offset by initialization
    off = sizeof(T);
   }
   else
   {
    // We have still empty room in the memory. Just copy the object...
    memcpy(mem+off,reinterpret_cast<char *>(&v),sizeof(T));
    // and update the offset by increment
    off += sizeof(T);
    if (deb)
     std::cout << "file " << filename << ": " << sizeof(T) << " bytes stored in buffer. Current offset: " << off << std::endl;
   }
  };

  // Function to write vectors of any fundamental data type
  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value, void>::type
  inline writeb(std::vector<T> &v)
  {
   size_t sv=v.size();
   // Total number of bytes we have to write: all the vector, plus a size_t to specify its length
   size_t total_size=sv*sizeof(T)+sizeof(size_t);
   // If we would be beyond the end of the chunk...
   if (off+total_size > memsize)
   {
    // Write the chunk to the disk up to the current position:
    f.write(mem,off);
    f.flush();
    // Copy the pending vector at the beginning of the chunk...
    memcpy(mem,reinterpret_cast<char *>(&sv),sizeof(size_t));
    // ... and update the offset by initialization
    off = sizeof(size_t);
    for (size_t i=0; i<sv; i++, off+=sizeof(T))
     memcpy(mem+off,reinterpret_cast<char *>(&(v[i])),sizeof(T));
    // off has already been incremented by the loop update
   }
   else
   {
    // We have still empty room in the memory. Just copy the vector...
    memcpy(mem+off,reinterpret_cast<char *>(&sv),sizeof(size_t));
    // ... and update the offset by increment
    off += sizeof(size_t);
    for (size_t i=0; i<sv; i++, off+=sizeof(T))
     memcpy(mem+off,reinterpret_cast<char *>(&(v[i])),sizeof(T));
    // off has already been incremented by the loop update
    if (deb)
     std::cout << "file " << filename << ": vector of " << total_size << " bytes (" << sv << " elements of " << sizeof(T) << " bytes each) stored in buffer. Current offset: " << off << std::endl;
   }
  };
 private:
  std::ofstream f;
};

#endif
// @endcond
