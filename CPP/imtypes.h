#ifndef _IMTYPES_H
#define _IMTYPES_H

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTIFFImageIO.h>
#include <itkNumericTraits.h>
#include <itkPoint.h>

constexpr unsigned int Planar = 2;

using BinPixelType = unsigned char;
using ImageBinType = itk::Image<BinPixelType, Planar>;
using ImageBinWriterType = itk::ImageFileWriter<ImageBinType>;
using ImageBinReaderType = itk::ImageFileReader<ImageBinType>;

using ColPixelType = itk::RGBPixel<unsigned char>;
using ImageColType = itk::Image<ColPixelType, Planar>;
using ImageColReaderType = itk::ImageFileReader<ImageColType>;
using ImageColWriterType = itk::ImageFileWriter<ImageColType>;

using TIFFIOType = itk::TIFFImageIO;

// These are valid for both types of images, since geometry is the same
using ImageSpType = ImageBinType::SpacingType;
using ImagePointType = ImageBinType::PointType;

using ImageRegionType = itk::ImageRegion<Planar>;
using ImageSizeType = ImageRegionType::SizeType;
using ImageIndexType = ImageRegionType::IndexType;

#endif
