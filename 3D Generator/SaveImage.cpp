/// \file SaveImage.cpp
/// \brief Source code for saving a noise cell as a PNG or TGA file.
///
/// Two versions of the Save2DNoise function, one for Windows and one for Unix.
/// The Windows version uses GDI+ code for saving a PNG image to a file. 
/// This code was mostly scored from the MSDN website, but getting the
/// details right was like pulling teeth. The Unix version saves the noise
/// tile as a TGA image because, frankly, the format is bog-easy and I could
/// write the code in under 30 seconds. No excuses.

// Copyright Ian Parberry, September 2013.
//
// This file is made available under the GNU All-Permissive License.
//
// Copying and distribution of this file, with or without modification,
// are permitted in any medium without royalty provided the copyright
// notice and this notice are preserved.  This file is offered as-is,
// without any warranty.
//
// Created by Ian Parberry, September 2013.
// Last updated May 7, 2014.

#include <stdio.h>

#include "defines.h"

#if defined(_MSC_VER) //Windows Visual Studio only

#include <gdiplus.h>
#include <GdiPlusImageCodec.h>

using namespace Gdiplus;

/// Get CLSID for image file codec.
/// \param format Color format.
/// \param pClsid Pointer to location to save CLSID into.
/// \return -1 for fail, non-negative number otherwise.

int GetEncoderClsid(const WCHAR* format, CLSID* pClsid){
  unsigned int num = 0;  // number of image encoders
  unsigned int size = 0; // size of the image encoder array in bytes

  ImageCodecInfo* pImageCodecInfo = NULL;

  GetImageEncodersSize(&num, &size);
  if(size == 0)return -1;  // Failure

  pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
  if(pImageCodecInfo == NULL)return -1;  // Failure

  GetImageEncoders(num, size, pImageCodecInfo);

  for(unsigned int j=0; j<num; ++j)
    if(wcscmp(pImageCodecInfo[j].MimeType, format) == 0){
      *pClsid = pImageCodecInfo[j].Clsid;
      free(pImageCodecInfo);
      return j;  // Success
    } //if

   free(pImageCodecInfo);
   return -1;  // Failure
} //GetEncoderClsid

#endif

/// Save a cell of 2D noise as a png file (Windows) or TGA file (Unix).
/// \param cell Pointer to array in which to store the noise.
/// \param scale Rescale by this to get in the range -1..1.
/// \param basefilename Root of the file name under which to store the noise.
/// \param n Cell size.

void Save2DNoise(float** cell, const int n, const float scale, const char* basefilename){  
  char filename[MAX_PATH];

#if defined(_MSC_VER) //Windows Visual Studio only

  //add file extension to name
  sprintf(filename, "%s.png", basefilename);

  //convert noise array to png image: preamble
  GdiplusStartupInput gdiplusStartupInput;
  ULONG_PTR gdiplusToken;
  GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

  Bitmap bitmap(n, n); 
  BitmapData bitmapData;

  //pixel data
  bitmap.LockBits(NULL, ImageLockModeWrite, PixelFormat32bppARGB, &bitmapData);
  unsigned int* pixels = (unsigned int*)bitmapData.Scan0;
  const int stride = bitmapData.Stride/4;
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j){
      unsigned int b = (unsigned int)(128.0f*(cell[i][j]*scale + 1.0f));
      pixels[i*stride + j] = 0xFF000000 | b<<16 | b<<8 | b;       
    } //for 
  bitmap.UnlockBits(&bitmapData);

  //convert file name to wide characters
  wchar_t* wcfilename;
  wcfilename = new wchar_t [strlen(filename)+1];
  size_t convertedChars = 0;
  mbstowcs_s(&convertedChars, wcfilename, strlen(filename)+1, filename, _TRUNCATE);

  //save to file
  CLSID pngClsid;
  GetEncoderClsid(L"image/png", &pngClsid);
  bitmap.Save(wcfilename, &pngClsid, NULL);

  delete [] wcfilename;

#else //other OS
  //add file extension to name
  sprintf(filename, "%s.tga", basefilename);
  
  const unsigned char hi = (n & 0xFF00)>>8, lo = n & 0x00FF;

  FILE* fd = fopen(filename, "wb");
  if(fd){
    //TGA header
    fputc(0x00, fd); fputc(0x00, fd); fputc(0x02, fd); 
    for(int i=0; i<9; i++)
   	  fputc(0x00, fd);
   	fputc(lo, fd); fputc(hi, fd);
   	fputc(lo, fd); fputc(hi, fd);
   	fputc(0x20, fd); fputc(0x00, fd);
    
    //pixel data
    for(int i=n-1; i>=0; --i)
      for(int j=0; j<n; ++j){
        unsigned char b = (unsigned char)(128.0f*(cell[i][j]*scale + 1.0f));
      	fputc(b, fd); fputc(b, fd); fputc(b, fd); //blue, green, red
      	fputc(0xFF, fd); //alpha
      } //for

    fclose(fd);
  } //if

#endif //other OS

} //Save2DNoise
