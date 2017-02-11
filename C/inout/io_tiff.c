/*
** Author: Pierre Soille 1988-2006
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef UNIX
#include "bytesex.h" /* <endian.h> only for linux ... */
#endif
#include "tiffio.h"
#include "xtiffio.h" 
#include "geotiff.h" 
#include "imio.h"

/*************************************************************************/
/*                                                                       */
ERROR_TYPE read_image_data(FILE *fp, IMAGE *im, int pc) {
  unsigned long int i, npix, nelem;
  int j, k, y, z, nx, ny, nz, nbpl, nplb, bpp;
  size_t szbuf;
  UCHAR data, temp, *p;
  
  p    = (UCHAR *)GetImPtr(im);
  nx   = GetImNx(im);
  npix = GetImNPix(im);
  
  switch(GetImDataType(im)) {
  case t_TIFFONEBITPERPIXEL: /* assumes 1 x-y plane */
    /* padding for 8 bits as in TIFF files                  */
    /* t_TIFFONEBITPERPIXEL just for I/O                    */
    /* the image is read into a t_UCHAR image               */
    /* note: image_create() takes this feature into account */
    
    SetImDataType(im,t_UCHAR); /* force to t_UCHAR */

    nbpl = nx/BITPERCHAR; /* nbr of bytes per line      */
    nplb = nx%BITPERCHAR; /* nbr of pixels in last byte */
    for (k=0; k<npix; k+=nx) {
      for (i=0; i<nbpl; i++) {
	if ((fread((char *)&data, (int)1, (int)1, fp)) != 1) {
          (void)sprintf(buf,"ERROR in read_image_data(): \
                   unable to read image data block\n"); errputstr(buf);
          return(ERROR);
        }
	temp = data;
	p[(BITPERCHAR*i)+k] = (temp & 0x80) / 0x80;
        for (j=1; j<8; j++){
          temp = data;
          p[(BITPERCHAR*i)+j+k] = ((temp<<j) & 0x80) / 0x80;
	}
      }
      if (nplb){
	if ((fread((char *)&data, (int)1, (int)1, fp)) != 1) {
          (void)sprintf(buf,"ERROR in read_image_data(): \
                   unable to read image data block\n"); errputstr(buf);
          return(ERROR);
	}
        temp = data;
        p[(BITPERCHAR*i)+k] = (temp & 0x80) / 0x80;
        for (j = 1; j < nplb; j++){
          temp = data;
          p[(BITPERCHAR*i)+j+k] = ((temp<<j) & 0x80) / 0x80;
	}
      }
    }
    return NO_ERROR;
  case t_FOURBITPERPIXEL: /* assumes 1 x-y plane */
    /* feature: t_FOURBITPERPIXEL just for I/O              */
    /* the image is read into a t_UCHAR image               */
    /* note: image_create() takes this feature into account */

    SetImDataType(im,t_UCHAR); /* force to t_UCHAR */

    nbpl = nx / 2; /* nbr of bytes per line */
    nplb = nx % 2; /* nbr of pixels in last byte */
    for (k=0; k < npix; k += nx) {
      for (i = 0; i < nbpl; i++) {
	if ((fread((char *)&data, (int)1, (int)1, fp)) != 1) {
	  (void)sprintf(buf,"ERROR in read_image_data(): \
                   unable to read image data block\n"); errputstr(buf);
	  return(ERROR);
	}
	temp = data;
	p[(2*i)+k] = temp>>4;
	p[(2*i)+1+k] = (data <<= 4)>>4;
      }
      if (nplb) {
        if ((fread((char *)&data, (int)1, (int)1, fp)) != 1) {
	  (void)sprintf(buf,"ERROR in read_image_data(): \
                   unable to read image data block\n"); errputstr(buf);
	  return(ERROR);
	}
	p[(2*i)+k] = data>>4;
      }
    }
    return(NO_ERROR);
  case t_UCHAR:
    szbuf=GetImNPix(im)*sizeof(UCHAR);
    break;
  case t_USHORT: /* BUG LITTLE vs BIG ENDIAN not taken into account */
  case t_SHORT:  /* BUG LITTLE vs BIG ENDIAN not taken into account */
    szbuf=GetImNPix(im)*sizeof(USHORT);
    break;
  case t_ULGINT: /* BUG LITTLE vs BIG ENDIAN not taken into account */
  case t_LGINT:  /* BUG LITTLE vs BIG ENDIAN not taken into account */
    szbuf=GetImNPix(im)*sizeof(ULGINT);
    break;
  case t_FLOAT:  /* BUG LITTLE vs BIG ENDIAN not taken into account */
    szbuf=GetImNPix(im)*sizeof(MIAFLOAT);
    break;
  case t_DOUBLE: /* BUG LITTLE vs BIG ENDIAN not taken into account */
    szbuf=GetImNPix(im)*sizeof(DOUBLE);
    break;
  default:
    (void)sprintf(buf,"ERROR in read_image_data(): \
                  invalid data type\n"); errputstr(buf);
    return(ERROR);
  }
  if (pc==PLANARCONFIG_CONTIG && GetImNz(im)>1){ /* interleaved x-y plane strorage (BIP) */
    nelem = GetImNPixPerPlane(im);
    p = (UCHAR *)GetImPtr(im);
    nz=GetImNz(im);
    bpp=GetImBitPerPixel(im)/8;
    for (i=0; i<nelem; i++) { /* read each successive channel values */
      for (z=0;z<nz;z++){
	fread((void *)p+(i+z*nelem) * bpp, bpp, 1, fp);
      }
    }
    /* SetImDataType(im,t_UCHAR);
    SetImNz(im,1);
    SetImNy(im,3*GetImNy(im)); */
    return(NO_ERROR);
  }
  else if (pc==3 && GetImNz(im)>1) { /* interleaved by line (BIL) */
    nelem = GetImNPixPerPlane(im);
    p = (UCHAR *)GetImPtr(im);
    nx=GetImNx(im);
    ny=GetImNy(im);
    nz=GetImNz(im);
    bpp=GetImBitPerPixel(im)/8;
    for(y=0;y<ny;y++){
      for(z=0;z<nz;z++){
	fread((void *)p+(y*nx+z*nelem) * bpp, bpp, nx, fp);
      }
    }
  } 
  else { /* assume sequential x-y plane storage (BSQ) */p = (UCHAR *)GetImPtr(im);
    p = (UCHAR *)GetImPtr(im);
    nx=GetImNx(im);
    ny=GetImNy(im);
    nz=GetImNz(im);
    bpp=GetImBitPerPixel(im)/8;
    for(z=0;z<nz;z++){
      for(y=0;y<ny;y++){
	fread((void *)p, (size_t)bpp, (size_t)nx, fp);
	p+=nx*bpp;
      }
    }

/* Before Feb 28 2006: */	
/*     if ((fread((void *)GetImPtr(im), (size_t)1, szbuf, fp)) !=szbuf) { */
/*       (void)sprintf(buf,"ERROR in read_image_data(): \ */
/*                    unable to read image data block\n"); errputstr(buf); */
/*       return(ERROR); */
/*     } */
  }
  return(NO_ERROR);
}

/*
**  Function to read TIFF, VISILOG, or KIFF image files.
*/
IMAGE *read_image(char *fn) {
  FILE *fp;
  IMAGE *im;
  TIFF *tiffp;
  unsigned long int nbyte=0, bread=0;
  uint32 *s_o_t, rps;
  uint16 spp=1, pc=1;
  uint32 fsot; /* offset to 1st pixel */
  unsigned short int *tred, *tgreen, *tblue;
  unsigned short int *red, *green, *blue;
  tstrip_t nstrip, strip;
  unsigned short int pmi;
  int cm=1; /* 1 if Color Map, 0 otherwise */
  unsigned short int *plut;

  unsigned char visilog = FALSE;
  LGINT visi[19];
  LGINT *pl1;
  long int i;
  int nx, ny, nz = 1, data_type=t_UNSUPPORTED;
  unsigned short int ctype; /* compression type */
  short int bitpp, sf;
  unsigned short kiff_head[6];
  IMAGE *imtmp=NULL;

  /* try first a tiff file using libtiff library */
  if ((tiffp = TIFFOpen(fn, "rc")) != NULL) {
    TIFFGetField(tiffp, TIFFTAG_IMAGEWIDTH, &nx);
    TIFFGetField(tiffp, TIFFTAG_IMAGELENGTH, &ny);
    TIFFGetField(tiffp, TIFFTAG_BITSPERSAMPLE, &bitpp);
    TIFFGetField(tiffp, TIFFTAG_SAMPLESPERPIXEL, &spp);
    if ( TIFFGetField(tiffp, TIFFTAG_SAMPLEFORMAT, &sf) != 1) 
      sf=1; /* 1 unsigned (default value in TIFF 6.0 specifications)
	       2 signed, 3 float, 4 undefined */
    if (spp>1)
      TIFFGetField(tiffp, TIFFTAG_PLANARCONFIG, &pc);

    if ((bitpp == BITPERCHAR) && (sf==1) )
      data_type = t_UCHAR;
    else if (bitpp == BITPERSHORT){
      if (sf==1)
	data_type = t_USHORT;
      else if (sf==2)
	data_type = t_SHORT;
    }
    else if ( (bitpp == BITPERLGINT) && (sf==3) )
      data_type = t_FLOAT;
    else if ( (bitpp == BITPERLGINT) && (sf==1) )
      data_type = t_ULGINT;
    else if ( (bitpp == BITPERLGINT) && (sf==2) )
      data_type = t_LGINT;
    else if (bitpp == BITPERDOUBLE)
      data_type = t_DOUBLE;
    if (data_type!=t_UNSUPPORTED) {
      TIFFGetField(tiffp, TIFFTAG_COMPRESSION, &ctype);
      nstrip = TIFFNumberOfStrips(tiffp);
      TIFFGetField(tiffp, TIFFTAG_ROWSPERSTRIP, &rps);	
    
      /*     if ((nstrip != 1 && ctype==1) || ctype==5) {*/
      if (nstrip != 1 || ctype != 1) {
#ifdef DEBUG
	(void)sprintf(buf,"TIFFTAG_ROWSPERSTRIP=%d\n", (int)rps); errputstr(buf);
	(void)sprintf(buf,"nstrip=%d\n", (int)nstrip); errputstr(buf);
	(void)sprintf(buf,"read_image(): number of strips exceeds 1 or compression, \
                           I shall try to read it using TIFFReadEncodedStrip\n");
	errputstr(buf);
#endif
	im = create_image(data_type, nx, ny, spp);
	if (im == NULL) {
	  (void)sprintf(buf,"read_image(): not enough memory!\n"); errputstr(buf);
	  TIFFClose(tiffp);
	  return(NULL);
	}
	for (strip = 0; strip < nstrip-1; strip++){
	  nbyte=TIFFReadEncodedStrip(tiffp, strip, GetImPtr(im)+bread, (tsize_t) -1);
	  bread+=nbyte;
	  if (nbyte==-1){
	    (void)sprintf(buf,"read_image(): problem when reading strip number %d\n", (int)strip); errputstr(buf);
	    TIFFClose(tiffp);
	    return(NULL);
	  }
	}
	/* read last strip (possibly not containing a full set of rows!)*/
	if (bread<nx*ny*spp*(int)(bitpp/8)) {
	  nbyte=TIFFReadEncodedStrip(tiffp, strip, GetImPtr(im)+bread,  (tsize_t) nx*ny*spp*(int)(bitpp/8)-bread);
	  bread+=nbyte;
	  if (nbyte==-1){
	    (void)sprintf(buf,"read_image(): problem when reading last strip\n"); errputstr(buf);
	    TIFFClose(tiffp);
	    return(NULL);
	  }
	}
#ifdef DEBUG
	(void)sprintf(buf,"strip= %d\t nbyte read in last strip=%d\n", (int) strip, (int)nbyte); errputstr(buf);
#endif
	TIFFClose(tiffp);
	if ((spp>1) && (pc==1)) { /* deinterleave */
	  imtmp= deinterleave(im);
	  if (imtmp == NULL) {
	    (void)sprintf(buf,"warning: read_image(): not enough memory for deinterleaving!\n"); errputstr(buf);
	    return im;
	  }
	  free_image(im);
	  return imtmp;
	}
	return im;
      }
    }
    TIFFGetField(tiffp, TIFFTAG_STRIPOFFSETS, &s_o_t);
    fsot = *s_o_t;
    TIFFGetField(tiffp, TIFFTAG_BITSPERSAMPLE, &bitpp);
    red = (uint16 *)malloc(1<<bitpp*sizeof(uint16));
    green = (uint16 *)malloc(1<<bitpp*sizeof(uint16));
    blue = (uint16 *)malloc(1<<bitpp*sizeof(uint16));
    TIFFGetField(tiffp, TIFFTAG_PHOTOMETRIC, &pmi);
    if (TIFFGetField(tiffp, TIFFTAG_COLORMAP, &tred, &tgreen, &tblue) != 1)
      cm=0;
    else {
      for (i=0; i<1<<bitpp; i++) {
	red[i]=tred[i];
	green[i]=tgreen[i];
	blue[i]=tblue[i];
      }
    }
          
    TIFFClose(tiffp);

    if (bitpp == 1)
      data_type = t_TIFFONEBITPERPIXEL;
    else if (bitpp == 4)
      data_type = t_FOURBITPERPIXEL;

    /* Allow for color images  */
    if (pmi == 2){
      /* data_type = t_RGB; t_RGB not used anymore for input */
      nz = spp;
    }
  
    if ((fp = fopen(fn, "rb")) == NULL) {
      free(red); free(green); free(blue);
      return(NULL);
    }
   
    (void)fseek(fp, fsot, 0);  /*  position file pointer to pixel map  */

    /* create output image */
    im = create_image(data_type, nx, ny, nz);
    if (im == NULL) {
      (void)sprintf(buf,"read_image(): not enough memory!\n"); errputstr(buf);
      free(red); free(green); free(blue);
      fclose(fp);
      return(NULL);
    }
    /*  Read image data  */
    if (read_image_data(fp, im, pc) != NO_ERROR){
      free_image(im);
      free(red); free(green); free(blue);
      fclose(fp);
      return(NULL);
    }
    fclose(fp);
    /* set Color Map if any */
    if (cm==1) {
      plut = (unsigned short int *)malloc(3*256*sizeof(short)); /* should be (1<<bitpp) instead of 256 */
      SetImLut(im,plut);
      if (plut!=NULL) {
	for (i=0; i<(1<<bitpp); i++) {
	  plut[i]=red[i];
	  plut[i+256]=green[i];
	  plut[i+2*256]=blue[i];
	}
      }
    }
    else
      SetImLut(im,NULL);
    free(red); free(green); free(blue);
    return(im);
  }
  
  /*  alternatively, open first file  */
  if ((fp = fopen(fn, "rb")) == NULL) {
    (void)sprintf(buf,"ERROR in read_image(\"%s\"): \
                   unable to read open file\n", fn); errputstr(buf);
    return(NULL);
  }
  
  /*  Check whether it is a VISILOG file  */
  if (fread((char *)visi, 76, 1, fp) != 1) {
    (void) fclose(fp);
    return(NULL);
  }
  if (visi[0] == 0x6931)  
    visilog = TRUE;
  else if (visi[0] == 0x31690000) {
    /*  Swap header  */
    pl1  = visi;
    for (i = 0; i < 19; ++i)
      swap_long(pl1++);
    visilog = TRUE;
  }

  /*  Else, check whether it is a TIFF or KIFF file  */  
  if (!visilog) {  /* try a KIFF file */
    (void)fseek(fp, 0L, 0);
    if (fread((char *)&kiff_head[0], sizeof(kiff_head), 1, fp) != 1) {
      (void) fclose(fp);
      return(NULL);
    }
    if (kiff_head[1] == 0x4712 && kiff_head[2] == 0x6DB0) {
      nx = ((kiff_head[3]&0xFF)<<8)|(kiff_head[3]>>8);
      ny = ((kiff_head[4]&0xFF)<<8)|(kiff_head[4]>>8);
      (void) fclose(fp);
      return(read_all(fn, nx, ny, 1, t_UCHAR, 128, 1));
    }
    else {
      (void)sprintf(buf,"read_file(): unable to read %s on disk\n", fn); errputstr(buf);
      im = NULL;
    }
  }    
  else { /* it is a VISILOG file */
    nx    = visi[1];
    ny    = visi[2];
    nz    = visi[3];
    bitpp = visi[9];
    if (bitpp == 8)
      data_type = t_UCHAR;
    else if (bitpp == 16)
      data_type = t_USHORT;
    else if (bitpp == 32)
      data_type = t_LGINT;
    (void) fclose(fp);
    return(read_all(fn, nx, ny, nz, data_type, 76, 1));
  }

  (void) fclose(fp);
  return(im);
}


/*****************************************************************************/

ERROR_TYPE write_image_data(FILE *fp, IMAGE *im, int pc)
{
  unsigned long int nelem, i, y, z;
  UCHAR *p, *plast;
  long int bpp, nx, ny, nz = GetImNz(im);
 
  switch (GetImDataType(im)){
  case t_TIFFONEBITPERPIXEL:
  case t_ONEBITPERPIXEL:
  case t_FOURBITPERPIXEL:
  case t_UCHAR:
  case t_USHORT:
  case t_SHORT:
  case t_ULGINT:
  case t_LGINT:
  case t_FLOAT:
  case t_DOUBLE:
    if (GetImNz(im)>1) {
      if (pc==PLANARCONFIG_CONTIG){ /* bip model (band interleaved by pixel) TIFF default */
	/* interleave x-y planes by pixel (bip) */  
	bpp=GetImBitPerPixel(im)/8;
	nelem=GetImNPixPerPlane(im);
	p=(UCHAR *)GetImPtr(im);
	for (i=0; i<nelem; i++) { /* write each successive channel values */
	  for (z=0;z<nz;z++){
	    fwrite(p+(i+z*nelem) * bpp, bpp, 1, fp);
	  }
	}
      }
/*	plast=p+nelem*bpp; */
/* 	for (; p<plast; p++){ */
/* 	  for (i=0; i<nz; i++){ */
/* 	    printf("val=%f\n", (float *)p); */
/* 	    fwrite(p + (i*nelem) * bpp, bpp, 1, fp); */
/* 	  }   */
/* 	} */
/*       } */
      else if (pc==PLANARCONFIG_SEPARATE){ /* bsq model (band sequential) */
	if (fwrite((char*)GetImPtr(im), 1, (int)GetImNByte(im), fp) != GetImNByte(im)) {
	  (void)sprintf(buf, "write_image_data(): cannot write data on disk\n"); errputstr(buf);
	  return(ERROR);
	}
      }
      else if (pc==PLANARCONFIG_BIL){ /* bil model (band interleaved by line) not supported by TIFF!!! */
	nelem = GetImNPixPerPlane(im);
	p = (UCHAR *)GetImPtr(im);
	nx=GetImNx(im);
	ny=GetImNy(im);
	nz=GetImNz(im);
	bpp=GetImBitPerPixel(im)/8;
	for(y=0;y<ny;y++){
	  for(z=0;z<nz;z++){
	    fwrite((char *)p+(y*nx+z*nelem) * bpp, bpp, nx, fp);
	  }
	}
      }
      else{
	(void)sprintf(buf, "write_image_data(): invalid planar configuration\n"); errputstr(buf);
	return(ERROR);
      }
    }
    else { /* bsq, bip, and bil models are identical */
      if (fwrite((char*)GetImPtr(im), 1, (int)GetImNByte(im), fp) != GetImNByte(im)) {
	(void)sprintf(buf, "write_image_data(): cannot write data on disk\n"); errputstr(buf);
	return(ERROR);
      }
    }
    break;
  case t_RGB:
    nelem=GetImNPixPerPlane(im);
    p=(RGB_TYPE *)GetImPtr(im);
    plast=p+nelem;
    for (; p<plast; p++){ /* bil interleaved (PlanarConfiguration=1) is default and recommended */
      for (i=0; i<GetImNz(im); i++)
	fwrite((char *)(p+i*nelem), sizeof(RGB_TYPE), 1, fp);
    }
    break;
  default:
    (void)sprintf(buf,"write_image_data(): invalid pixel type\n"); errputstr(buf);
    return(ERROR);
  }
  return(NO_ERROR);
}


/*************************************************************************/
/*                                                                       */

ERROR_TYPE write_tiff(IMAGE *im, char *fn)
{
  size_t nbyte=0;
  FILE *fp;
  IMAGE *imtmp;
  UCHAR uc_maxi;
  char doc_name[32];
  char doc_info[256] = "File created using Pierre Soille's library of independent image analysis routines";
  short int nbr_tags = NBR_TAGS;
  struct mytiff tiff_head;
  struct tag *tag_info;
  int pmi=1, spp=1, maxi=255;
  int i, sf;
  int xr_n = 1, xr_d = 1, yr_n = 1, yr_d = 1, swapflag=0;
  short int eight=8;
  G_TYPE *pg=NULL;

  if (GetImNz(im)==1 && GetImDataType(im)==t_UCHAR && GetImLut(im) != NULL) { /* palette COLOUR TIFF */
      return(write_ColorMap_tiff(im, fn));
  }

  /*   if (GetImNz(im)==1 || GetImDataType(im)==t_RGB) */
    nbr_tags -=1;
  
  /*  Allocate memory for header  */  
  if ((tag_info = (struct tag *)calloc((unsigned)nbr_tags, (unsigned)sizeof(struct tag))) == NULL) {
    (void)sprintf(buf,"write_tiff_file(): not enough memory\n"); errputstr(buf);
    return(ERROR);
  }
    
  /*  Header initialisation  */
#if BYTE_ORDER==BIG_ENDIAN
  tiff_head.byte_order  = 0x4d4d; /* big-endian byte ordering */
#else /* must be little endian */
  tiff_head.byte_order  = 0x4949; /* little-endian byte ordering */
#endif
  tiff_head.version     = 42;     /*  TIFF version  0x2a            */
  tiff_head.ptr_dir1    = 8;      /*  Pointer to first directory    */
  
  /*  Tag type initialisation */
  tag_info[WIDTH].type  = 256;  /*  ImageWidth tag or number of columns  */
  tag_info[LENGTH].type = 257;  /*  ImageLength tag or number of lines */
  tag_info[BPS].type    = 258;  /*  BitPerSamples tag          */
  tag_info[PMI].type    = 262;  /*  PhotometricInterpretation  tag */
  tag_info[NAME].type   = 269;  /*  DocumentNameTag           */
  tag_info[DESC].type   = 270;  /*  ImageDescription Tag      */
  tag_info[SOT].type    = 273;  /*  StripOffsets  tag         */
  tag_info[SPP].type    = 277;  /*  SamplesPerPixel  tag      */
  tag_info[SBC].type    = 279;  /*  StripByteCounts */
  tag_info[XR].type     = 282;  /*  XResolution  tag          */
  tag_info[YR].type     = 283;  /*  YResolution  tag          */
  tag_info[SF].type     = 339;  /*  SampleFormat tag          */

/*   if (GetImNz(im)!=1 && GetImDataType(im)!=t_RGB) */
/*     tag_info[N_PLANE].type = 32768;  /\*  ZResolution  tag  (for 3D images) *\/ */
/*                                      /\*  NON-STANDART  TAG  !!!!           *\/ */
  
  /*  Data Type initialisation */
  tag_info[WIDTH].data_type  = 4;  /* 4  ==>  LGINT            */
  tag_info[LENGTH].data_type = 4;  
  tag_info[BPS].data_type    = 3;  /* 3  ==>  SHORT           */    
  tag_info[PMI].data_type    = 3;  
  tag_info[NAME].data_type   = 2;  /* 2  ==>  8-bits ASCII codes */  
  tag_info[DESC].data_type   = 2; 
  tag_info[SOT].data_type    = 4;
  tag_info[SPP].data_type    = 3;  
  tag_info[SBC].data_type    = 4;
  tag_info[XR].data_type     = 5;  /* 5  ==>  RATIONAL           */
  tag_info[YR].data_type     = 5;
  tag_info[SF].data_type     = 3;
/*   if (GetImNz(im)!=1 && GetImDataType(im)!=t_RGB) */
/*     tag_info[N_PLANE].data_type  = 3;   */
  
  /*  Length initialisation  */
  tag_info[WIDTH].length  = 1;    
  tag_info[LENGTH].length = 1;    
  tag_info[BPS].length    = 1;    
  tag_info[PMI].length    = 1;  
  tag_info[NAME].length   = 32;    /*  Document name      */  
  tag_info[DESC].length   = 256;   /*  Document information  */
  tag_info[SOT].length    = 1;
  tag_info[SPP].length    = 1;    
  tag_info[SBC].length    = 1;
  tag_info[XR].length     = 1;    
  tag_info[YR].length     = 1;   
/*   if (GetImNz(im)!=1 && GetImDataType(im)!=t_RGB) */
/*     tag_info[N_PLANE].length  = 1;     */

  /* take minimum and maximum image values into account */
  nbyte=(GetImBitPerPixel(im)/8) * GetImNPix(im);  //GetImNByte(im);
  if (GetImDataType(im)==t_UCHAR){
    pg = min_max(im);
    if (pg != NULL){
      uc_maxi = pg[1].uc_val;
      free(pg);
      if (uc_maxi < 2) {
	printf("converting to 1 bit per pixel\n");
	imtmp = to_tiff1bitpp(im);
	if (imtmp == NULL) {
	  free((char *)tag_info);
	  return(ERROR);
	}
	im = imtmp;
	swapflag = 1;
	maxi = 1;
	nbyte=GetImNy(im)*(GetImNx(im)/8+(GetImNx(im)%8 ? 1: 0))*sizeof(UCHAR);
	printf("nbyte=%d\n", nbyte);
      }
      else if (uc_maxi < 16) {
	printf("converting to 4 bits per pixel\n");
	imtmp = to_tiff4bitpp(im);
	if (imtmp == NULL) {
	  free((char *)tag_info);
	  return(ERROR);
	}
	im = imtmp;
	swapflag = 1;
	maxi = 15;
	nbyte=GetImNy(im)*(GetImNx(im)/2+(GetImNx(im)%2))*sizeof(UCHAR);
	printf("nbyte=%d\n", nbyte);
      }
    }
  }
 
/*   if (GetImDataType(im)==t_RGB){ /\* full color image *\/    */
/*     tag_info[BPS].length  = GetImNz(im); */
/*     pmi=2; */
/*     spp=GetImNz(im); /\* may be more for multispectral images *\/ */
/*   } */
    
  if (GetImNz(im)>1){ /* multichannel or 3-D image */   
    tag_info[BPS].length  = GetImNz(im);
    pmi=2;
    spp=GetImNz(im);
  }
     
  tag_info[SF].length     = spp;

  switch (GetImDataType(im)) { /* set SampleFormat field value */
      case t_UCHAR:
      case t_USHORT:
      case t_ULGINT:
      case t_RGB:
	  sf=1; /* unsigned integer data */
	  break;
      case t_SHORT:
      case t_LGINT:
	  sf=2; /* 2's complement signed integer data */
	  break;
      case t_FLOAT:
      case t_DOUBLE:
	  sf=3; /* IEEE floating point data */
	  break;
      default:
	  sf=4; /* undefined data format */
  }

  /*  Values fitting in tag initialisation  */
#if BYTE_ORDER==BIG_ENDIAN
  tag_info[WIDTH].ValOrPoint  = GetImNx(im);    
  tag_info[LENGTH].ValOrPoint = GetImNy(im);
  tag_info[BPS].ValOrPoint    = GetImBitPerPixel(im)<<16;
  tag_info[PMI].ValOrPoint    = pmi<<16;
  tag_info[SOT].ValOrPoint    = 10+nbr_tags*12+4+32+256+16;
  tag_info[SPP].ValOrPoint    = spp<<16;
  tag_info[SF].ValOrPoint     = sf<<16;
  tag_info[SBC].ValOrPoint    = nbyte; //GetImNByte(im); //GetImBitPerPixel(im)/8 * GetImNPix(im);
/*   if (GetImNz(im)!=1 && GetImDataType(im)!=t_RGB) */
/*     tag_info[N_PLANE].ValOrPoint = GetImNz(im)<<16; */

  /*  Offset of values not fitting in tag initialisation  */  
  tag_info[NAME].ValOrPoint = 10+nbr_tags*12+4;
  tag_info[DESC].ValOrPoint = 10+nbr_tags*12+4+32;
  tag_info[XR].ValOrPoint   = 10+nbr_tags*12+4+32+256;
  tag_info[YR].ValOrPoint   = 10+nbr_tags*12+4+32+256+8;
  if (pmi==2){ /* reset tag 2 and 6 */
   tag_info[BPS].ValOrPoint = 10+nbr_tags*12+4+32+256+8+8;
   tag_info[SOT].ValOrPoint = 10+nbr_tags*12+4+32+256+8+8+GetImNz(im)*2;
  }
  if (spp>1){ 
   tag_info[SF].ValOrPoint  = 10+nbr_tags*12+4+32+256+8+8+GetImNz(im)*2;
   tag_info[SOT].ValOrPoint += GetImNz(im)*2;
  }
#elif BYTE_ORDER==LITTLE_ENDIAN
  tag_info[WIDTH].ValOrPoint  = GetImNx(im);    
  tag_info[LENGTH].ValOrPoint = GetImNy(im);
  tag_info[BPS].ValOrPoint    = (short)GetImBitPerPixel(im);
  tag_info[PMI].ValOrPoint    = (short)pmi;
  tag_info[SOT].ValOrPoint    = 10+nbr_tags*12+4+32+256+16;
  tag_info[SPP].ValOrPoint    = (short)spp;
  tag_info[SF].ValOrPoint     = (short)sf;
  tag_info[SBC].ValOrPoint    = nbyte; //GetImNByte(im); //GetImBitPerPixel(im)/8 * GetImNPix(im);
/*   if (GetImNz(im)!=1 && GetImDataType(im)!=t_RGB) */
/*     tag_info[N_PLANE].ValOrPoint = (short)GetImNz(im); */
  
  /*  Offset of values not fitting in tag initialisation  */  
  tag_info[NAME].ValOrPoint = 10+nbr_tags*12+4;
  tag_info[DESC].ValOrPoint = 10+nbr_tags*12+4+32;
  tag_info[XR].ValOrPoint   = 10+nbr_tags*12+4+32+256;
  tag_info[YR].ValOrPoint   = 10+nbr_tags*12+4+32+256+8;
  if (pmi==2){ /* reset tag 2 and 6 */
   tag_info[BPS].ValOrPoint = 10+nbr_tags*12+4+32+256+8+8; /* reset */
   tag_info[SOT].ValOrPoint = 10+nbr_tags*12+4+32+256+8+8+GetImNz(im)*2;
  }
  if (spp>1){ 
   tag_info[SF].ValOrPoint  = 10+nbr_tags*12+4+32+256+8+8+GetImNz(im)*2;
   tag_info[SOT].ValOrPoint += GetImNz(im)*2;
  }

#else
 #error BYTE_ORDER must be either BIG_ENDIAN or LITTLE_ENDIAN
#endif
  /*  Offset to next directory (if any) initialisation  */
  ptr_dir2 = 0;  

  /*  Open output file  */
  if ((fp = fopen(fn, "wb")) == NULL){
    (void)sprintf(buf, "write_tiff(): unable to open output file\n"); errputstr(buf);
    if (swapflag)
      free_image(im);
    return(ERROR);
  }
  
  /*  Write header  */
  (void) fwrite((char *)&tiff_head, sizeof(tiff_head), 1, fp);
  (void) fwrite((char *)&nbr_tags, 2, 1, fp);
  (void) fwrite((char *)tag_info, 12*nbr_tags, 1, fp);
  free((char *)tag_info);
  (void) fwrite((char *)&ptr_dir2, 4, 1, fp);
  (void) fwrite((char *)doc_name, sizeof(doc_name), 1, fp);
  (void) fwrite((char *)doc_info, sizeof(doc_info), 1, fp);
  (void) fwrite((char *)&xr_n, 4, 1, fp);
  (void) fwrite((char *)&xr_d, 4, 1, fp);
  (void) fwrite((char *)&yr_n, 4, 1, fp);
  (void) fwrite((char *)&yr_d, 4, 1, fp);
  if (pmi==2) {
    for (i=0;i<GetImNz(im); i++)
      (void) fwrite((char *)&eight, 2, 1, fp);
  }
  if (spp>1) {
    for (i=0;i<GetImNz(im); i++)
      (void) fwrite((char *)&sf, 2, 1, fp);
  }

  /*  Write image data  */
  if (write_image_data(fp, im, PLANARCONFIG_CONTIG) != NO_ERROR){
    (void)sprintf(buf, "write_tiff(): unable to write \"%s\" on disk\n", fn); errputstr(buf);
    (void)fclose(fp);
    if (swapflag)
      free_image(im);
    return(ERROR);
  }

  (void) fclose(fp);
  if (swapflag)
    free_image(im);

  return(NO_ERROR);
}

