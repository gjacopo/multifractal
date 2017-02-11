

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <utils.h>

#include <utl_alloc.h>	
#include <utl_char.h>	
#include <utl_stats.h>
#include <utl_operator.h>	
	
#include <flt_stats1d.h>
#include <flt_stats2d.h>

#include <inout.h>  	
#include <io_raw.h>  	


/***************************************************************************/
int read_dimensiones_raw( char *name,
			  int *bd, int *dimx, int *dimy, int *dimv, int *dimz ) {
  /***************************************************************************/
  FILE *canal;
  char *strhdr;
  int il;
  
  for( il=strlen(name)-1; (il>=0)&&(name[il]!='.'); il-- );
  strhdr[il]='\0';
  sprintf(strhdr,"%s.hdr",strncpy( strhdr, name, il));
  
  TrackNullAlloc( canal=fopen(strhdr,"rt") );
  fscanf(canal,"%d",bd);
  fscanf(canal,"%d",dimx);
  fscanf(canal,"%d",dimy);
  fscanf(canal,"%d",dimv);
  fscanf(canal,"%d",dimz);
  fclose(canal);
  
  return OK;
} // end of read_dimensiones_raw


/***************************************************************************/
int read_raw( char *name, 
	      int bd, int dimx, int dimy, int dimv, int dimz, int iz, 
	      double ***data) {
  /***************************************************************************/
  
  int levels;
  
  TrackError( (levels=read_winraw( name, bd, dimx, dimy, dimv, dimz, iz,
				   0, 0, dimx, dimy, data )),
	      "Error reading image window" );
  
  return levels;
} // end of read_raw


/***************************************************************************/
int read_winraw( char *name, 
		 int bd, int dimx, int dimy, int dimv, int dimz, int iz, 
		 int ix0, int iy0, int ndimx, int ndimy, 
		 double ***data) {
  /***************************************************************************/
  
  FILE *canal;
  float value;
  int dat,levels;
  int ix,iy,iv,ib;
  
  TrackNull( (canal=fopen(name,"rb")), "Error opening file in 'r' mode" );
  
  TrackError( fseek(canal,(iz-dimz)*dimx*dimy*dimv*bd,SEEK_END), 
	      "Error seeking file" );
  
  for( iy=0; iy<dimy; iy++ )
    for( ix=0; ix<dimx; ix++ )     
      for( iv=0; iv<dimv; iv++ )	{
	
	if(bd<4)	  {
	  value=0.;
	  for(ib=0;ib<bd;ib++)	    {
	    TrackEOF( dat=(int)getc(canal), "Error reading file input" );
	    if(dat<0) dat+=256;
	    if(LITTLE_INDIAN) value+=pow(256.,(float)ib)*dat;
	    else value=256*value+dat;
	  }
	  
	} else fread(&value,sizeof(float),1,canal);
	
	if((ix>=ix0)&&(iy>=iy0)&&(ix<ix0+ndimx)&&(iy<iy0+ndimy))
	  data[iv][iy-iy0][ix-ix0]=value;
      }
  
  fclose(canal);
  if(bd<4) levels=(int) (pow(256.,(double)bd));
  else levels=1;
  
  return levels;
} // end of read_winraw


/***************************************************************************/
int read_en_double_int( int xmax, int ymax, int ix0, int iy0, 
			int dimx, int dimy,
			int littleindian, char* nombre, 
			double **data )/*lee_double_int */ {
/***************************************************************************/

  FILE* canal;
  double value;
  int dat1,dat2,ix,iy;
  
  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );

  for( iy=0; iy<ymax; iy++ )    
    for( ix=0; ix<xmax; ix++ )	{
      dat1=(int) getc(canal);
      if(dat1<0) dat1=dat1+256;
      dat2=(int) getc(canal);
      if(dat2<0) dat2=dat2+256;
      
      if(littleindian) value=(double)(dat1+256*dat2);
      else value=(double)(dat2+256*dat1);
      
      if((ix>=ix0)&&(ix<ix0+dimx)&&(iy>=iy0)&&(iy<iy0+dimy))
	data[dimy-1-(iy-iy0)][ix-ix0]=value;
    }
  
  fclose(canal);
  
  return 32768 ;
} // end of read_en_double_int


/***************************************************************************/
int check_unformatted( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		       char *nombre, Read2D *p_lee ) {
/***************************************************************************/

  FILE *canal;
  long int length;
  int error=OK;

  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );

  fseek(canal,0,SEEK_SET);
  length=-ftell(canal);
  fseek(canal,0,SEEK_END);
  length+=ftell(canal);
  fclose(canal);
  if(length==3145728)     {
    *dimx=XMAXVH;
    *dimy=YMAXVH;
    *dimv=1;
    *dimz=1;
    *bd=2;
  }  else if(length==1242000)    {
    *dimx=XMAXINRIA;
    *dimy=YMAXINRIA;
    *dimv=1;
    *dimz=1;
    *bd=1;
  }
  else error=ERROR;
  if(!error) *p_lee=&read_unformatted;
  
  return(error);
} // end of check_unformatted


/***************************************************************************/
int read_unformatted( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
		      char *nombre, double ***data ) {
  /***************************************************************************/
  
  int levels;
  
  if((dimx==XMAXVH)&&(dimy==YMAXVH))
    levels = read_en_double_int(XMAXVH,YMAXVH,0,0,dimx,dimy,0,nombre,data[0]);
  else 
    levels = read_en_double_int(XMAXINRIA,YMAXINRIA,0,0,dimx,dimy,1,
			     nombre,data[0]);
  
  return levels;
} // end of read_unformatted


/***************************************************************************/
int ext_unformatted( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		     char *nombre, Read2D *p_lee ) {
  /***************************************************************************/
  char ext[MAXNAMELENGTH];
  int error=OK;
  int lext=extract_extension(nombre,'.',ext);
  
  if(lext!=3) error=ERROR;
  else if(
	  ((ext[0]=='i')||(ext[0]=='I'))&&
	  ((ext[1]=='m')||(ext[1]=='M'))&&
	  ((ext[2]=='c')||(ext[2]=='C'))
	  ) error=OK;
  else if(
	  ((ext[0]!='r')&&(ext[0]!='R'))||
	  ((ext[1]!='a')&&(ext[1]!='A'))||
	  ((ext[2]!='w')&&(ext[2]!='W'))
	  ) error=ERROR;
  
  if(!error) error=check_unformatted(dimx,dimy,dimv,dimz,bd,nombre,p_lee);

  return error;
} // end of ext_unformatted


/***************************************************************************/
int read_flow( int dimx, int dimy, char *nombre_out, double **data )/*lee*/ {
  /***************************************************************************/
  FILE* canal;
  int ix,iy;
  
  TrackNull( canal=fopen(nombre_out,"rb"), "Error opening file in 'r' mode" );

  for(iy=0;iy<dimy;iy++) fread(data[iy],sizeof(double),dimx,canal);
  fclose(canal);

  return OK;
} // end of read_flow


/***************************************************************************/
int write_flow( int dimx, int dimy, char* nombre_out,double **data) {
/***************************************************************************/
  FILE* canal;
  int ix,iy;

  TrackNull( canal=fopen(nombre_out,"wb"), "Error opening file in 'w' mode" );

  for(iy=0;iy<dimy;iy++) fwrite(data[iy],sizeof(double),dimx,canal);
  fclose(canal);

  return OK;
} // end of write_flow


/***************************************************************************/
int read_file( int leff, int dimy, FILE *canal, double **signal) {
  /***************************************************************************/
  int iy;
  for(iy=0;iy<dimy;iy++) fread(signal[iy],sizeof(double),leff,canal);
  return OK;
} // end of read_file


/***************************************************************************/
int read_data( int leff, int dimy, char *nombre_in, double **signal)/*lee_datos*/ {
  /***************************************************************************/
  FILE *canal;
  char nombre[MAXNAMELENGTH];
  int iy;
  
  sprintf(nombre,"%s.dat",nombre_in);
  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );
  
  for(iy=0;iy<dimy;iy++) fread(signal[iy],sizeof(double),leff,canal);
  fclose(canal);

  return OK;
} // end of read_data


/***************************************************************************/
int write_data( int leff, int dimy, char *nombre_in, double **signal)
  /*graba_datos*/ {
  /***************************************************************************/
  FILE *canal;
  char nombre[MAXNAMELENGTH];
  int iy;
  
  sprintf(nombre,"%s.dat",nombre_in);
  TrackNull( canal=fopen(nombre,"wb"), "Error opening file in 'w' mode" );
  
  for(iy=0;iy<dimy;iy++) fwrite(signal[iy],sizeof(double),leff,canal);
  fclose(canal);

  sprintf(nombre,"%s.hdr",nombre_in);
  TrackNull( canal=fopen(nombre,"wt"), "Error opening file in 'w' mode" );
  
  fprintf(canal,"%s.dat\n%d",nombre_in,leff);
  fclose(canal);

  return OK;
} // end of write_data


/***************************************************************************/
int write_file( int leff, int dimy, FILE *canal, double **signal) {
  /***************************************************************************/
  int iy;
  for(iy=0;iy<dimy;iy++) fwrite(signal[iy],sizeof(double),leff,canal);
  return OK;
} // end of write_file


/***************************************************************************/
int read_data_infloat( int leff, int dimy, char *nombre_in, double **signal)
  /*lee_datos_float*/ {
/***************************************************************************/
  FILE *canal;
  char nombre[MAXNAMELENGTH];
  float val;
  int ix,iy;

  sprintf(nombre,"%s.fdat",nombre_in);
  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );

  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<leff;ix++)	{ 
      fread(&val,sizeof(float),1,canal);
      signal[iy][ix]=(double)val;
    }
  
  fclose(canal);

  return OK;
} // end of read_data_infloat


/***************************************************************************/
int write_data_infloat( int leff, int dimy, char *nombre_in, double **signal)
  /*graba_datos_float*/ {
/***************************************************************************/
  FILE *canal;
  char nombre[MAXNAMELENGTH];
  float val;
  int ix,iy;

  sprintf(nombre,"%s.fdat",nombre_in);
  TrackNull( canal=fopen(nombre,"wb"), "Error opening file in 'w' mode" );
  
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<leff;ix++)	{ 
      val=(float)signal[iy][ix];
      fwrite(&val,sizeof(float),1,canal);
    }
  
  fclose(canal);

  sprintf(nombre,"%s.hdr",nombre_in);
  canal=fopen(nombre,"wt");
  fprintf(canal,"%s.dat\n%d",nombre_in,leff);
  fclose(canal);

  return OK;
} // end of write_data_infloat


/***************************************************************************/
int write_serie( int leff, char *nombre_in, double *datos)/*graba_serie*/ {
/***************************************************************************/
  FILE *canal;
  char nombre[MAXNAMELENGTH];
  int ix;

  sprintf(nombre,"%s.txt",nombre_in);
  TrackNull( canal=fopen(nombre,"wt"), "Error opening file in 'w' mode" );
  
  for(ix=0;ix<leff;ix++)
    fprintf(canal,"%f\n",datos[ix]);
  fclose(canal);
  
  return OK;
} // end of write_serie
