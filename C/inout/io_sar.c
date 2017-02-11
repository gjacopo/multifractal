#include <stdio.h>
#include <math.h>

#include <utl_alloc.h>
#include <utl_stats.h>

#include <utl_operator.h>

#include <io_sar.h>


void read_sar( int dimx, int dimy, int x0, int y0, char *nombre, double ***data)
{
     FILE *canal;
     double norm;
     char val;
     int value;
     int bx,by;
     int ix,iy;
     int lx,ly;

     fill0(dimx,dimy,data[0],NULL);
     fill0(dimx,dimy,data[1],NULL);

     canal=fopen(nombre,"rb");
     fseek(canal,Y0SAR*(X0SAR+2*XMAXSAR),SEEK_SET);
     if((x0<0)||(y0<0))
     {
           bx=XMAXSAR/dimx;
	   by=YMAXSAR/dimy;
	   norm=(double)(bx*by);
           for(iy=0;iy<dimy;iy++)
	   {
	     for(ly=0;ly<by;ly++)
	     {
	            fseek(canal,X0SAR,SEEK_CUR);
		    for(ix=0;ix<dimx;ix++)
		    {
		      for(lx=0;lx<bx;lx++)
		      {
			fread(&val,sizeof(char),1,canal);
			value=(int) val;
			if(value<0) value+=256;
			data[0][iy][ix]+=((double)value)/norm;
			fread(&val,sizeof(char),1,canal);
			value=(int) val;
			if(value<0) value+=256;
			data[1][iy][ix]+=((double)value)/norm;			
		      }
		    }
		    fseek(canal,2*(XMAXSAR-dimx*bx),SEEK_CUR);
	     }
	   }	   
     }
     else
     {
           fseek(canal,y0*(X0SAR+2*XMAXSAR),SEEK_CUR);
	   for(iy=0;iy<dimy;iy++)
	   {
	          fseek(canal,X0SAR+2*x0,SEEK_CUR);
		  for(ix=0;ix<dimx;ix++)
		  {
		    fread(&val,1,1,canal);
		    value=(int) val;
		    if(value<0) value+=256;
		    data[0][iy][ix]=(double)value;
		    fread(&val,1,1,canal);
		    value=(int) val;
		    if(value<0) value+=256;
		    data[1][iy][ix]=(double)value;
		  }
		  fseek(canal,2*(XMAXSAR-x0-dimx),SEEK_CUR);
	   }
     }
     fclose(canal);

}


void write_sar( int dimx, int dimy, int x0, int y0, char *nombre, double ***data)
{
     FILE *canal;
     double norm;
     char val[2];
     int value;
     int bx,by;
     int ix,iy;
     int lx,ly;

     canal=fopen(nombre,"r+b");
     fseek(canal,(Y0SAR+y0)*(X0SAR+2*XMAXSAR),SEEK_SET);
     for(iy=0;iy<dimy;iy++)
     {
       fseek(canal,X0SAR+2*x0,SEEK_CUR);
       for(ix=0;ix<dimx;ix++)
       {
	 val[0]=(char)(int)data[0][iy][ix];
	 val[1]=(char)(int)data[1][iy][ix];
	 fwrite(&(val[0]),1,2,canal);
       }
       fseek(canal,2*(XMAXSAR-x0-dimx),SEEK_CUR);
     }
     fclose(canal);

}
