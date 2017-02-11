
/*	Ficheros de manipulacion de formatos graficos	*/
/*	Version del 3 de Julio de 2002		*/

/*	Necesita de las funciones de reserva de memoria	*/
/*	Esta disennado para matrices de datos double	*/


/*	Parametros de dimension (formato VH)		*/

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
#include <io_grafic.h>  	


#ifdef _PARSE_IO_PARAMETERS_
#include <io_parse.h>
extern ParIO *p_io;
#endif /*!_PARSE_IO_PARAMETERS_*/

/*		Parametros GIF			*/
#ifndef DELAY
#define DELAY 50 // centesimas de segundo de retardo entre frames
#endif
#ifndef WR
#define WR 0.3
#endif
#ifndef WG
#define WG 0.59 // pesos de las componentes cromaticas en las calibraciones
#endif
#ifndef WB
#define WB 0.11
#endif


/***************************************************************************/
int ext_grafic( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		char *nombre, Read2D *p_lee ) {
  /***************************************************************************/
  int error=ERROR;
  
  if(error) error=ext_gif(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
  if(error) error=ext_ppm(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
  // if(error) error=ext_unformatted(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
  
  return error;
} // end of ext_grafic


/***************************************************************************/
int check_grafic( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		  char *nombre, Read2D *p_lee) {
/***************************************************************************/
  int error=ERROR;

  if(error) error=check_gif(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
  if(error) error=check_ppm(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
  // if(error) error=check_unformatted(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
    
  return error;
} // end of check_grafic


/***************************************************************************/
int check_gif( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	       char *nombre,  Read2D *p_lee) {
  /***************************************************************************
   * Called by check_grafic
   *           ext_gif
   ***************************************************************************/

  FILE *canal;
  char R,G,B;
  unsigned char buffer[MAXCHARBUFFLENGTH];
  unsigned char quest;
  unsigned char b0;
  int error;
  int gdimv;
  int xmax,ymax;
  int x0,y0, ix,iy;
  int ic,ie,id,id0,it,in;
  int lt,last_cd;
  unsigned int color_res,color_size;
  int bits,bit0,w0,w1,l_code,letra;
  int verbose=0; // For debugging purposes;
  int global_ct,local_ct,fin,apren;
	
  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );
 
  fread(buffer,sizeof(char),6,canal);
  if(feof(canal)) return ERROR;

  buffer[6]='\0';
  if(verbose) fprintf(stderr,"Version: %s\n",buffer);
  fread(buffer,sizeof(char),7,canal);
  if(feof(canal)) return ERROR;
   
  *dimx=256*(int)buffer[1]+buffer[0];
  *dimy=256*(int)buffer[3]+buffer[2];
  gdimv=*dimv=1; // unless we see color we asume that images are monochrome

  if(verbose) fprintf(stderr,"Dimensiones: %d x %d\n",*dimx,*dimy);
  if((buffer[4]&0x80)&&verbose) 
    fprintf(stderr,"Color de fondo: %d\n",(int)buffer[5]);
  color_res=1+(int)(buffer[4]&0x70)/16;
  if(verbose) fprintf(stderr,"Profundidad de color: %d bits\n",color_res);
  if(buffer[4]&0x80)   {
    global_ct=1;
    color_size=1+(int)(buffer[4]&0x07);
    color_size=(int)pow(2.,(double)color_size);
    if((buffer[4]&0x2)&&verbose)
      fprintf(stderr,"Tabla de color ordenada\n");
    if(verbose)	{
      fprintf(stderr,"Tamanno de la tabla de color: %d\n",
	      color_size);
      fprintf(stderr,"Razon de aspecto del pixel: %f\n",
	      ((double)buffer[6]+15.)/64);
    }
  } else    {
    global_ct=0;
    if(verbose) fprintf(stderr,"Imagen Gif sin tabla de color global\n");
  }
  if(global_ct==1)    {
    for(ic=0;ic<color_size;ic++)	{
      R=getc(canal);
      G=getc(canal);
      B=getc(canal);
      if((R!=G)||(R!=B)||(G!=B)) gdimv=3;
    }
  }
  
  for(quest=getc(canal),in=0;quest!=0x3B;quest=getc(canal))    {
    if(quest==0x21) interpreta_extension_gif(canal,verbose);
    else if(quest!=0x2c)      {
      ErrorV("Identificador desconocido!! Codigo: %d\n",
	       (int)quest);
    } else      {
      fread(buffer,sizeof(char),9,canal);
      if(feof(canal)) return(-1);
      x0=256*(int)buffer[1]+(int)buffer[0];
      y0=256*(int)buffer[3]+(int)buffer[2];
      xmax=256*(int)buffer[5]+(int)buffer[4];
      ymax=256*(int)buffer[7]+(int)buffer[6];
      local_ct=(int)(buffer[8]&0x80);
      if(local_ct)	{
	color_size=1+(int)(buffer[8]&0x07);
	color_size=(int)pow(2.,(double)color_size);
      }
      else *dimv=Max(*dimv,gdimv);
      if(verbose)	{
	fprintf(stderr,"Imagen numero: %d\n",in);
	fprintf(stderr,"Posicion de la ventana: (%d,%d)\n",
	       x0,y0);
	fprintf(stderr,"Dimensiones: %d x %d\n",xmax,ymax);
	if(local_ct)	  {
	  fprintf(stderr,"Tamanno de la tabla de color local: %d\n",
		 color_size);
	  if(buffer[8]&0x20) fprintf(stderr,"Tabla ordenada\n");
	}
	if(buffer[8]&0x40) fprintf(stderr,"Imagen entrelazada\n");
      }
      if(local_ct)	    {
	for(ic=0;ic<color_size;ic++)	  {
	  R=getc(canal);
	  G=getc(canal);
	  B=getc(canal);
	  if((R!=G)||(R!=B)||(G!=B)) *dimv=3;
	}
      }
      
      if((local_ct==0)&&(global_ct==0))	{
	Error("Image without table of colors: decodificacion impossible");
      }
      
      /*		lectura y decodificacion de la imagen		*/      
      bit0=(int)getc(canal);
      error=read_comment_gif(canal,0);
      in++;
    }
  }
  fclose(canal);
  
  if(!error)    {
    if(verbose) printf("Numero de imagenes: %d\n",in);
    
    *dimz=in;
    if(*dimv==3) *bd=3;
    else *bd=1;
    
    *p_lee=&read_gif;
  }
  
  return error;
} // end of check_gif


/***************************************************************************/
int check_ppm( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	       char *nombre, Read2D *p_lee ) {
  /***************************************************************************
   * Called by check_grafic
   *           ext_ppm
   ***************************************************************************/
  FILE* canal;
  char iden[3],buffer[MAXCHARBUFFLENGTH];
  double value;
  int error=0;
  int levels;
  int ix,iy;

  *dimx=*dimy=*dimz=*dimv=*bd=0;

  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode");
  fread(iden,sizeof(char),3,canal);
  if(feof(canal)) return ERROR;

  if((iden[0]!='P')||
     ((iden[1]!='2')&&(iden[1]!='3')&&(iden[1]!='5')&&(iden[1]!='6'))
     ) return ERROR;
    
  *dimx=read_en_ascii(canal);  
  *dimy=read_en_ascii(canal);
  if((iden[1]=='3')||(iden[1]=='6')) *dimv=3;
  else *dimv=1;
  *dimz=1;
  levels=1+read_en_ascii(canal);
  *bd=(int)(log((double)levels)/log(256.));

  error=OK;
  switch(iden[1])    {
    
  case '2': // ASCII graylevel
    for(iy=0;(iy<*dimy)&&(!error);iy++)
      for(ix=0;(ix<*dimx)&&(!error);ix++)	    {
	read_en_ascii(canal);
	error=feof(canal);
      }
    break;
  case '3': // ASCII color
    for(iy=0;(iy<*dimy)&&(!error);iy++)
      for(ix=0;(ix<*dimx)&&(!error);ix++)	  {
	read_en_ascii(canal);
	read_en_ascii(canal);
	read_en_ascii(canal);
	error=feof(canal);
      }
    break;
  case '5': // Binary graylevel
    for(iy=0;(iy<*dimy)&&(!error);iy++)
      for(ix=0;(ix<*dimx)&&(!error);ix++) {
	error=(fscanf(canal,"%c",&buffer)==EOF)?ERROR:OK;
      }
    break;
  case '6': // Binary color
    for(iy=0;(iy<*dimy)&&(!error);iy++)
      for(ix=0;(ix<*dimx)&&(!error);ix++)	  {
	error=(fscanf(canal,"%c",&buffer)==EOF)?ERROR:OK;
	error=(fscanf(canal,"%c",&buffer)==EOF)?ERROR:OK;
	error=(fscanf(canal,"%c",&buffer)==EOF)?ERROR:OK;
      }
    break;
  }
  
  if((ix==*dimx)&&(iy==*dimy))    {	
    if(fscanf(canal,"%c",buffer)==EOF) error=OK; // File is over, ok!
    else error=ERROR;
  }
  else error=ERROR;
  fclose(canal);

  if(!error) *p_lee=&read_ppm;

  return error;
} // end of check_ppm


/***************************************************************************/
int read_gif( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
	      char *nombre, double ***data ) {
  /***************************************************************************
   * Called by check_gif
   ***************************************************************************/

  FILE *canal;
  double value;
  unsigned char buffer[MAXCHARBUFFLENGTH];
  unsigned char quest;
  unsigned char *Rg,*Gg,*Bg;
  unsigned char *Rl,*Gl,*Bl;
  unsigned char *encoded;
  unsigned char b0;
  int *decoded,*tab0,*tab1;
  int xmax,ymax;
  int x0,y0;
  int ibe,base;
  int ic,ie,id,id0,it,in;
  int lt,last_cd;
  int ix,iy;
  unsigned int color_res,color_size;
  int bits,bit0,w0,w1,l_code,letra;
  int verbose=0,global_ct,local_ct,fin,apren;

  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );
  fread(buffer,sizeof(char),6,canal);

  buffer[6]='\0';
  if(verbose) fprintf(stderr,"Version: %s\n",buffer);

  fread(buffer,sizeof(char),7,canal);
  xmax=256*(int)buffer[1]+buffer[0];
  ymax=256*(int)buffer[3]+buffer[2];
  if(verbose) fprintf(stderr,"Dimensiones: %d x %d\n",xmax,ymax);
  if((xmax>dimx)||(ymax>dimy))   Error("Image too big");

  if((buffer[4]&0x80)&&verbose) 
    fprintf(stderr,"Color de fondo: %d\n", (int)buffer[5]);
  color_res=1+(int)(buffer[4]&0x70)/16;
  if(verbose) fprintf(stderr,"Profundidad de color: %d bits\n",color_res);

  if(buffer[4]&0x80)    {
    global_ct=1;
    color_size=1+(int)(buffer[4]&0x07);
    color_size=(int)pow(2.,(double)color_size);
    if((buffer[4]&0x2)&&verbose)
      fprintf(stderr,"Tabla de color ordenada\n");
    if(verbose)	{
      fprintf(stderr,"Tamanno de la tabla de color: %d\n",
	      color_size);
      fprintf(stderr,"Razon de aspecto del pixel: %f\n",
	      ((double)buffer[6]+15.)/64);
    }

  }  else    {
    global_ct=0;
    if(verbose) Warning("Imagen Gif sin tabla de color global\n");
  }
  
  if(global_ct==1)    {
      Rg=(char *) calloc(color_size,sizeof(unsigned char));
      Gg=(char *) calloc(color_size,sizeof(unsigned char));
      Bg=(char *) calloc(color_size,sizeof(unsigned char));

      for(ic=0;ic<color_size;ic++)	{
	Rg[ic]=getc(canal);
	Gg[ic]=getc(canal);
	Bg[ic]=getc(canal);
      }
  }
  
  for(quest=getc(canal),in=0;quest!=0x3B;quest=getc(canal))    {
    if(quest==0x21) interpreta_extension_gif(canal,verbose);
    
    else if(quest!=0x2c)	{
      ErrorV("Unknown identifier: code %d", (int)quest);
    }
    
    else	{
      fread(buffer,sizeof(char),9,canal);
      x0=256*(int)buffer[1]+(int)buffer[0];
      y0=256*(int)buffer[3]+(int)buffer[2];
      xmax=256*(int)buffer[5]+(int)buffer[4];
      ymax=256*(int)buffer[7]+(int)buffer[6];
      local_ct=(int)(buffer[8]&0x80);
      if(local_ct)	{
	color_size=1+(int)(buffer[8]&0x07);
	color_size=(int)pow(2.,(double)color_size);
      }
      if(verbose)	{
	fprintf(stderr,"Imagen numero: %d\n",in);
	fprintf(stderr,"Posicion de la ventana: (%d,%d)\n",x0,y0);
	fprintf(stderr,"Dimensiones: %d x %d\n",xmax,ymax);
	if(local_ct) {
	  fprintf(stderr,"Tamanno de la tabla de color local: %d\n",
		  color_size);
	  if(buffer[8]&0x20) fprintf(stderr,"Tabla ordenada\n");
	}
	if(buffer[8]&0x40) fprintf(stderr,"Imagen entrelazada\n");
      }
      if(local_ct)	    {
	Rl=(char *) calloc(color_size,sizeof(unsigned char));
	Gl=(char *) calloc(color_size,sizeof(unsigned char));
	Bl=(char *) calloc(color_size,sizeof(unsigned char));
	
	for(ic=0;ic<color_size;ic++)	  {
	  Rl[ic]=getc(canal);
	  Gl[ic]=getc(canal);
	  Bl[ic]=getc(canal);
	}
      }
      
      if((local_ct==0)&&(global_ct==0)&&(in==iz))
	Error("Image without table of colors: decodificacion impossible");
      
      /*  lectura y decodificacion de la imagen */
      
      bit0=(int)getc(canal);
      if(in==iz)	{
	if(verbose)	  {
	  fprintf(stderr,"\nEmpieza la decodificacion...\n");
	  fprintf(stderr,"Tamanno inicial de palabra: %d bits\n",bit0);
	}
	l_code=interpreta_stream_gif(canal,&encoded);
	if(verbose) fprintf(stderr,"Tamanno del codigo: %d bytes\n",
			   l_code);
	decoded=(int *) calloc(xmax*ymax,sizeof(int));
	tab0=(int *) calloc(4096,sizeof(int));
	tab1=(int *) calloc(4096,sizeof(int));
		
	w0=(int) pow(2.,(double)bit0);
	w1=2*w0;
	bits=bit0+1;
	
	lt=0;
	last_cd=w0;
	id0=0;
	for( id=0,ie=0,b0=0x01,fin=1; fin==1; )	  {
	  if(id>=xmax*ymax) fin=3;
	  /* Interpretando el caracter vigente  */	  
	  
	  base=1;
	  letra=0;
	  for(ibe=0;(ibe<bits)&&(fin!=2);ibe++)	      {
	    if(encoded[ie]&b0) letra+=base;
	    base=2*base;
	    b0=2*b0;
	    if(b0==0x00)	      {
	      b0=0x01;
	      ie++;
	      if(ie>=l_code) fin=2;
	    }
	  }
	  
	  /* Produciendo la salida asociada  */
	  
	  if(letra==w0)		    {
	    lt=0;
	    bits=bit0+1;
	    w1=2*w0;
	  }
	  else if(letra==w0+1) fin=0;
	  else if(fin==1)	      {
	    if(letra>lt+w0+2) 
	      fin=4;
	    else if(letra==lt+w0+2) {
	      /*  Aprendizaje de un codigo no visto */
	      
	      apren=1;
	      tab0[lt]=last_cd;
	      tab1[lt]=primer_car_gif(w0,last_cd,tab0);
	      lt++;
	    }

	    /* Reproduccion del codigo	*/	    
	    id=annade_decoded_gif(w0,letra,tab0,tab1,id,decoded);	    
	  }

	  /* Aprendizaje de codigos vistos */
	  if((letra==w0)||(letra==w0+1)||(fin!=1)) apren=1;
	  if((apren!=1)&&(last_cd!=w0))	    {
	    
	    tab0[lt]=last_cd;
	    tab1[lt]=decoded[id0];
	    
	    for(ic=0;(ic<lt)&&(id0<id);ic++)	      {
	      
	      if((tab0[ic]==tab0[lt])&&(tab1[ic]==tab1[lt]))		{
		tab0[lt]=ic+w0+2;
		id0++;
		if(id0<id) tab1[lt]=decoded[id0];
	      }
	    }
	    if(id0<id) lt++;
	  }	  
	  
	  /*		Actualizacion de variables		*/	  
	  if((lt+w0+2==w1)&&(letra!=w0))		    {
	    w1=2*w1;
	    bits++;
	    if(bits>12) bits=12;
	  }
	  last_cd=letra;
	  apren=0;
	  id0=id;
	  
	}
	if((fin==0)&&verbose) fprintf(stderr,"Lectura exitosa!\n");
	if(fin==2) 
	  fprintf(stderr,"Error: sennal fin de imagen no encontrada\n");
	if(fin==3) fprintf(stderr,"Error: codigo desborda imagen\n");
	if(fin==4) fprintf(stderr,"Error de lectura\n");
	
	/*		Fin de la decodificacion		*/	
	free(encoded);
	free(tab0);
	free(tab1);
	
	for(id=0,ix=0,iy=0;id<xmax*ymax;id++,ix++)		{
	  if(ix >= xmax)	    {
	    ix-=xmax;
	    iy++;
	  }
	  if(dimv==3)		    {
	    if(local_ct)	      {
	      data[0][dimy-1-iy][ix]=(double)Rl[decoded[id]];
	      data[1][dimy-1-iy][ix]=(double)Gl[decoded[id]];
	      data[2][dimy-1-iy][ix]=(double)Bl[decoded[id]];
	    }  else   {
	      data[0][dimy-1-iy][ix]=(double)Rg[decoded[id]];
	      data[1][dimy-1-iy][ix]=(double)Gg[decoded[id]];
	      data[2][dimy-1-iy][ix]=(double)Bg[decoded[id]];
	    }
	  }  else  {
	    if(local_ct) 
	      data[0][dimy-1-iy][ix]=
		WR*(double)Rl[decoded[id]]+
		WG*(double)Gl[decoded[id]]+
		WB*(double)Bl[decoded[id]];
	    else
	      data[0][dimy-1-iy][ix]=
		WR*(double)Rg[decoded[id]]+
		WG*(double)Gg[decoded[id]]+
		WB*(double)Bg[decoded[id]];
	  }
	}
	free(decoded);
      }
      else read_comment_gif(canal,0);
      in++;
      
      if(local_ct)	{
	free(Rl);
	free(Gl);
	free(Bl);
      }
      
    }
  }
  fclose(canal);
  if(global_ct)    {
    free(Rg);
    free(Gg);
    free(Bg);
  }
  
  if(verbose) fprintf(stderr,"Numero de imagenes: %d\n",in);
  
  if(in<iz)    {
    printf("Imagen %d de %d no encontrada\n",iz,in);
    return(-1);
  }
  else return(color_size);
} // end of read_gif


/***************************************************************************/
int read_ppm( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
	      char *nombre, double ***data ) {
/***************************************************************************/

  FILE* canal;
  char iden[3];
  double value;
  int xmax,ymax,levels;
  int ix,iy;

  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );
  fread(iden,sizeof(char),3,canal);
  if(feof(canal)) return ERROR;

  xmax=read_en_ascii(canal);  
  ymax=read_en_ascii(canal);
  levels=1+read_en_ascii(canal);

  switch((int)iden[1])    {

    case '2': // ASCII graylevel
      for(iy=0;iy<dimy;iy++)
	for(ix=0;ix<dimx;ix++)
	  data[0][dimy-1-iy][ix]=(double) read_en_ascii(canal);
      break;
  case '3': // ASCII color
    for(iy=0;iy<ymax;iy++)
      for(ix=0;ix<xmax;ix++)	{
	data[0][dimy-1-iy][ix]=(double) read_en_ascii(canal);
	data[1][dimy-1-iy][ix]=(double) read_en_ascii(canal);
	data[2][dimy-1-iy][ix]=(double) read_en_ascii(canal);
      }
    break;
  case '5': // Binary graylevel
    for(iy=0;iy<dimy;iy++)  
      for(ix=0;ix<dimx;ix++)	{
	value=(double)fgetc(canal);
	if(value<0) value+=256.;
	data[0][dimy-1-iy][ix]=value;
      }
    break;
  case '6': // Binary color
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)	    {
	value=(double)fgetc(canal);
	if(value<0) value+=256.;
	data[0][dimy-1-iy][ix]=value;
	value=(double)fgetc(canal);
	if(value<0) value+=256.;
	data[1][dimy-1-iy][ix]=value;
	value=(double)fgetc(canal);
	if(value<0) value+=256.;
	data[2][dimy-1-iy][ix]=value;
      }
    break;
  }
  
  fclose(canal);

  return levels;
} // end of read_ppm


/***************************************************************************/
int read_dimensiones_foto( char *nombre, int *dimx, int *dimy ) {
/***************************************************************************/
  FILE *canal;
  char buffer[3];
  int res;

  TrackNullAlloc( canal=fopen(nombre,"rt") );
  fread(buffer,sizeof(char),3,canal);
  fclose(canal);

  if((buffer[0]=='G')&&(buffer[1]=='I')&&(buffer[2]=='F'))
    res = read_dimensiones_gif( nombre, dimx, dimy );
  else if(buffer[0]=='P')
    res = read_dimensiones_ppm( nombre, dimx, dimy );
  else {
    WarnMsg("Unknown format. Impossible to provide dimensions");
    res = ERROR;
  }

  return res;
} // end of read_dimensiones_foto


/***************************************************************************/
int read_dimensiones_ppm( char *nombre, int *dimx, int *dimy)  {
/***************************************************************************/
  FILE* canal;
  char iden[3];

  TrackNull(canal=fopen(nombre,"rt"),"Error opening file in 'r' mode");
  fread(iden,sizeof(char),3,canal);
  if((iden[0]!='P')||((iden[1]!='2')&&(iden[1]!='3')&&
		      (iden[1]!='5')&&(iden[1]!='6')))
      return ERROR;

  *dimx=read_en_ascii( canal );  
  *dimy=read_en_ascii( canal );

  return OK;
} // end of read_dimensiones_ppm


/***************************************************************************/
int read_dimensiones_gif( char *nombre, int *dimx, int *dimy ) {
  /***************************************************************************/
  FILE *canal;
  unsigned char buffer[MAXCHARBUFFLENGTH];
  
  TrackNull(canal=fopen(nombre,"rb"),"Error opening file in 'r' mode");

  fread(buffer,sizeof(char),6,canal);
  fread(buffer,sizeof(char),7,canal);
  *dimx=256*(int)buffer[1]+buffer[0];
  *dimy=256*(int)buffer[3]+buffer[2];
  
  return OK;
} // end of read_dimensiones_gif


/***************************************************************************/
int ext_gif( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	     char *nombre, Read2D *p_lee ) {
  /***************************************************************************/
  char ext[90];
  int error=OK;
  int lext=extract_extension(nombre,'.',ext);

  if(lext!=3) error=ERROR;
  else if(
	  ((ext[0]!='g')&&(ext[0]!='G'))||
	  ((ext[1]!='i')&&(ext[1]!='I'))||
	  ((ext[2]!='f')&&(ext[2]!='F'))
	  ) error=ERROR;
  
  if(!error) error=check_gif(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
  
  return error;
} // end of ext_gif


/***************************************************************************/
int ext_ppm( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	     char *nombre, Read2D *p_lee ) {
  /***************************************************************************/
  char ext[MAXCHARLENGTH];
  int error=OK;
  int lext=extract_extension(nombre,'.',ext);

  if(lext!=3) error=ERROR;
  else if(
	  ((ext[0]!='p')&&(ext[0]!='P'))||
	  ((ext[1]!='p')&&(ext[1]!='P')&&(ext[1]!='g')&&(ext[1]!='G'))||
	  ((ext[2]!='m')&&(ext[2]!='M'))
	  ) error=ERROR;
 
  if(!error) error=check_ppm(dimx,dimy,dimv,dimz,bd,nombre,p_lee);

  return error;
} // end of ext_ppm


/***************************************************************************/
int read_foto_gris( int dimx, int dimy, char *nombre, double **data ) {
  /***************************************************************************/
  FILE *canal_test;
  char buffer[3];
  int res;
  
  TrackNull( canal_test=fopen(nombre,"rt"), "Error opening file in 'r' mode" );
  fread(buffer,sizeof(char),3,canal_test);
  fclose(canal_test);
  
  if((buffer[0]=='G') && (buffer[1]=='I') && (buffer[2]=='F'))
    res = read_gif_gris( dimx, dimy, nombre, data );
  else if(buffer[0]=='P') {	
    switch(buffer[1])      {
    case '6':
    case '3':
      res=read_ppm_gris( dimx, dimy, nombre, data );
      break;
    case '5':
    case '2':
      res=read_pgm( dimx, dimy, nombre, data );
      break;
    default : res=ERROR;
    }
  } else res=read_ppm_raw( dimx, dimy, nombre, data );
  
  if(res == ERROR) Exit("Unknown format");

  return res;
} // end of read_foto_gris


/***************************************************************************/
int read_foto_color( int dimx, int dimy, char *nombre, char **Red, 
		     char **Green, char **Blue ) {
  /***************************************************************************/

  FILE *canal_test;
  char buffer[3];
  double **aux;
  int ix,iy,res;
  
  TrackNull( canal_test=fopen(nombre,"rt"), "Error opening file in 'r' mode" );
  fread(buffer,sizeof(char),3,canal_test);
  fclose(canal_test);
  
  if((buffer[0]=='G') && (buffer[1]=='I') && (buffer[2]=='F'))
    res = read_gif_rgb( dimx, dimy, 0, nombre, Red, Green, Blue );
  else if(buffer[0] == 'P')    {	
    switch(buffer[1])	{
    case '6':
    case '3':
      res = read_ppm_rgb( dimx, dimy, nombre, Red, Green, Blue );
      break;
    case '5':
    case '2':
      aux = matrix2D( dimy, dimx ); 
      res = read_pgm( dimx, dimy, nombre, aux );
      for( iy=0; iy<dimy; iy++ )  
	for( ix=0; ix<dimx; ix++ )
	  Red[iy][ix] = (char)(int)aux[iy][ix];
      ccopy( dimx, dimy, Red, Green,NULL );
      ccopy( dimx, dimy, Red, Blue,NULL );
      free_matrix2D( aux, dimy );
      break;
    default : res = ERROR;
    }
  }
  else res = ERROR;
  
  if(res == ERROR) Exit("Unknown format");

  return res;
} // end of read_foto_color


/***************************************************************************/
int read_color_block( int dimx, int dimy, int block, char *nombre,
		      double ***data, double *med_cr) {
/***************************************************************************/
  char ***color;
  int sizex,sizey;
  int ic,ix,iy,ibx,iby;
  int dcr;

  TrackNullAlloc( color=cmatrix3D(3,dimy,dimx) );

  read_foto_color(dimx,dimy,nombre,color[0],color[1],color[2]);

  sizex=dimx/block, sizey=dimy/block;

  for(iy=0;iy<sizey;iy++)    
    for(ix=0;ix<sizex;ix++)	{
      for(iby=0;iby<block;iby++)	    
	for(ibx=0;ibx<block;ibx++)	  
	  for(ic=0;ic<3;ic++)	    {
	    dcr=(int)color[ic][iy*block+iby][ix*block+ibx];
	    if(dcr<0) dcr+=256;
	    data[ic][iy][ix]+=dcr;
	  }
      for(ic=0;ic<3;ic++) data[ic][iy][ix]=
			    data[ic][iy][ix]/((double)(block*block));
    }
  
  for(ic=0;ic<3;ic++)    {
    med_cr[ic]=0.;
    for(iy=0;iy<dimy;iy++)     	
      for(ix=0;ix<dimx;ix++)   
	med_cr[ic]+=data[ic][iy][ix];
    med_cr[ic]=med_cr[ic]/((double)(dimx*dimy));
    for(iy=0;iy<dimy;iy++)	
      for(ix=0;ix<dimx;ix++)
	  data[ic][iy][ix]-=med_cr[ic];
  }
  
  free_cmatrix3D(color,3,dimy);
  
  return(256);
} // end of read_color_block


/***************************************************************************/
int read_ppm_gris( int dimx, int dimy, char *nombre, double **data) {
  /***************************************************************************/
  char **Red,**Green,**Blue;
  int ix,iy;
  int dR,dG,dB;
  
  TrackNullAlloc( Red=cmatrix2D(dimy,dimx) );
  TrackNullAlloc( Green=cmatrix2D(dimy,dimx) );
  TrackNullAlloc( Blue=cmatrix2D(dimy,dimx) );

  read_ppm_rgb(dimx,dimy,nombre,Red,Green,Blue);
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      dR=(int)Red[iy][ix];
      if(dR<0) dR+=256;
      dG=(int)Green[iy][ix];
      if(dG<0) dG+=256;
      dB=(int)Blue[iy][ix];
      if(dB<0) dB+=256;
      data[iy][ix]=WR*(double)dR+WG*(double)dG+WB*(double)dB;
    }

  free_cmatrix2D(Red,dimy);
  free_cmatrix2D(Green,dimy);
  free_cmatrix2D(Blue,dimy);
  
  return(256);
} // end of read_ppm_gris


/***************************************************************************/
int read_ppm_raw( int dimx, int dimy, char* nombre, double **data) {
  /***************************************************************************/
  int lx,ly;
  int ix0,iy0;
  FILE* canal;
  double value;
  int dat1,dat2,ix,iy;
  
  ix0=(XMAXVH-dimx)/2, iy0=(YMAXVH-dimy)/2;

  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );

  if(ix0<0)    {
    ix0=0,    lx=XMAXVH;
  }  else lx=dimx;

  if(iy0<0)    {
    iy0=0,    ly=YMAXVH;
  }  else ly=dimy;
  
  for( iy=0; iy<YMAXVH; iy++ )    
    for( ix=0; ix<XMAXVH; ix++ )	{
      dat1=(int) getc(canal);
      if(dat1<0) dat1=dat1+256;
      dat2=(int) getc(canal);
      if(dat2<0) dat2=dat2+256;      
      value=(double)(dat2+256*dat1);      
      if((ix>=ix0)&&(ix<ix0+lx)&&(iy>=iy0)&&(iy<iy0+ly))
	data[ly-1-(iy-iy0)][ix-ix0]=value;
    }
  
  fclose(canal);

  return 32768;
} // end of read_ppm_raw


/***************************************************************************/
int read_ppm_rgb( int dimx, int dimy, char *nombre, char **Red, char **Green,
		  char **Blue)  {
/***************************************************************************/
  FILE* canal_test;
  char iden[3];
  double ajuste;
  int xmax,ymax,ncolor;
  int ix,iy;
  int dato,verbose=0;

  TrackNull( canal_test=fopen(nombre,"rt"), "Error opening file in 'r' mode" );
  fread(iden,sizeof(char),3,canal_test);
  if((iden[0]!='P')||((iden[1]!='3')&&(iden[1]!='6')))
      Error("Image not in PPM");

  xmax=read_en_ascii(canal_test);  
  ymax=read_en_ascii(canal_test);
  ncolor=1+read_en_ascii(canal_test);
  if(verbose)    {
    fprintf(stderr,"\n");
    fprintf(stderr,"Leyendo la imagen: %s\n",nombre);
    fprintf(stderr,"Dimensiones: %d x %d\n",xmax,ymax);
    fprintf(stderr,"Niveles por color: %d\n",ncolor);
  }
  if((xmax>dimx)||(ymax>dimy))
    Error("Image too big");
  

  if(iden[1]=='6')
    for(iy=0;iy<ymax;iy++)
      for(ix=0;ix<xmax;ix++)	{
	fread(&Red[iy][ix],sizeof(char),1,canal_test);
	fread(&Green[iy][ix],sizeof(char),1,canal_test);
	fread(&Blue[iy][ix],sizeof(char),1,canal_test);
      }
  else    {
    if(ncolor>256)      {
      ajuste=255/((double)ncolor-1.);
      ncolor=256;
    }
    
    for(iy=0;iy<ymax;iy++)
      for(ix=0;ix<xmax;ix++)	  {
	dato=read_en_ascii(canal_test);
	if(ncolor>256) dato=(int) (ajuste*(double)dato);
	Red[iy][ix]=(char) dato;
	
	dato=read_en_ascii(canal_test);
	if(ncolor>256) dato=(int) (ajuste*(double)dato);
	Green[iy][ix]=(char) dato;
	
	dato=read_en_ascii(canal_test);
	if(ncolor>256) dato=(int) (ajuste*(double)dato);
	Blue[iy][ix]=(char) dato;
      }
  }

  fclose(canal_test);

  return(ncolor);
} // end of read_ppm_rgb


/***************************************************************************/
int read_pgm( int dimx, int dimy, char *nombre_in, double **cont) {
  /***************************************************************************
   * Called by read_foto_gris
   *           read_foto_color
   ***************************************************************************/
  FILE* canal_test;
  char iden[3];
  int ix,iy;
  int xmax,ymax,ngris,dato;
  int verbose=0;

  TrackNull( canal_test=fopen(nombre_in,"rt"), "Error opening file in 'r' mode" );
  fread(iden,sizeof(char),3,canal_test);
  if((iden[0]!='P')||((iden[1]!='2')&&(iden[1]!='5')))
      Error("Image not in PGM");

  xmax=read_en_ascii(canal_test);
  ymax=read_en_ascii(canal_test);
  ngris=1+read_en_ascii(canal_test);

  if(verbose)    {
    fprintf(stderr,"\n");
    fprintf(stderr,"Leyendo la imagen: %s\n",nombre_in);
    fprintf(stderr,"Dimensiones: %d x %d\n",xmax,ymax);
    fprintf(stderr,"Niveles de gris: %d\n",ngris);
  }
  if((xmax>dimx)||(ymax>dimy))
    Error("Image too big");
  
  if(iden[1]=='5')    
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	  cont[iy][ix]=(double) getc(canal_test);
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) {
	dato=read_en_ascii(canal_test);
	cont[iy][ix]=(double) dato;
      }

  fclose(canal_test);

  return(ngris);
} // end of read_pgm


/***************************************************************************/
int read_en_ascii/*lee_en_ascii*/( FILE* canal ) {
  /***************************************************************************
   * Called by check_ppm
   *           read_ppm
   *           read_ppm_rgb
   *           read_pgm
   *           read_dimensiones_ppm 
   ***************************************************************************/
  int comienzo=0,fin=0,comentario=0;
  int leo;
  int salida=0;

  while( !fin )    {
      leo = getc(canal);
      if(leo == '#') comentario=1;
      if(leo == '\n') comentario=0;
      if( !comentario )	{
	if((leo>='0')&&(leo<='9')) 	  {
	  salida = salida*10 + (leo-'0');
	  comienzo = 1;
	}
	else if(comienzo) fin=1;
      }
  }
  
  return salida;
} // end of read_en_ascii


/***************************************************************************/
int read_gif_gris( int dimx, int dimy, char *nombre, double **data ) {
  /***************************************************************************
   * Called by read_foto_gris
   ***************************************************************************/
  char **Red,**Green,**Blue;
  int ix,iy;
  int dR,dG,dB;

  TrackNullAlloc( Red=cmatrix2D(dimy,dimx) );
  TrackNullAlloc( Green=cmatrix2D(dimy,dimx) );
  TrackNullAlloc( Blue=cmatrix2D(dimy,dimx) );

  read_gif_rgb(dimx,dimy,0,nombre,Red,Green,Blue);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)	{
      dR=(int)Red[iy][ix];
      if(dR<0) dR+=256;
      dG=(int)Green[iy][ix];
      if(dG<0) dG+=256;
      dB=(int)Blue[iy][ix];
      if(dB<0) dB+=256;
      data[iy][ix]=WR*(double)dR+WG*(double)dG+WB*(double)dB;
    }

  free_cmatrix2D(Red,dimy);
  free_cmatrix2D(Green,dimy);
  free_cmatrix2D(Blue,dimy);

  return 256;
} // end of read_gif_gris


/***************************************************************************/
int read_gif_rgb( int dimx, int dimy, int nimag, char *nombre, char **Red, 
		  char **Green, char **Blue ) {
  /***************************************************************************
   * Called by read_foto_color
   *           read_gif_gris
   ***************************************************************************/
  FILE *canal;
  unsigned char buffer[MAXCHARBUFFLENGTH];
  unsigned char quest;
  char *Rg,*Gg,*Bg;
  char *Rl,*Gl,*Bl;
  unsigned char *encoded;
  unsigned char b0;
  int *decoded,*tab0,*tab1;
  int xmax,ymax;
  int x0,y0;
  int ibe,base;
  int ic,ie,id,id0,it,in;
  int lt,last_cd;
  int ix,iy;
  unsigned int color_res,color_size;
  int bits,bit0,w0,w1,l_code,letra;
  int verbose=0,global_ct,local_ct,fin,apren;

  TrackNull( canal=fopen(nombre,"rb"), "Error opening file in 'r' mode" );
  fread(buffer,sizeof(char),6,canal);

  buffer[6]='\0';
  if(verbose) fprintf(stderr,"Version: %s\n",buffer);
  fread(buffer,sizeof(char),7,canal);
  xmax=256*(int)buffer[1]+buffer[0];
  ymax=256*(int)buffer[3]+buffer[2];
  if(verbose) fprintf(stderr,"Dimensiones: %d x %d\n",xmax,ymax);
  if((xmax>dimx)||(ymax>dimy))
    Error("Image too big");
  
  if((buffer[4]&0x80)&&verbose) fprintf(stderr,"Color de fondo: %d\n",
					(int)buffer[5]);
  color_res=1+(int)(buffer[4]&0x70)/16;
  if(verbose) fprintf(stderr,"Profundidad de color: %d bits\n",color_res);
  
  if(buffer[4]&0x80)    {
    global_ct=1;
    color_size=1+(int)(buffer[4]&0x07);
    color_size=(int)pow(2.,(double)color_size);
    if((buffer[4]&0x2)&&verbose)
      fprintf(stderr,"Tabla de color ordenada\n");
    if(verbose)      {
      fprintf(stderr,"Tamanno de la tabla de color: %d\n",
	      color_size);
      fprintf(stderr,"Razon de aspecto del pixel: %f\n",
	      ((double)buffer[6]+15.)/64);
    }
    
  }  else    {
    global_ct=0;
    if(verbose) Warning("Imagen Gif sin tabla de color global");
  }
  
  if(global_ct==1)    {
    TrackNullAlloc( Rg=(char*)calloc(color_size,sizeof(char)) );
    TrackNullAlloc( Gg=(char*)calloc(color_size,sizeof(char)) );
    TrackNullAlloc( Bg=(char*)calloc(color_size,sizeof(char)) );
    
    for(ic=0;ic<color_size;ic++)	{
      Rg[ic]=getc(canal);
      Gg[ic]=getc(canal);
      Bg[ic]=getc(canal);
    }
  }
  
  for(quest=getc(canal),in=0;quest!=0x3B;quest=getc(canal))    {

    if(quest==0x21) interpreta_extension_gif(canal,verbose);

    else if(quest!=0x2c) {
      ErrorV("Unknown identifier: code %d", (int)quest);

    } else	{
      
      fread(buffer,sizeof(char),9,canal);
      x0=256*(int)buffer[1]+(int)buffer[0];
      y0=256*(int)buffer[3]+(int)buffer[2];
      xmax=256*(int)buffer[5]+(int)buffer[4];
      ymax=256*(int)buffer[7]+(int)buffer[6];
      local_ct=(int)(buffer[8]&0x80);
      if(local_ct)	{
	color_size=1+(int)(buffer[8]&0x07);
	color_size=(int)pow(2.,(double)color_size);
      }
      
      if(verbose)	    {
	fprintf(stderr,"Imagen numero: %d\n",in);
	fprintf(stderr,"Posicion de la ventana: (%d,%d)\n",
		x0,y0);
	fprintf(stderr,"Dimensiones: %d x %d\n",xmax,ymax);
	if(local_ct)	  {
	  fprintf(stderr,"Tamanno de la tabla de color local: %d\n",
		  color_size);
	  if(buffer[8]&0x20) fprintf(stderr,"Tabla ordenada\n");
	}
	if(buffer[8]&0x40) fprintf(stderr,"Imagen entrelazada\n");
      }
      
      if(local_ct)	{
	Rl=(char *) calloc(color_size,sizeof(char));
	Gl=(char *) calloc(color_size,sizeof(char));
	Bl=(char *) calloc(color_size,sizeof(char));
	
	for(ic=0;ic<color_size;ic++)	  {
	  Rl[ic]=getc(canal);
	  Gl[ic]=getc(canal);
	  Bl[ic]=getc(canal);
	}
      }
      
      if((local_ct==0)&&(global_ct==0)&&(in==nimag))
	Error("Image without table of colors: decodificacion impossible");
      
      /*  lectura y decodificacion de la imagen	 */
      
      bit0=(int)getc(canal);
      if(in==nimag)	{
	if(verbose)	  {
	  fprintf(stderr,"\nEmpieza la decodificacion...\n");
	  fprintf(stderr,"Tamanno inicial de palabra: %d bits\n",bit0);
	}
	l_code=interpreta_stream_gif(canal,&encoded);
	if(verbose) fprintf(stderr,"Tamanno del codigo: %d bytes\n",
			    l_code);
	decoded=(int *) calloc(xmax*ymax,sizeof(int));
	tab0=(int *) calloc(4096,sizeof(int));
	tab1=(int *) calloc(4096,sizeof(int));
		
	w0=(int) pow(2.,(double)bit0);
	w1=2*w0;
	bits=bit0+1;	
	
	lt=0;
	last_cd=w0;
	id0=0;
	for(id=0,ie=0,b0=0x01,fin=1;fin==1;)		{
	  if(id>=xmax*ymax) fin=3;
	  /*   Interpretando el caracter vigente       	*/	  
	  
	  base=1;
	  letra=0;
	  for(ibe=0;(ibe<bits)&&(fin!=2);ibe++)	    {
	    if(encoded[ie]&b0) letra+=base;
	    base=2*base;
	    b0=2*b0;
	    if(b0==0x00)	      {
	      b0=0x01;
	      ie++;
	      if(ie>=l_code) fin=2;
	    }
	  }
	  
	  /*   Produciendo la salida asociada	       	*/
	  
	  if(letra==w0)	    {
	    lt=0;
	    bits=bit0+1;
	    w1=2*w0;
	  }
	  else if(letra==w0+1) fin=0;
	  else if(fin==1)	    {
	    if(letra>lt+w0+2)	   
	      fin=4;
	    else if(letra==lt+w0+2)	      {
	      /* Aprendizaje de un codigo no visto  */   
	      apren=1;
	      tab0[lt]=last_cd;
	      tab1[lt]=primer_car_gif(w0,
				      last_cd,tab0);
	      lt++;
	    }

	    /*		Reproduccion del codigo			*/	    
	    id=annade_decoded_gif(w0,letra,
				  tab0,tab1,id,decoded);
	    
	  }
	  
	  /*		Aprendizaje de codigos vistos		*/	  	  
	  if((letra==w0)||(letra==w0+1)||(fin!=1))
	    apren=1;

	  if((apren!=1)&&(last_cd!=w0))	    {	      
	    tab0[lt]=last_cd;
	    tab1[lt]=decoded[id0];
	    
	    for(ic=0;(ic<lt)&&(id0<id);ic++)			{
	      
	      if((tab0[ic]==tab0[lt]) && (tab1[ic]==tab1[lt])) {
		tab0[lt]=ic+w0+2;
		id0++;
		if(id0<id) tab1[lt]=decoded[id0];
	      }
	    }
	    if(id0<id) lt++;
	  }
	  
	  
	  /*		Actualizacion de variables		*/
	  
	  if((lt+w0+2==w1)&&(letra!=w0))	    {
	    w1=2*w1;
	    bits++;
	    if(bits>12) bits=12;
	  }
	  
	  last_cd=letra;
	  apren=0;
	  id0=id;
	  
	}
	if((fin==0)&&verbose)
	  fprintf(stderr,"Lectura exitosa!\n");
	if(fin==2) fprintf(stderr,"Error: sennal fin de imagen no encontrada\n");
	if(fin==3) fprintf(stderr,"Error: codigo desborda imagen\n");
	if(fin==4) fprintf(stderr,"Error de lectura\n");
	
	/*		Fin de la decodificacion		*/	
	free(encoded);
	free(tab0);
	free(tab1);
	
	for(id=0,ix=0,iy=0;id<xmax*ymax;id++,ix++)	  {
	  if(ix>=xmax)	    {
	    ix-=xmax;
	    iy++;
	  }
	  
	  if(local_ct)	    {
	    Red[iy][ix]=Rl[decoded[id]];
	    Green[iy][ix]=Gl[decoded[id]];
	    Blue[iy][ix]=Bl[decoded[id]];
	  } else {
	    Red[iy][ix]=Rg[decoded[id]];
	    Green[iy][ix]=Gg[decoded[id]];
	    Blue[iy][ix]=Bg[decoded[id]];
	  }
	}
	free(decoded);
      }
      else read_comment_gif(canal,0);
      in++;
      
      if(local_ct)	{
	free(Rl);
	free(Gl);
	free(Bl);
      }      
    }
  }

  fclose(canal);
  if(global_ct)    {
    free(Rg);
    free(Gg);
    free(Bg);
  }
  
  if(verbose) fprintf(stderr,"Numero de imagenes: %d\n",in);
  
  if(in<nimag) {
    Error("Image do not exist");
  } else return(color_size);
  
} // end of read_gif_rgb


/***************************************************************************/
void interpreta_extension_gif( FILE *canal, int verbose) {
  /***************************************************************************
   * Called by check_gif
   *           read_gif
   *           read_gif_rgb
  ***************************************************************************/
  unsigned char decide;

  decide=getc(canal);

  switch(decide)    {
  case 0xF9:
    read_graphic_ce_gif(canal,verbose);
    break;
  case 0xFE:
    read_comment_gif(canal,verbose);
    break;
  case 0xFF:
    read_application_gif(canal,verbose);
    break;
  default:
    read_comment_gif(canal,0);
  }

} // end of interpreta_extension_gif


/***************************************************************************/
int read_graphic_ce_gif( FILE *canal, int verbose) {
  /***************************************************************************
   * Called by interpreta_extension_gif 
   ***************************************************************************/
  unsigned char buffer[6];
  unsigned int disposal,delay;

  fread(buffer,sizeof(char),6,canal);

  if(verbose)    {
    fprintf(stderr,"\nGraphical control extension parameters:\n");
    fprintf(stderr,"=======================================\n");
    disposal=(int) (buffer[1]&0x70)/16;
    switch(disposal)	{
    case 0:
      fprintf(stderr,"Disposicion: ninguna\n");
      break;
    case 1:
      fprintf(stderr,"Disposicion: ultimo\n");
      break;
    case 2:
      fprintf(stderr,"Disposicion: fondo\n");
      break;
    case 3:
      fprintf(stderr,"Disposicion: previo\n");
      break;
    default:
      fprintf(stderr,"Disposicion: indefinida\n");
    }
    if(buffer[1]&0x02) fprintf(stderr,"Interaccion con el usuario\n");
    if(buffer[1]&0x01)
      fprintf(stderr,"Color transparente: %d\n",(int) buffer[4]);
    delay=10*(256*(int) buffer[3]+(int)buffer[2]);
    fprintf(stderr,"Tiempo de retardo: %d ms\n",delay);
  }
  
  return OK;
} // end of read_graphic_ce_gif


/***************************************************************************/
int read_comment_gif( FILE *canal, int verbose) {
  /***************************************************************************
   * Called by check_gif
   *           read_gif
   *           read_gif_rgb
   *           interpreta_extension_gif 
  ***************************************************************************/
  unsigned char *buff;
  unsigned char quest;
  unsigned int size;
  int error=OK;

  if(verbose) Warning("Coment");

  for(quest=getc(canal);quest!=0x00;quest=getc(canal))    {
    size=(int)quest;
    TrackNullAlloc( buff=(unsigned char*)calloc(size,sizeof(unsigned char)) );
    fread(buff,sizeof(unsigned char),size,canal);
    if(feof(canal)) error=ERROR;
    if(verbose) fprintf(stderr,"%s",buff);
    free(buff);
  }
  
  if(verbose) Warning("");

  return error;
} // end of read_comment_gif


/***************************************************************************/
int read_application_gif(FILE *canal, int verbose ) {
  /***************************************************************************
   * Called by interpreta_extension_gif
  ***************************************************************************/
  unsigned char *buff;
  unsigned char quest;
  int size,ic;
  int error=OK;
  
  TrackNullAlloc( buff=(unsigned char*)calloc(12,sizeof(unsigned char)) );
  fread(buff,sizeof(unsigned char),1,canal);
  if(feof(canal)) return ERROR;
  fread(buff,sizeof(unsigned char),8,canal);
  if(feof(canal)) return ERROR;
  buff[8]='\0';
  fread(buff+9,sizeof(unsigned char),3,canal);
  if(feof(canal)) return ERROR;

  if(verbose)    {
    WarningV("Application: %s\n",buff);
    fprintf(stderr,"Bytes de autentificacion: %d %d %d\n",
	    (int)buff[9],(int)buff[10],(int)buff[11]);
  }
  free(buff);
  
  for(quest=getc(canal),ic=0;quest!=0x00;quest=getc(canal))    {
    size=(int)quest;
    ic+=size;
    TrackNullAlloc( buff=(unsigned char*)calloc(size,sizeof(unsigned char)) );
    fread(buff,sizeof(unsigned char),size,canal);
    if(feof(canal)) error=ERROR;
    free(buff);
  }
  
  if(verbose) WarningV("Tamanno de application: %d\n",ic);
  
  return error;
} // end of read_application_gif


/***************************************************************************/
int interpreta_stream_gif( FILE *canal, unsigned char **encoded) {
  /***************************************************************************
   * Called by read_gif
   *           read_gif_rgb
   ***************************************************************************/
  char debuff[MAXCHARBUFFLENGTH];
  unsigned char quest;
  long int pos;
  int l_code,ib;

  pos = ftell(canal);
  l_code = 0;
  for(quest=getc(canal);quest!=0x00;quest=getc(canal))    {
    l_code+=(int) quest;
    TrackEOF( fseek(canal,(int)quest,SEEK_CUR), 
	      "Error with fseek positioning");
  }
  encoded[0]=(unsigned char *)calloc(l_code,sizeof(unsigned char));
  
  TrackEOF( fseek(canal,pos,SEEK_SET), "Error with fseek positioning");
  for(quest=getc(canal),ib=0;quest!=0x00;quest=getc(canal))    {
    fread(encoded[0]+ib,sizeof(unsigned char),(int)quest,canal);
    ib+=(int)quest;
  }

  return(l_code);
} // end of interpreta_stream_gif


/***************************************************************************/
int write_stream_gif( int l_code, unsigned char *encoded, FILE *canal) {
  /***************************************************************************
   * Called by old_write_gif_rgb
   *           old_write_gif_rgb_animado
   *           write_gif_rgb
   *           write_gif_rgb_animado
   ***************************************************************************/
  unsigned char bsize;
  int ib;
  
  bsize=0xFF;
  for(ib=0;ib<l_code-(int)bsize;ib+=(int)bsize)    {
    fwrite(&bsize,sizeof(char),1,canal);
    fwrite(&(encoded[ib]),sizeof(char),(int)bsize,canal);
  }

  bsize=(char)(l_code-ib);
  fwrite(&bsize,sizeof(char),1,canal);
  fwrite(&(encoded[ib]),sizeof(unsigned char),(int)bsize,canal);
  
  bsize=0x00;
  fwrite(&bsize,sizeof(char),1,canal);

  return OK;
} // end of write_stream_gif


/***************************************************************************/
int annade_decoded_gif(int w0, int letra, int *tab0, int *tab1, int id0, 
		       int *decoded ) {
  /***************************************************************************
   * Called by annade_decoded_gif (recursive call)
   *           read_gif
   *           read_gif_rgb
   ***************************************************************************/
  int letra0,letra1;
  int id=id0;

  if(letra<w0)
      decoded[id++]=letra;

  else    {
    letra0=tab0[letra-w0-2];
    letra1=tab1[letra-w0-2];
    id=annade_decoded_gif(w0,letra0,tab0,tab1,id,decoded);
    id=annade_decoded_gif(w0,letra1,tab0,tab1,id,decoded);
  }

  return id;
} // end of annade_decoded_gif


/***************************************************************************/
int primer_car_gif( int w0, int letra, int *tab0 ) {
 /***************************************************************************/

  if(letra<w0) 
    return letra;
  
  else 
    return primer_car_gif(w0,tab0[letra-w0-2],tab0);

} // end of primer_car_gif


/***************************************************************************/
int old_write_gif_rgb( int dimx, int dimy, char *nombre, char **Red, char **Green,
		       char **Blue ) {
  /***************************************************************************/
  FILE *canal;
  unsigned char buffer[MAXCHARBUFFLENGTH];
  unsigned char desize;
  int *table_c,*freq_c;
  int color_size;
  int aux;

  unsigned char *encoded;
  unsigned char b0;
  int *decoded,*tab0,*tab1;
  int x0,y0;
  int ibe;
  int ic,ie,id;
  int lt,last_cd;
  int ix,iy;
  int bits,bit0,w0,w1,letra;
  int verbose=0,local_ct;

  /*	Creacion de la tabla de color asociada		*/
  color_size=create_table_color_gif(dimx,dimy,Red,Green,Blue,
				      &table_c,&freq_c);

  if(color_size>256) bit0=8;
  else if(color_size==1) bit0=1;
  else bit0=adimensiona_pos(color_size-1);
  desize=(char) (bit0-1);

  TrackNull( canal=fopen(nombre,"wb"), "Error opening file in 'w' mode" );
  sprintf(buffer,"GIF87a");
  fwrite(buffer,sizeof(char),6,canal);

  buffer[1]=(char)(dimx/256);
  buffer[0]=(char)(dimx-256*(int)buffer[1]);
  buffer[3]=(char)(dimy/256);
  buffer[2]=(char)(dimy-256*(int)buffer[3]);
  buffer[4]=0x70|desize; 
  // lo cual significa: no hay tabla de color global, profundidad de color
  // maxima, tabla no ordenada; se introduce, no obstante, el tamanno de
  // la tabla de color para que las imagenes sean reconocibles al xv.

  buffer[5]=0x00; // color de fondo global, no dado
  buffer[6]=0x00; //pixel aspect ratio

  fwrite(buffer,sizeof(char),7,canal);

  //	Comienza la grabacion de la imagen

  buffer[0]=0x2C; //identificador de imagen
  buffer[1]=0x00;
  buffer[2]=0x00; //Left position
  buffer[3]=0x00;
  buffer[4]=0x00; //Top position
  buffer[6]=(char)(dimx/256);
  buffer[5]=(char)(dimx-256*(int)buffer[6]);
  buffer[8]=(char)(dimy/256);
  buffer[7]=(char)(dimy-256*(int)buffer[8]);
  buffer[9]=0xA0|desize; // que significa: tabla local si,
  // entrelazado no, tabla ordenada y desize codifica el tamanno de la tabla

  fwrite(buffer,sizeof(char),10,canal);

  for(ic=0;ic<Min(color_size,256);ic++)    {
    ix=Mod(table_c[ic],dimx);
    iy=table_c[ic]/dimx;
    putc(Red[iy][ix],canal);
    putc(Green[iy][ix],canal);
    putc(Blue[iy][ix],canal);
  }
  if((color_size>1)&&(color_size<256)) aux=dimensiona(color_size);
  else aux=2;
  if(color_size<aux)    
    for(ic=color_size;ic<aux;ic++)	{
      putc(0,canal);
      putc(0,canal);
      putc(0,canal);
    }
    
  /*		Codificacion de acuerdo a la tabla		*/
  decoded=(int *) calloc(dimx*dimy,sizeof(int));
  codify_imagen_con_table_pos_gif(dimx,dimy,Red,Green,Blue,
				    color_size,table_c,freq_c,decoded);

  /*		Codificacion de la imagen		*/
  encoded=(unsigned char *) calloc((3*dimx*dimy)/2,sizeof(char));
  TrackNullAlloc( tab0=(int*)calloc(4096,sizeof(int)) );
  TrackNullAlloc( tab1=(int*)calloc(4096,sizeof(int)) );

  /*		Registro del tamanno inicial de palabra	*/
  if(bit0<2) bit0=2;
  buffer[0]=(char) bit0;
  fwrite(buffer,sizeof(unsigned char),1,canal);
  
  /*		Inicializacion			*/
  w0=(int) pow(2.,(double)bit0);
  encoded[0]=0;
  last_cd=w0;
  w1=2*w0;
  bits=bit0+1;
  
  /*		Graba un reset como primer codigo	*/
  ie=0;
  b0=0x01;
  letra=w0;
  for(ibe=0;ibe<bits;ibe++,letra=letra/2)    {
    if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
    b0=2*b0;
    if(b0==0x00)	{
      b0=0x01;
      ie++;
      encoded[ie]=0;
    }
  }

  /*		Bucle codificador		*/
  for(id=0;id<dimx*dimy;)    {    
    /*		Examinando la cadena		*/
    if(bits==13)	{
      bits=12;
      letra=w0;
    } else {
      if(last_cd==w0)	{
	lt=0;
	bits=bit0+1;
	w1=2*w0;
      }     
      letra=decoded[id++];
      for(ic=0;(ic<lt)&&(id<dimx*dimy);ic++)	{
	if((tab0[ic]==letra)&&(tab1[ic]==decoded[id]))	  {
	  letra=w0+2+ic;
	  id++;
	}
      }
      
      //	Las cadenas aprendidas corresponden siempre a la ultima cadena 
      //		registrada mas el caracter siguiente.
      if(id<dimx*dimy)	{
	tab0[lt]=letra;
	tab1[lt++]=decoded[id];
      }
    }

    /*		Actualizacion de variables		*/
    last_cd=letra;
    
    /*		Registrando el caracter vigente		*/
    
    for(ibe=0;ibe<bits;ibe++,letra=letra/2)      {
      if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
      b0=2*b0;
      if(b0==0x00)	    {
	b0=0x01;
	ie++;
	encoded[ie]=0;
      }
    }
    
    if(lt+w0+1==w1)	{
      w1=2*w1;
      bits++;
    }
    
  }

  if(bits==13) bits=12;
  letra=w0+1;
  for(ibe=0;ibe<bits;ibe++,letra=letra/2)    {
    if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
    b0=2*b0;
    if(b0==0x00)      {
      b0=0x01;
      ie++;
      encoded[ie]=0;
    }
  }
  ie++; //ie representa ahora la longitud del codificado
  
  /*		Fin de la codificacion		*/
  free(decoded);
  free(tab0);
  free(tab1);

  /*	    Grabacion del codigo resultante	*/
  write_stream_gif(ie,encoded,canal);
  free(encoded);

  /*		Trailer y cierre		*/
  buffer[0]=0x3B;
  fwrite(buffer,sizeof(char),1,canal);
  fclose(canal);

  return OK;
} // end of old_write_gif_rgb 


/***************************************************************************/
int old_write_gif_rgb_animado( int dimx, int dimy, int modo, int colort, 
			       char *nombre, char **Red, char **Green, char **Blue) {
  /***************************************************************************/
  FILE *canal;
  char cabecera[]="FIREFOX";
  unsigned char buffer[MAXCHARBUFFLENGTH];
  unsigned char desize;
  int *table_c,*freq_c;
  int aux,color_size;

  unsigned char *encoded;
  unsigned char b0;
  int *decoded,*tab0,*tab1;
  int x0,y0;
  int ibe;
  int ic,ie,id;
  int lt,last_cd;
  int ix,iy;
  int bits,bit0,w0,w1,letra;
  int verbose=0;

  if(modo==0) {
    
    //	MODO 0: Grabacion de la cabecera
    
    TrackNull( canal=fopen(nombre,"wb"), "Error opening file in 'w' mode" );
    sprintf(buffer,"GIF87a");
    fwrite(buffer,sizeof(char),6,canal);
    
    buffer[1]=(char)(dimx/256);
    buffer[0]=(char)(dimx-256*(int)buffer[1]);
    buffer[3]=(char)(dimy/256);
    buffer[2]=(char)(dimy-256*(int)buffer[3]);
    
    buffer[4]=0xFF; 
    // lo cual significa: tabla de color global, profundidad de color
    // maxima, tabla ordenada y de tamanno maximo
    
    buffer[5]=0x00; // color de fondo global: el negro
    buffer[6]=0x00; //pixel aspect ratio
    fwrite(buffer,sizeof(char),7,canal);
    
    for(ic=0;ic<256;ic++)      {
      putc((char)ic,canal);
      putc((char)ic,canal);
      putc((char)ic,canal);
    }
    
    //	Insercion de application asociada a la animacion
    
    buffer[0]=0x21; //extension
    buffer[1]=0xFF; //de application
    buffer[2]=0x0B; //que mide 11 bytes
    buffer[3]='\0';
    strcat(buffer,cabecera); //identificada como NETSCAPE
    buffer[14]=0x03;
    buffer[15]=0x01;
    buffer[16]=0x00;
    buffer[17]=0x00;
    buffer[18]=0x00;
    fwrite(buffer,sizeof(unsigned char),19,canal);
     
    fclose(canal);

  }  else if(modo==1)    {
    //	MOD0 1: Insercion de un frame
      /*	Creacion de la tabla de color asociada		*/
    if(colort)	{
      color_size=create_table_color_gif(dimx,dimy,Red,Green,Blue,
				      &table_c,&freq_c);
      
      if(color_size>256) bit0=8;
      else if(color_size==1) bit0=1;
      else bit0=adimensiona_pos(color_size-1);
      desize=(char) (bit0-1);
    }      else bit0=8;
    
    TrackNull( canal=fopen(nombre,"ab"), "Error opening file in 'a' mode" );
    TrackEOF( fseek(canal,0,SEEK_END), "Error with fseek positioning");
    
    //	Grabacion de la extension grafica de control
    
    buffer[0]=0x21; //extension
    buffer[1]=0xF9; //grafica de control
    buffer[2]=0x04; //que mide 4 bytes
    buffer[3]=0x00; // ni disposicion ni usuario ni transparencia
    buffer[4]=Mod(DELAY,256);
    buffer[5]=DELAY/256;  // tiempo de retardo entre frames
    buffer[6]=0x00; // color transparente
    buffer[7]=0x00; //terminador
    fwrite(buffer,sizeof(unsigned char),8,canal);
    
    //	Comienza la grabacion de la imagen
    buffer[0]=0x2C; //identificador de imagen
    buffer[1]=0x00;
    buffer[2]=0x00; //Left position
    buffer[3]=0x00;
    buffer[4]=0x00; //Top position
    buffer[6]=(char)(dimx/256);
    buffer[5]=(char)(dimx-256*(int)buffer[6]);
    buffer[8]=(char)(dimy/256);
    buffer[7]=(char)(dimy-256*(int)buffer[8]);
    
    if(colort) buffer[9]=0xAF; // que significa:
    // tabla local si,  entrelazado no, tabla ordenada de tamanno maximo, forzado
    else buffer[9]=0x0;
    
    fwrite(buffer,sizeof(char),10,canal);
    
    if(colort)	{
      for(ic=0;ic<Min(color_size,256);ic++)	{
	ix=Mod(table_c[ic],dimx);
	iy=table_c[ic]/dimx;
	putc(Red[iy][ix],canal);
	putc(Green[iy][ix],canal);
	putc(Blue[iy][ix],canal);
      }
      if((color_size>1)&&(color_size<256)) aux=dimensiona(color_size);
      else aux=2;
      
      if(color_size<256)		  
	for(ic=color_size;ic<256;ic++)		{
	  putc(0,canal);
	  putc(0,canal);
	  putc(0,canal);
	}
      
    }
    color_size=Max(color_size,256);
    bit0=8; // La tabla pasa a tener al menos 256 colores, 
    // forzado por compatibilidad con xv.
    
    /*		Codificacion de acuerdo a la tabla		*/
    decoded=(int *) calloc(dimx*dimy,sizeof(int));
    
    if(colort) 
      codify_imagen_con_table_pos_gif(dimx,dimy,
					Red,Green,Blue,color_size,
					table_c,freq_c,decoded);
    else 
      codify_imagen_gris_gif(dimx,dimy,Red,Green,Blue,
			       decoded);
    
    /*		Codificacion de la imagen		*/
    encoded=(unsigned char *) calloc((3*dimx*dimy)/2,sizeof(char));
    TrackNullAlloc( tab0=(int*)calloc(4096,sizeof(int)) );
    TrackNullAlloc( tab1=(int*)calloc(4096,sizeof(int)) );
    
    /*		Registro del tamanno inicial de palabra	*/
    if(bit0<2) bit0=2;
    buffer[0]=(char) bit0;
    fwrite(buffer,sizeof(unsigned char),1,canal);
    
    /*		Inicializacion			*/
    w0=(int) pow(2.,(double)bit0);
    encoded[0]=0;
    last_cd=w0;
    w1=2*w0;
    bits=bit0+1;
    
    /*		Graba un reset como primer codigo	*/    
    ie=0;
    b0=0x01;
    letra=w0;
    for(ibe=0;ibe<bits;ibe++,letra=letra/2)	{
      if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
      b0=2*b0;
      if(b0==0x00)	    {
	b0=0x01;
	ie++;
	encoded[ie]=0;
      }
    }
    
    /*		Bucle codificador		*/
    for(id=0;id<dimx*dimy;)	{      
      /*		Examinando la cadena		*/      
      if(bits==13)	    {
	bits=12;
	letra=w0;
      }	  else	    {
	if(last_cd==w0)	  {
	  lt=0;
	  bits=bit0+1;
	  w1=2*w0;
	}
	
	letra=decoded[id++];
	for(ic=0;(ic<lt)&&(id<dimx*dimy);ic++)	  {
	  if((tab0[ic]==letra)&&(tab1[ic]==decoded[id]))	    {
	    letra=w0+2+ic;
	    id++;
	  }
	}
	
	//	Las cadenas aprendidas corresponden siempre a la ultima cadena 
	//		registrada mas el caracter siguiente.	
	if(id<dimx*dimy)	  {
	  tab0[lt]=letra;
	  tab1[lt++]=decoded[id];
	}
	
      }
      
      /*		Actualizacion de variables		*/      
      last_cd=letra;
      
      /*		Registrando el caracter vigente		*/
      for(ibe=0;ibe<bits;ibe++,letra=letra/2)	    {
	if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
	b0=2*b0;
	if(b0==0x00)		{
	  b0=0x01;
	  ie++;
	  encoded[ie]=0;
	}
      }
      
      if(lt+w0+1==w1)	{
	w1=2*w1;
	bits++;
      }  
    }
     
    if(bits==13) bits=12;
    letra=w0+1;
    for(ibe=0;ibe<bits;ibe++,letra=letra/2)	{
      if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
      b0=2*b0;
      if(b0==0x00)	    {
	b0=0x01;
	ie++;
	encoded[ie]=0;
      }
    }
    ie++; //ie representa ahora la longitud del codificado
    
    /*		Fin de la codificacion		*/    
    free(decoded);
    free(tab0);
    free(tab1);
    
    /*	    Grabacion del codigo resultante	*/    
    write_stream_gif(ie,encoded,canal);
    free(encoded);
    fclose(canal);

  }  else    {    
    //	Modo por defecto: inserta el trailer y libera la memoria
    //	asociada a la tabla de color
    /*		Trailer y cierre		*/
    
    TrackNull( canal=fopen(nombre,"ab"), "Error opening file in 'a' mode" );
    TrackEOF( fseek(canal,0,SEEK_END), "Error with fseek positioning");
    buffer[0]=0x3B;
    fwrite(buffer,sizeof(char),1,canal);
    fclose(canal);    
  }
  
  return OK;
} // end of old_write_gif_rgb_animado


/***************************************************************************/
int write_gif_rgb( int dimx, int dimy, char *nombre, char **Red, char **Green,
		   char **Blue ) {
/***************************************************************************/

  FILE *canal;
  unsigned char buffer[MAXCHARBUFFLENGTH];
  unsigned char desize;
  char auxR,auxG,auxB;

  int *table_ICY;
  int color_size;
  int aux;

  unsigned char *encoded;
  unsigned char b0;
  int *decoded,*tab0,*tab1;
  int dimI,dimC,dimY;
  int x0,y0;
  int ibe;
  int ic,ie,id;
  int lt,last_cd;
  int ix,iy;
  int iI,iC,iY;
  int bits,bit0,w0,w1,letra;
  int verbose=0,local_ct;

  /*	Creacion de la tabla de color asociada		*/

  color_size=create_table_color_gif_rgb_icy(dimx,dimy,Red,Green,Blue,
				  &dimI,&dimC,&dimY,&table_ICY);

  if(color_size>256) bit0=8;
  else if(color_size==1) bit0=1;
  else bit0=adimensiona_pos(color_size-1);
  desize=(char) (bit0-1);

  TrackNull( canal=fopen(nombre,"wb"), "Error opening file in 'w' mode" );
  sprintf(buffer,"GIF87a");
  fwrite(buffer,sizeof(char),6,canal);

  buffer[1]=(char)(dimx/256);
  buffer[0]=(char)(dimx-256*(int)buffer[1]);
  buffer[3]=(char)(dimy/256);
  buffer[2]=(char)(dimy-256*(int)buffer[3]);
  buffer[4]=0x70|desize; 
  // lo cual significa: no hay tabla de color global, profundidad de color
  // maxima, tabla no ordenada; se introduce, no obstante, el tamanno de
  // la tabla de color para que las imagenes sean reconocibles al xv.

  buffer[5]=0x00; // color de fondo global, no dado
  buffer[6]=0x00; //pixel aspect ratio
	
  fwrite(buffer,sizeof(char),7,canal);

  //	Comienza la grabacion de la imagen
  buffer[0]=0x2C; //identificador de imagen
  buffer[1]=0x00;
  buffer[2]=0x00; //Left position
  buffer[3]=0x00;
  buffer[4]=0x00; //Top position
  buffer[6]=(char)(dimx/256);
  buffer[5]=(char)(dimx-256*(int)buffer[6]);
  buffer[8]=(char)(dimy/256);
  buffer[7]=(char)(dimy-256*(int)buffer[8]);
  buffer[9]=0xA0|desize; // que significa: tabla local si,
  // entrelazado no, tabla ordenada y desize codifica el tamanno de la tabla

  fwrite(buffer,sizeof(char),10,canal);

  for(ic=0;ic<color_size;ic++)    {
    iI=table_ICY[ic]%dimI;
    iC=((table_ICY[ic]-iI)/dimI)%dimC;
    iY=(table_ICY[ic]-iI-dimI*iC)/(dimI*dimC);
    RGB_ICY(dimI,dimC,dimY,1,&auxR,&auxG,&auxB,&iI,&iC,&iY);
    
    putc(auxR,canal);
    putc(auxG,canal);
    putc(auxB,canal);
  }
  if((color_size>1)&&(color_size<256)) aux=dimensiona(color_size);
  else aux=2;
  if(color_size<aux)    {
    for(ic=color_size;ic<aux;ic++)	{
      putc(0,canal);
      putc(0,canal);
      putc(0,canal);
    }
  }
  
  /*		Codificacion de acuerdo a la tabla		*/
  decoded=(int *) calloc(dimx*dimy,sizeof(int));
  codify_imagen_con_table_ICY(dimx,dimy,dimI,dimC,dimY,Red,Green,Blue,
				color_size,table_ICY,decoded);
  free(table_ICY);

  /*		Codificacion de la imagen		*/
  encoded=(unsigned char *) calloc((3*dimx*dimy)/2,sizeof(char));
  TrackNullAlloc( tab0=(int*)calloc(4096,sizeof(int)) );
  TrackNullAlloc( tab1=(int*)calloc(4096,sizeof(int)) );
  
  /*		Registro del tamanno inicial de palabra	*/
  if(bit0<2) bit0=2;
  buffer[0]=(char) bit0;
  fwrite(buffer,sizeof(unsigned char),1,canal);
  
  /*		Inicializacion			*/
  w0=(int) pow(2.,(double)bit0);
  encoded[0]=0;
  last_cd=w0;
  w1=2*w0;
  bits=bit0+1;
  
  /*		Graba un reset como primer codigo	*/
  ie=0;
  b0=0x01;
  letra=w0;
  for(ibe=0;ibe<bits;ibe++,letra=letra/2)    {
    if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
    b0=2*b0;
    if(b0==0x00)      {
      b0=0x01;
      ie++;
      encoded[ie]=0;
    }
  }
  
  /*		Bucle codificador		*/
  for(id=0;id<dimx*dimy;)    {
    /*		Examinando la cadena		*/
    if(bits==13)      {
      bits=12;
      letra=w0;
    }      else      {
      if(last_cd==w0)	    {
	lt=0;
	bits=bit0+1;
	w1=2*w0;
      }
      
      letra=decoded[id++];
      for(ic=0;(ic<lt)&&(id<dimx*dimy);ic++)	{
	if((tab0[ic]==letra)&&(tab1[ic]==decoded[id]))	  {
	  letra=w0+2+ic;
	  id++;
	}
      }
      
      //	Las cadenas aprendidas corresponden siempre a la ultima cadena 
      //		registrada mas el caracter siguiente.      
      if(id<dimx*dimy)	    {
	tab0[lt]=letra;
	tab1[lt++]=decoded[id];
      }
      
    }
    
    /*		Actualizacion de variables		*/
    last_cd=letra;
    

    /*		Registrando el caracter vigente		*/
    for(ibe=0;ibe<bits;ibe++,letra=letra/2)	{
      if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
      b0=2*b0;
      if(b0==0x00)	{
	b0=0x01;
	ie++;
	encoded[ie]=0;
      }
    }
    
    if(lt+w0+1==w1)	{
      w1=2*w1;
      bits++;
    }
    
  }
  if(bits==13) bits=12;
  letra=w0+1;
  for(ibe=0;ibe<bits;ibe++,letra=letra/2)    {
    if(letra%2) encoded[ie]=encoded[ie]|b0;
    b0=2*b0;
    if(b0==0x00)	{
      b0=0x01;
      ie++;
      encoded[ie]=0;
    }
    }
  ie++; //ie representa ahora la longitud del codificado
    
  /*		Fin de la codificacion		*/
  
  free(decoded);
  free(tab0);
  free(tab1);


  /*	    Grabacion del codigo resultante	*/
  write_stream_gif(ie,encoded,canal);
  free(encoded);

  /*		Trailer y cierre		*/
  buffer[0]=0x3B;
  fwrite(buffer,sizeof(char),1,canal);
  fclose(canal);

  return OK;
} // end of write_gif_rgb


/***************************************************************************/
int write_gif_rgb_animado( int dimx, int dimy, int modo, int colort, 
			   char *nombre, char **Red, char **Green, char **Blue ) {
  /***************************************************************************/
  
  FILE *canal;
  char cabecera[7]="FIREFOX";
  unsigned char buffer[MAXCHARBUFFLENGTH];
  unsigned char desize;
  char auxR,auxG,auxB;
  
  int *table_ICY;
  int aux,color_size;

  unsigned char *encoded;
  unsigned char b0;
  unsigned char bsize;
  int *decoded,*tab0,*tab1;
  int dimI,dimC,dimY;
  int x0,y0;
  int ib, ibe;
  int ic,ie,id;
  int lt,last_cd;
  int ix,iy;
  int iI,iC,iY;
  int bits,bit0,w0,w1,letra;
  int verbose=0;

  switch(modo)    {
 
  case 0:    
    //	MODO 0: Grabacion de la cabecera
    TrackNull( canal=fopen(nombre,"wb"), "Error opening file in 'w' mode" );
    sprintf(buffer,"GIF87a");
    fwrite(buffer,sizeof(char),6,canal);
    
    buffer[1]=(char)(dimx/256);
    buffer[0]=(char)(dimx-256*(int)buffer[1]);
    buffer[3]=(char)(dimy/256);
    buffer[2]=(char)(dimy-256*(int)buffer[3]);
    
    buffer[4]=0xFF; 
    // lo cual significa: tabla de color global, profundidad de color
    // maxima, tabla ordenada y de tamanno maximo    
    
    buffer[5]=0x00; // color de fondo global: el negro
    buffer[6]=0x00; //pixel aspect ratio
    fwrite(buffer,sizeof(char),7,canal);    
    
    for(ic=0;ic<256;ic++)	{
      putc((char)ic,canal);
      putc((char)ic,canal);
      putc((char)ic,canal);
    }    
    
    //	Insercion de application asociada a la animacion   
    buffer[0]=0x21; //extension
    buffer[1]=0xFF; //de application
    buffer[2]=0x0B; //que mide 11 bytes
    buffer[3]='\0';
    strcat(buffer,cabecera); //identificada como FIREFOX
    buffer[10]=0x03;
    buffer[11]=0x01;
    buffer[12]=0x00;
    buffer[13]=0x00;
    buffer[14]=0x00;
    fwrite(buffer,sizeof(unsigned char),15,canal);
    fclose(canal);
    break;

  case 1:
    //	MOD0 1: Insercion de un frame   
    /*	Creacion de la tabla de color asociada		*/
    if(colort)	{
      color_size=create_table_color_gif_rgb_icy(dimx,dimy,Red,Green,Blue,
				      &dimI,&dimC,&dimY,
				      &table_ICY);
      
      if(color_size>256) bit0=8;
      else if(color_size==1) bit0=1;
      else bit0=adimensiona_pos(color_size-1);
      desize=(char) (bit0-1);
    }  else bit0=8;
    
    TrackNull( canal=fopen(nombre,"ab"), "Error opening file in 'a' mode" );
    TrackEOF( fseek(canal,0,SEEK_END), "Error with fseek positioning");
    
    //	Grabacion de la extension grafica de control
    buffer[0]=0x21; //extension
    buffer[1]=0xF9; //grafica de control
    buffer[2]=0x04; //que mide 4 bytes
    buffer[3]=0x00; // ni disposicion ni usuario ni transparencia
    buffer[4]=Mod(DELAY,256);
    buffer[5]=DELAY/256;  // tiempo de retardo entre frames
    buffer[6]=0x00; // color transparente
    buffer[7]=0x00; //terminador
    fwrite(buffer,sizeof(unsigned char),8,canal);
    
    //	Comienza la grabacion de la imagen
    buffer[0]=0x2C; //identificador de imagen
    buffer[1]=0x00;
    buffer[2]=0x00; //Left position
    buffer[3]=0x00;
    buffer[4]=0x00; //Top position
    buffer[6]=(char)(dimx/256);
    buffer[5]=(char)(dimx-256*(int)buffer[6]);
    buffer[8]=(char)(dimy/256);
    buffer[7]=(char)(dimy-256*(int)buffer[8]);
    
    if(colort) buffer[9]=0xAF; // que significa:
    // tabla local si,  entrelazado no, tabla ordenada de tamanno maximo, forzado
    else buffer[9]=0x0;
    
    fwrite(buffer,sizeof(char),10,canal);
    
    if(colort)      {
      for(ic=0;ic<color_size;ic++)	    {
	iI=table_ICY[ic]%dimI;
	iC=((table_ICY[ic]-iI)/dimI)%dimC;
	iY=(table_ICY[ic]-iI-dimI*iC)/(dimI*dimC);
	RGB_ICY(dimI,dimC,dimY,1,&auxR,&auxG,&auxB,&iI,&iC,&iY);
	
	putc(auxR,canal);
	putc(auxG,canal);
	putc(auxB,canal);
      }
      if((color_size>1)&&(color_size<256)) 
	aux=dimensiona(color_size);
      else aux=2;
      
      if(color_size<256)	{
	for(ic=color_size;ic<256;ic++)	  {
	  putc(0,canal);
	  putc(0,canal);
	  putc(0,canal);
	}
      }
    }
    color_size=Max(color_size,256);
    bit0=8; // La tabla pasa a tener al menos 256 colores, 
    // forzado por compatibilidad con xv.
    
    /*		Codificacion de acuerdo a la tabla		*/
    decoded=(int *) calloc(dimx*dimy,sizeof(int));
    
    codify_imagen_con_table_ICY(dimx,dimy,dimI,dimC,dimY,
				  Red,Green,Blue,
				  color_size,table_ICY,
				  decoded);
    free(table_ICY);

    /*		Codificacion de la imagen		*/
    encoded=(unsigned char *) calloc((3*dimx*dimy)/2,sizeof(char));
    TrackNullAlloc( tab0=(int*)calloc(4096,sizeof(int)) ); 
    TrackNullAlloc( tab1=(int*)calloc(4096,sizeof(int)) );
    
    /*		Registro del tamanno inicial de palabra	*/
    if(bit0<2) bit0=2;
    buffer[0]=(char) bit0;
    fwrite(buffer,sizeof(unsigned char),1,canal);
    
    /*		Inicializacion			*/    
    w0=(int) pow(2.,(double)bit0);
    encoded[0]=0;
    last_cd=w0;
    w1=2*w0;
    bits=bit0+1;
    
    /*		Graba un reset como primer codigo	*/
    ie=0;
    b0=0x01;
    letra=w0;
    for(ibe=0;ibe<bits;ibe++,letra=letra/2)	{
      if(letra%2) encoded[ie]=encoded[ie]|b0;
      b0=2*b0;
      if(b0==0x00)	{
	b0=0x01;
	ie++;
	encoded[ie]=0;
      }
    }
    
    /*		Bucle codificador		*/
    for(id=0;id<dimx*dimy;)      {
      
      /*		Examinando la cadena		*/      
      if(bits==13)	    {
	bits=12;
	letra=w0;
      }      else	{
	if(last_cd==w0)	  {
	  lt=0;
	  bits=bit0+1;
	  w1=2*w0;
	}
	
	letra=decoded[id++];
	for(ic=0;(ic<lt)&&(id<dimx*dimy);ic++)	  {
	  if((tab0[ic]==letra)&&(tab1[ic]==decoded[id]))	    {
	    letra=w0+2+ic;
	    id++;
	  }
	}
	
	//	Las cadenas aprendidas corresponden siempre a la ultima cadena 
	//		registrada mas el caracter siguiente.
	if(id<dimx*dimy)		{
	  tab0[lt]=letra;
	  tab1[lt++]=decoded[id];
	}
	
      }
      
      /*		Actualizacion de variables		*/           
      last_cd=letra;
      
      /*		Registrando el caracter vigente		*/
      for(ibe=0;ibe<bits;ibe++,letra=letra/2)	{
	if(letra%2) encoded[ie]=encoded[ie]|b0;
	b0=2*b0;
	if(b0==0x00)	  {
	  b0=0x01;
	  ie++;
	  encoded[ie]=0;
	}
      }
      
      if(lt+w0+1==w1)       {
	w1=2*w1;
	bits++;
      }
      
    }
    if(bits==13) bits=12;
    letra=w0+1;
    for(ibe=0;ibe<bits;ibe++,letra=letra/2)	{
      if(letra%2) encoded[ie]=encoded[ie]|b0;
      b0=2*b0;
      if(b0==0x00)	{
	b0=0x01;
	ie++;
	encoded[ie]=0;
      }
    }
    ie++; //ie representa ahora la longitud del codificado
    
    /*		Fin de la codificacion		*/    
    free(decoded);
    free(tab0);
    free(tab1);
    
    /*	    Grabacion del codigo resultante	*/
    write_stream_gif(ie,encoded,canal);
    free(encoded);
    fclose(canal);
    break;

  case 2:
  default:
    //	Modo por defecto: inserta el trailer y libera la memoria
    //	asociada a la tabla de color
    /*		Trailer y cierre		*/

    TrackNull( canal=fopen(nombre,"ab"), "Error opening file in 'a' mode" );
    TrackEOF( fseek(canal,0,SEEK_END), "Error with fseek positioning");
    buffer[0]=0x3B;
    fwrite(buffer,sizeof(char),1,canal);
    fclose(canal);
    break;
  }
  
  return OK;
} // end of write_gif_rgb_animado


/***************************************************************************/
int create_table_color_gif_rgb_icy( int dimx, int dimy, char **Red, char **Green, 
			  char **Blue, int *dimI, int *dimC, int *dimY,
			  int **table_ICY ) {
  /***************************************************************************
   * Called by write_gif_rgb
   *           write_gif_rgb_animado
   ***************************************************************************/
  char auxR,auxG,auxB;
  double *freq;
  int *histo;
  int color_size;
  int ix,iy,iI,iC,iY,it=0;
  
  /*   Initializing Intensity-Chrominance-Yellow space   */
  *dimI=256; // 2**8
  *dimC=64;  // 2**6
  *dimY=64;  // 2**6. With these choices, the ICY space is described with
  // 2**20 values, 16 times less than the maximum possible.
  // This represents 4 Mb in memory, reasonable nowadays.

  TrackNullAlloc( histo=(int*)calloc(*dimI**dimC**dimY,sizeof(int)) );

  /*	Obtaining the color table already reduced in ICY space	*/

  color_size=0;
  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<dimx;ix++)      {
      RGB_ICY(*dimI,*dimC,*dimY,0,&Red[iy][ix],&Green[iy][ix],&Blue[iy][ix],
	      &iI,&iC,&iY);
      
      if(histo[iI+*dimI*(iC+*dimC*iY)]==0) color_size++;
      histo[iI+*dimI*(iC+*dimC*iY)]+=1;
    }

  for(;color_size>256;) { // performing chromatical reduction
    color_size=0;
    *dimI/=2;
    *dimC/=2;
    *dimY/=2;
    for(iY=0;iY<*dimY;iY++)
      for(iC=0;iC<*dimC;iC++)
	for(iI=0;iI<*dimI;iI++)	  {
	  histo[iI+*dimI*(iC+*dimC*iY)]=
	    histo[2*iI+2**dimI*(2*iC+2**dimC*2*iY)]+
	    histo[2*iI+1+2**dimI*(2*iC+2**dimC*2*iY)]+
	    histo[2*iI+2**dimI*(2*iC+1+2**dimC*2*iY)]+
	    histo[2*iI+1+2**dimI*(2*iC+1+2**dimC*2*iY)]+
	    histo[2*iI+2**dimI*(2*iC+2**dimC*(2*iY+1))]+
	    histo[2*iI+1+2**dimI*(2*iC+2**dimC*(2*iY+1))]+
	    histo[2*iI+2**dimI*(2*iC+1+2**dimC*(2*iY+1))]+
	    histo[2*iI+1+2**dimI*(2*iC+1+2**dimC*(2*iY+1))];
	  
	  if(histo[iI+*dimI*(iC+*dimC*iY)]>0) color_size++;
	}
      }

  /*   Creating the color table and its frequency   */
  *table_ICY=(int *) calloc(*dimI**dimC**dimY,sizeof(int));
  TrackNullAlloc( freq=(double*)calloc(color_size,sizeof(double)) );
  
  for(iY=0;iY<*dimY;iY++)    
    for(iC=0;iC<*dimC;iC++)      
      for(iI=0;iI<*dimI;iI++)	
	if(histo[iI+*dimI*(iC+*dimC*iY)]>0)		{
	  table_ICY[0][it]=iI+*dimI*(iC+*dimC*iY);
	  freq[it++]=(double)histo[iI+*dimI*(iC+*dimC*iY)];
	}
  
  /*   Ordering table according to the frequency value   */
  quicksort_ref(0,color_size-1,table_ICY[0],freq);

  /*        Freeing memory before finishing              */
  free(histo);
  free(freq);

  return(color_size);
} // end of create_table_color_gif_rgb_icy


/***************************************************************************/
int RGB_ICY( int dimI, int dimC, int dimY, int mode, char *Red, char *Green,
	     char *Blue, int *iI, int *iC, int *iY ) {
/***************************************************************************/
  double dR,dG,dB;
  double dI,dC,dY;
  double norm;
  
  if(mode==0)    {
    dR=((double)*Red)/256.;      if(dR<0.) dR+=1.;
    dG=((double)*Green)/256.;    if(dG<0.) dG+=1.;
    dB=((double)*Blue)/256.;     if(dB<0.) dB+=1.;
    
    dI=WR*dR+WG*dG+WB*dB;
    dC=0.5*(dG-dR+1.);
    dY=0.25*(dG+dR-2.*dB+2);
    
    *iI=(int)(dI*(double)dimI);    if(*iI>dimI) *iI=dimI-1;
    *iC=(int)(dC*(double)dimC);    if(*iC>dimC) *iC=dimC-1;
    *iY=(int)(dY*(double)dimY);    if(*iY>dimI) *iY=dimY-1;
    
  }  else    {
    dI=((double)*iI)/((double)dimI);
    dC=2.*((double)*iC)/((double)dimC)-1.;
    dY=4.*((double)*iY)/((double)dimY)-2.;
    
    norm=2.*(WR+WG+WB);
    dR=(2.*dI+WB*dY-(2.*WG+WB)*dC)/norm;
    dG=(2.*dI+WB*dY+(2.*WR+WB)*dC)/norm;
    dB=(2.*dI-(WR+WG)*dY+(WR-WG)*dC)/norm;
    
    if(dR<0) dR=0.;
    if(dG<0) dG=0.;
    if(dB<0) dB=0.;
    
    dR=256*dR;
    dG=256*dG;
    dB=256*dB;
    *Red=(((int)dR)>255)?(char)255:(char)(int)dR;
    *Green=(((int)dG)>255)?(char)255:(char)(int)dG;
    *Blue=(((int)dB)>255)?(char)255:(char)(int)dB;  
  }
  
  return OK;
} // end of RGB_ICY


/***************************************************************************/
void codify_imagen_con_table_ICY( int dimx, int dimy, int dimI, int dimC,
				    int dimY, char **Red, char **Green, 
				    char **Blue, int color_size, 
				    int *table_ICY, int *decoded ) {
  /***************************************************************************
   * Called by write_gif_rgb
   *           write_gif_rgb_animado
   ***************************************************************************/

  int apren;
  int ix,iy,ic;
  int iI,iC,iY,iI0,iC0,iY0;

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      RGB_ICY(dimI,dimC,dimY,0,&Red[iy][ix],&Green[iy][ix],&Blue[iy][ix],
	      &iI0,&iC0,&iY0);
      
      for(ic=0,apren=0;(ic<color_size)&&(apren==0);ic++)	    {
	iI=table_ICY[ic]%dimI;
	iC=((table_ICY[ic]-iI)/dimI)%dimC;
	iY=(table_ICY[ic]-iI-dimI*iC)/(dimI*dimC);		
	if((iI==iI0)&&(iC==iC0)&&(iY==iY0))	  {
	  apren=1;
	  decoded[ix+dimx*(dimy-1-iy)]=ic;
	}	
      }
      
      if(apren==0)
	fprintf(stderr,"Error: Color no encontrado en la tabla: %d %d\n",ix,iy);
    }
  
} // end of codify_imagen_con_table_ICY




/***************************************************************************/
int create_table_color_gif( int dimx, int dimy, char **Red, char **Green, 
			  char **Blue, int **table_c, int **freq_c ) {
  /***************************************************************************
   * Called by old_write_gif_rgb
   *           old_write_gif_rgb_animado
   *           reduccion_cromatica_gif
   ***************************************************************************/
  int color_size=0;
  int ic,ic0,ic1,ix,iy,ix0,iy0;
  int t0,f0;
  double d_cr,d_cr0;
  double dR,dG,dB,dR0,dG0,dB0,dl;
  int apren;

  TrackNullAlloc( table_c[0]=(int*)calloc(dimx*dimy,sizeof(int)) );
  TrackNullAlloc( freq_c[0]=(int*)calloc(dimx*dimy,sizeof(int)) );
  for(ic=0;ic<dimx*dimy;ic++) freq_c[0][ic]=0;

  /*		Obteniendo la tabla de color		*/
  for(ic=0;ic<dimx*dimy;ic++)    {
    ix = Mod(ic,dimx);
    iy = ic/dimx;
    for(ic0=0,apren=1;(ic0<color_size)&&(apren==1);ic0++)      {
      ix0=Mod(table_c[0][ic0],dimx);
      iy0=table_c[0][ic0]/dimx;
      if((Red[iy0][ix0]==Red[iy][ix])&&
	 (Green[iy0][ix0]==Green[iy][ix])&&
	 (Blue[iy0][ix0]==Blue[iy][ix])) apren=0;
    }
    if(apren==1)	{
      table_c[0][color_size]=ic;
      freq_c[0][color_size++]=1;
    }      else freq_c[0][ic0]++;
  }
  
  /*		Ordenando la tabla de color		*/
  for(ic0=color_size-1,apren=1;(ic0>0)&&(apren==1);ic0--)    {
    apren=0;
    for(ic=0;ic<ic0;ic++)      
      if(freq_c[0][ic+1]>freq_c[0][ic])	{
	apren=1;
	f0=freq_c[0][ic];
	t0=table_c[0][ic];
	freq_c[0][ic]=freq_c[0][ic+1];
	table_c[0][ic]=table_c[0][ic+1];
	freq_c[0][ic+1]=f0;
	table_c[0][ic+1]=t0;
      }
  }
  
  /*	Truncando la tabla de color si es preciso	*/
  if(color_size>256)    
    for(ic=256;ic<color_size;ic++)	{
      ix0=Mod(table_c[0][ic],dimx);
      iy0=table_c[0][ic]/dimx;
      
      dR0=WR*(double)Mod((int)Red[iy0][ix0],256);
      dG0=WG*(double)Mod((int)Green[iy0][ix0],256);
      dB0=WB*(double)Mod((int)Blue[iy0][ix0],256);
      d_cr0=1e30;
      ic0=0;
      for(ic0=0;ic1<256;ic1++)	{
	ix=Mod(table_c[0][ic1],dimx);
	iy=table_c[0][ic1]/dimx;
	
	dR=WR*(double)Mod((int)Red[iy][ix],256);
	dG=WG*(double)Mod((int)Green[iy][ix],256);
	dB=WB*(double)Mod((int)Blue[iy][ix],256);
	dl=dR+dG+dB-dR0-dG0-dB0;
	d_cr=(dR-dR0)*(dR-dR0)+(dG-dG0)*(dG-dG0)+
	  (dB-dB0)*(dB-dB0)+dl*dl;
	if(d_cr<d_cr0)	  {
	  ic0=ic1;
	  d_cr0=d_cr;
	}
      }
      freq_c[0][ic0]+=freq_c[0][ic];
      freq_c[0][ic]=ic0;
    }
  
  return(color_size);
} // end of create_table_color_gif

/***************************************************************************/
int reduccion_cromatica_gif( int dimx, int dimy, int factor0, char **Red,
			      char **Green, char **Blue ) {
 /***************************************************************************
  * Called by write_video_block
  *           write_video_RGB_block
 ***************************************************************************/
 char **R, **G, **B;
  double d_cr0,d_cr;
  double dR,dG,dB,dR0,dG0,dB0,dl;
  int *t_c,*f_c;
  int idR,idG,idB;
  int dbl=256,dblx,dbly,dimbx,dimby;
  int ix,iy,ibx,iby;
  int ic0,ic1,c_s;
  int factor;
  int verbose=0;

  if(dimx>=dbl)    {
    dimbx=dimx/dbl;
    dblx=dbl;
  }  else    {
    dimbx=1;
    dblx=dimx;
  }
  if(dimy>=dbl)    {
    dimby=dimy/dbl;
    dbly=dbl;
  }  else    {
    dimby=1;
    dbly=dimy;
  }
  
  TrackNullAlloc( R=cmatrix2D(dbly,dblx) );
  TrackNullAlloc( G=cmatrix2D(dbly,dblx) );
  TrackNullAlloc( B=cmatrix2D(dbly,dblx) );
  
  for( c_s=0,factor=dimensiona(factor0); (c_s<256)&&(factor>0);
       factor=factor/2)    {
    for(iy=0;iy<dbly;iy++)	
      for(ix=0;ix<dblx;ix++)	    {
	idR = idG = idB = 0;
	for(iby=0;iby<dimby;iby++)	
	  for(ibx=0;ibx<dimbx;ibx++)	    {
	    idR+=Mod((int)Red[iy*dimby+iby][ix*dimbx+ibx],256);
	    idG+=Mod((int)Green[iy*dimby+iby][ix*dimbx+ibx],256);
	    idB+=Mod((int)Blue[iy*dimby+iby][ix*dimbx+ibx],256);
	  }
	idR=idR/(factor*dimbx*dimby);
	idG=idG/(factor*dimbx*dimby);
	idB=idB/(factor*dimbx*dimby);
	R[iy][ix]=(char)(factor*idR+factor/2);
	G[iy][ix]=(char)(factor*idG+factor/2);
	B[iy][ix]=(char)(factor*idB+factor/2);
      }
    c_s=create_table_color_gif(dblx,dbly,R,G,B,&t_c,&f_c);
    if((c_s<256)&&(factor>1))      {
      free(t_c);
      free(f_c);
    }
  }
  
  if(verbose)
    fprintf(stderr,"Colores obtenidos en la reduccion cromatica: %d\n",
	   c_s);
  c_s=Min(c_s,256);

  for(iy=0;iy<dimy;iy++)  
    for(ix=0;ix<dimx;ix++)	{
      dR0=WR*(double)Mod((int)Red[iy][ix],256);
      dG0=WG*(double)Mod((int)Green[iy][ix],256);
      dB0=WB*(double)Mod((int)Blue[iy][ix],256);
      ic0=0;
      d_cr0=1e30;
      for(ic1=0;(ic1<c_s)&&(d_cr0>0.);ic1++)	    {
	ibx=Mod(t_c[ic1],dblx);
	iby=t_c[ic1]/dblx;
	
	dR=WR*(double)Mod((int)R[iby][ibx],256);
	dG=WG*(double)Mod((int)G[iby][ibx],256);
	dB=WB*(double)Mod((int)B[iby][ibx],256);
	dl=dR+dG+dB-dR0-dG0-dB0;
	d_cr=(dR-dR0)*(dR-dR0)+(dG-dG0)*(dG-dG0)+
	  (dB-dB0)*(dB-dB0)+dl*dl;
	
	if(d_cr<d_cr0)		{
	  ic0=ic1;
	  d_cr0=d_cr;
	}
      }
      ibx=Mod(t_c[ic0],dblx);
      iby=t_c[ic0]/dblx;
      Red[iy][ix]=R[iby][ibx];
      Green[iy][ix]=G[iby][ibx];
      Blue[iy][ix]=B[iby][ibx];
    }
  
  free(t_c);
  free(f_c);

  free_cmatrix2D(R,dbly);
  free_cmatrix2D(G,dbly);
  free_cmatrix2D(B,dbly);

  return OK;
} // end of reduccion_cromatica_gif


/***************************************************************************/
void codify_imagen_con_table_pos_gif( int dimx, int dimy,
					char **Red, char **Green, char **Blue,
					int color_size, int *table_c, 
					int *freq_c, int *decoded ) {
  /***************************************************************************
   * Called by old_write_gif_rgb
   *           old_write_gif_rgb_animado
   ***************************************************************************/
  int ix,iy,ic,ix0,iy0,apren;

  for(iy=0;iy<dimy;iy++)  
    for(ix=0;ix<dimx;ix++)	{
      for(ic=0,apren=0;(ic<256)&&(apren==0);ic++)	{
	ix0=Mod(table_c[ic],dimx);
	iy0=table_c[ic]/dimx;
	if((Red[iy][ix]==Red[iy0][ix0])&&
	   (Green[iy][ix]==Green[iy0][ix0])&&
	   (Blue[iy][ix]==Blue[iy0][ix0]))	  {
	  apren=1;
	  decoded[ix+dimx*iy]=ic;
	}
      }
      for(ic=256;(ic<color_size)&&(apren==0);ic++)	{
	ix0=Mod(table_c[ic],dimx);
	iy0=table_c[ic]/dimx;
	if((Red[iy][ix]==Red[iy0][ix0])&&
	   (Green[iy][ix]==Green[iy0][ix0])&&
	     (Blue[iy][ix]==Blue[iy0][ix0]))    {
	  apren=1;
	  decoded[ix+dimx*iy]=freq_c[ic];
	}
      }
      if(apren==0)
	WarningVV("Error: Color no encontrado en la tabla: %d %d\n",ix,iy);
    }

  free(table_c);
  free(freq_c);
} // end of codify_imagen_con_table_pos_gif


/***************************************************************************/
void codify_imagen_con_table_color_gif( int dimx, int dimy,
					  char **Red, char **Green, char **Blue,
					  int color_size, char *Rg, char *Gg, char *Bg, 
					  int *freq_c, int *decoded ) {
  /***************************************************************************/
  int ix,iy,ic,ic0,d_cr,d_cr0;

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)	{
      
      d_cr0=(Red[iy][ix]-Rg[0])*(Red[iy][ix]-Rg[0])+
	(Green[iy][ix]-Gg[0])*(Green[iy][ix]-Gg[0])+
	(Blue[iy][ix]-Bg[0])*(Blue[iy][ix]-Bg[0]); 
      ic0=0;
      for(ic=1;(ic<color_size)&&(d_cr0>0);ic++)	    {
	d_cr=(Red[iy][ix]-Rg[ic])*(Red[iy][ix]-Rg[ic])+
	  (Green[iy][ix]-Gg[ic])*(Green[iy][ix]-Gg[ic])+
	  (Blue[iy][ix]-Bg[ic])*(Blue[iy][ix]-Bg[ic]); 
	if(d_cr<d_cr0)	  {
	  ic0=ic;
	  d_cr0=d_cr;
	}
      }
      if(d_cr0>0) fprintf(stderr,"Color no presente en la tabla\n");
      if(ic<256) decoded[ix+iy*dimx]=ic;
      else decoded[ix+iy*dimx]=freq_c[ic];
    }

} // end of codify_imagen_con_table_color_gif


/***************************************************************************/
void codify_imagen_gris_gif(int dimx, int dimy, char **Red, char **Green, 
			      char **Blue, int *decoded ) {
  /***************************************************************************
   * Called by old_write_gif_rgb_animado
   ***************************************************************************/
  int ix,iy,ic,ic0,d_cr,d_cr0;
  int dR,dG,dB;

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      dR=(int)Red[iy][ix];
      if(dR<0) dR+=256;
      dG=(int)Green[iy][ix];
      if(dG<0) dG+=256;
      dB=(int)Blue[iy][ix];
      if(dB<0) dB+=256;
      d_cr0=3*65536;
      ic0=0;
      for(ic=0;(ic<256)&&(d_cr0>0);ic++,dR--,dG--,dB--)	{
	d_cr=dR*dR+dG*dG+dB*dB;
	if(d_cr<d_cr0)	  {
	  ic0=ic;
	  d_cr0=d_cr;
	}
      }
      decoded[ix+iy*dimx]=ic0;
    }

} // end of codify_imagen_gris_gif


/***************************************************************************/
int table_pos_table_color_gif(int color_size, int *table_c,
			       char **Rg, char **Gg, char **Bg, int dimx, int dimy, 
			       char **Red, char **Green, char **Blue ) {
/***************************************************************************/
  int ix,iy,ic;

  TrackNullAlloc( Rg[0]=(char*)calloc(color_size,sizeof(char)) );
  TrackNullAlloc( Gg[0]=(char*)calloc(color_size,sizeof(char)) );
  TrackNullAlloc( Bg[0]=(char*)calloc(color_size,sizeof(char)) );

  for(ic=0;ic<color_size;ic++)    {
    ix=Mod(table_c[ic],dimx);
    iy=table_c[ic]/dimx;
    Rg[0][ic]=Red[iy][ix];
    Gg[0][ic]=Green[iy][ix];
    Bg[0][ic]=Blue[iy][ix];
  }
      
  // free(table_c);

  return OK;
} // end of table_pos_table_color_gif

/* OLD:
   void write_foto(int dimx, int dimy, char* nombre, double **data) {
   write_foto_block(dimx,dimy,1.,nombre,data);
   }
*/

/***************************************************************************/
int write_foto_block_limites(int dimx, int dimy, 
			      char* nombre, double min_c, double max_c, 
			      double **data ) {
/***************************************************************************/
  char **foto;
  int sizex,sizey;
  int ix,iy,kx,ky;
  double block;

#ifdef _PARSE_IO_PARAMETERS_
  block = p_io->block_out,
#else
    block = BLOCKOUT;
#endif /*!_PARSE_IO_PARAMETERS_*/
  
  sizex=dimx*block, sizey=dimy*block;
  TrackNullAlloc( foto=cmatrix2D(sizey,sizex) );

  prepare_foto_block_limites(dimx,dimy,min_c,max_c,block,data,foto);

  write_gif_rgb(sizex,sizey,nombre,foto,foto,foto);

  free_cmatrix2D(foto,sizey);

  return OK;
} // end of write_foto_block_limites


/***************************************************************************/
int write_foto_vec_block( int dimx, int dimy, char* nombre,
			  double **vx, double **vy ) {
/***************************************************************************/
  double **data;
  int ix,iy;
  double block;
  
#ifdef _PARSE_IO_PARAMETERS_
  block=p_io->block_out;
#else
  block=BLOCKOUT;
#endif /*!_PARSE_IO_PARAMETERS_*/

  TrackNullAlloc( data=matrix2D(dimy,dimx) );

  for(iy=0;iy<dimy;iy++)   
    for(ix=0;ix<dimx;ix++)	{
      data[iy][ix]=sqrt(vx[iy][ix]*vx[iy][ix]+
			vy[iy][ix]*vy[iy][ix]);
      if(data[iy][ix]>1e-30) data[iy][ix]=log(data[iy][ix]);
      else data[iy][ix]=-30.*log(10.);
    }
  write_foto_block(dimx,dimy,block,nombre,data);

  free_matrix2D(data,dimy);
  
  return OK;
} // end of write_foto_vec_block


/***************************************************************************/
int write_foto_4(int dimx, int dimy, char* nombre, double **data ) {
/***************************************************************************/
  char **foto;

  TrackNullAlloc( foto=cmatrix2D(dimy,dimx) );

  prepare_foto_4(dimx,dimy,data,foto);
  write_gif_rgb(dimx,dimy,nombre,foto,foto,foto);

  free_cmatrix2D(foto,dimy);

  return OK;
} // end of write_foto_4


/***************************************************************************/
int write_binary_block( int dimx, int dimy, double block, char* name,
			char **bin ) {
/***************************************************************************
 * Called by write_binary_foto
***************************************************************************/

  char **foto;
  int sizex,sizey;
  int mm[2], min, max;
  int new, nval;
  int flag_change=FALSE;
  int flag_unary, flag_bin;

  sizex=dimx*block, sizey=dimy*block;

  nval=cop_find(dimx,dimy,bin,2,mm);
  
  /*  if((nval=cop_find(dimx,dimy,bin,2,mm)) == FALSE) { 
      WarningV("Image stored in %s is not a binary image",name);  
      }  else {
  */
  
  /* Binary image */
  
  flag_bin = TRUE;
  if(nval == 2) {
    min = Min(mm[0],mm[1]);
    max = Max(mm[0],mm[1]);
    flag_unary = FALSE;
    
  } else { /* Particular case of an unary image */
    flag_unary = TRUE;
    min = max = mm[0];
    if(abs(max) <= 1)
      if((max=255*abs(max)) > 128) new = 255;
      else new = 0;
    else if(abs(max) < 128) new = 0;
    else if(abs(max) >= 128) new = 255;
  }
  
  if(min!=0 || max!=255.)   flag_change = TRUE;
  /* } */
  
  if(block!=1. || flag_bin==TRUE || flag_change==TRUE) {
    
    TrackNullAlloc( foto=cmatrix2D(sizey,sizex) );
    prepare_char_block( dimx,  dimy,  block, bin, foto );
    if(flag_change == TRUE) { /* that can also happen only with a binary
			       * image (i.e. when flag_bin is TRUE) */
      if(flag_unary == FALSE) {
	cop_change( sizex, sizey, foto, min, 0, NULL );
	cop_change( sizex, sizey, foto, max, 255, NULL );
      } else 
	cop_change( sizex, sizey, foto, min, new, NULL );
    }
    
  } else foto = bin;
  
  write_gif_rgb(sizex,sizey,name,foto,foto,foto);
  //  write_pgm(sizex,sizey,FALSE,name,foto);

  if(block!=1. || flag_bin==TRUE || flag_change==TRUE)     
    free_cmatrix2D(foto,sizey);

  return OK;
} // end of write_binary_block


/***************************************************************************/
int write_RGB_block( int dimx, int dimy, char* nombre,
		     char **Red, char **Green, char **Blue) {
  /***************************************************************************/
  char **fR,**fG,**fB;
  int sizex,sizey;
  int ix,iy;
  double block;
  
#ifdef _PARSE_IO_PARAMETERS_
  block=p_io->block_out;
#else
  block=BLOCKOUT;
#endif /*!_PARSE_IO_PARAMETERS_*/
   sizex=dimx*block, sizey=dimy*block;

  TrackNullAlloc( fR=cmatrix2D(sizey,sizex) );
  TrackNullAlloc( fG=cmatrix2D(sizey,sizex) );
  TrackNullAlloc( fB=cmatrix2D(sizey,sizex) );

  prepare_char_block(dimx,dimy,block,Red,fR);
  prepare_char_block(dimx,dimy,block,Green,fG);
  prepare_char_block(dimx,dimy,block,Blue,fB);

  write_gif_rgb(sizex,sizey,nombre,fR,fG,fB);

  free_cmatrix2D(fR,sizey);
  free_cmatrix2D(fG,sizey);
  free_cmatrix2D(fB,sizey);

  return OK;
} // end of write_RGB_block


/***************************************************************************/
int write_foto_color_block( int dimx, int dimy, int n_cr, double block,
			    char *nombre, double ***data) {
  /***************************************************************************
   * Called by write_foto_color
   ***************************************************************************/
  char ***foto;
  int dimxf,dimyf;
  int icr;

  dimxf=dimx*block, dimyf=dimy*block;

  TrackNullAlloc( foto=cmatrix3D(n_cr,dimyf,dimxf) );

  for(icr=0;icr<n_cr;icr++) 
    prepare_foto_block(dimx,dimy,block,
		       data[icr],foto[icr]);

  if(n_cr==1) write_gif_rgb(dimxf,dimyf,nombre,foto[0],foto[0],foto[0]);
  else write_ppm(dimxf,dimyf,1,nombre,foto[0],foto[1],foto[2]);

  free_cmatrix3D(foto,n_cr,dimyf);

  return OK;
} // end of write_foto_color_block


/***************************************************************************/
int write_video_block( int dimx, int dimy, int n_cr, double block, 
			int modo, char *nombre, double ***data) {
  /***************************************************************************
   * Called by write_foto
  ***************************************************************************/
  char ***video;
  int dimxf,dimyf;
  int icr;
  
  dimxf=dimx*block, dimyf=dimy*block;
  TrackNullAlloc( video=cmatrix3D(n_cr,dimyf,dimxf) );

  if(modo==1) {
      
    for(icr=0;icr<n_cr;icr++) prepare_foto_block(dimx,dimy,block,
						 data[icr],video[icr]);
    
    if(n_cr==1) write_gif_rgb_animado(dimxf,dimyf,modo,1,nombre,
				  video[0],video[0],video[0]);
    else {
      reduccion_cromatica_gif(dimxf,dimyf,8,	
			      video[0],video[1],video[2]);
      write_gif_rgb_animado(dimxf,dimyf,modo,1,nombre,
			video[0],video[1],video[2]);
    }
  }  else write_gif_rgb_animado(dimxf,dimyf,modo,1,nombre,
			    video[0],video[0],video[0]);
  
  free_cmatrix3D(video,n_cr,dimyf);

  return OK;
} // end of write_video_block


/***************************************************************************/
int write_video_RGB_block( int dimx, int dimy, double block, int modo, 
			   char *nombre, char **Red, char **Green, char **Blue) {
  /***************************************************************************
   * Called by write_foto_RGB
   ***************************************************************************/
  
  char ***video;
  int dimxf,dimyf;
  
  dimxf=dimx*block, dimyf=dimy*block;  
  TrackNullAlloc( video=cmatrix3D(3,dimyf,dimxf) );

  if(modo==1)   {
    prepare_char_block(dimx,dimy,block,Red,video[0]);
    prepare_char_block(dimx,dimy,block,Green,video[1]);
    prepare_char_block(dimx,dimy,block,Blue,video[2]);
    
    //reduccion_cromatica_gif(dimxf,dimyf,8,	
    //video[0],video[1],video[2]);
    
    write_gif_rgb_animado(dimxf,dimyf,modo,1,nombre,
		      video[0],video[1],video[2]);

  }  else if(modo==0) write_gif_rgb_animado(dimxf,dimyf,modo,0,nombre,
					video[0],video[1],video[2]);
  else write_gif_rgb_animado(dimxf,dimyf,modo,1,nombre,
			 video[0],video[1],video[2]);
  
  free_cmatrix3D(video,3,dimyf);

  return OK;
} // end of write_video_RGB_block


/***************************************************************************/
int prepare_foto_block( int dimx, int dimy, double block, double **data,
			char **foto) {
  /***************************************************************************
   * Called by write_foto_color_block
   *           write_video_block
   *           prepare_foto_block_limites
   *           write_foto_block
   ***************************************************************************/
  char buff;
  int xmax,ymax;
  int ix,iy,ibx,iby,beff;
  double maximo,minimo,prov;

  if(block >= 1.)    {
    beff = (int) block;
    maximo = minimo = data[0][0];
    for(ix=0;ix<dimx;ix++)      
      for(iy=0;iy<dimy;iy++)	{
	maximo = fMax(maximo,data[iy][ix]);
	minimo = fMin(minimo,data[iy][ix]);
      }
    
    if(maximo-minimo > 1e-30)     
      for(ix=0;ix<dimx;ix++)   
	for(iy=0;iy<dimy;iy++)		{
	  buff=(char)((int)(255*(data[iy][ix]-minimo)
			    /(maximo-minimo)));
	  for(iby=0;iby<beff;iby++)	    
	    for(ibx=0;ibx<beff;ibx++)			
	      foto[iy*beff+iby][ix*beff+ibx] = buff;	   	 
	}
    
    else       
      for(iy=0;iy<dimy*beff;iy++)      
	for(ix=0;ix<dimx*beff;ix++)
	  foto[iy][ix]=0x7f;
    
  }  else    {
    beff = (int)(1./block);
    xmax = dimx/beff;
    ymax = dimy/beff;
    maximo = -1e30;
    minimo = 1e30;
    for(iy=0;iy<ymax;iy++)
      for(ix=0;ix<xmax;ix++)	    {
	prov = 0.;
	for(iby=0;iby<beff;iby++)	
	  for(ibx=0;ibx<beff;ibx++)
	    prov += data[iy*beff+iby][ix*beff+ibx];
	maximo = fMax(maximo,prov);
	minimo = fMin(minimo,prov);
      }
    
    if(maximo > minimo) 
      for(iy=0;iy<ymax;iy++)
	for(ix=0;ix<xmax;ix++)	    {
	  prov = 0.;
	  for(iby=0;iby<beff;iby++)
	    for(ibx=0;ibx<beff;ibx++)
	      prov += data[iy*beff+iby][ix*beff+ibx];
	  foto[iy][ix] = (char)((int)(255*(prov-minimo)
				    /(maximo-minimo)));
	}
    else
      for(iy=0;iy<ymax;iy++)
	for(ix=0;ix<xmax;ix++)
	  foto[iy][ix] = 0x7f;
  }
  
  return OK;
} // end of prepare_foto_block


/***************************************************************************/
int prepare_foto_fijo_block( int dimx, int dimy, int levels, 
			     double **data, char **foto ) {
  /***************************************************************************/
  char buff;
  int xmax,ymax;
  int ix,iy,ibx,iby,beff;
  double block;
  double prov;
  
#ifdef _PARSE_IO_PARAMETERS_
  block=p_io->block_out;
#else
  block=BLOCKOUT;
#endif /*!_PARSE_IO_PARAMETERS_*/

  if(block>=1.)    {
    beff=(int) block;
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)	{
	prov=255.*data[iy][ix]/((double)levels-1.);
	if(prov>255) buff=(char) 255;
	else if(prov<0) buff=(char) 0;
	else buff=(char)(int)prov;
	for(iby=0;iby<beff;iby++)
	  for(ibx=0;ibx<beff;ibx++)
	    foto[iy*beff+iby][ix*beff+ibx]=buff;
      }

  }  else    {
    beff=(int)(1./block);
    xmax=dimx/beff;
    ymax=dimy/beff;
    for(iy=0;iy<ymax;iy++)
      for(ix=0;ix<xmax;ix++)	{
	prov=0.;
	for(iby=0;iby<beff;iby++)
	  for(ibx=0;ibx<beff;ibx++)
	    prov+=data[iy*beff+iby][ix*beff+ibx];
	prov=255.*prov/((double)(beff*beff*levels));
	if(prov>255) buff=(char) 255;
	else if(prov<0) buff=(char) 0;
	else buff=(char)(int)prov;
	foto[iy][ix]=buff;
      }
  }
  
  return OK;
} // end of prepare_foto_fijo_block


/***************************************************************************/
int prepare_foto_block_limites( int dimx, int dimy, double min_c, 
				double max_c, double block, double **data, 
				char **foto) {
  /***************************************************************************
  * Called by write_foto_block_limites
  *           write_image_reference
  ***************************************************************************/
  int ix,iy;
  double **buff;

  TrackNullAlloc( buff=matrix2D(dimy,dimx) );

  for(ix=0;ix<dimx;ix++)
    for(iy=0;iy<dimy;iy++)      {
      if(data[iy][ix]>max_c) buff[iy][ix]=max_c;
      else if(data[iy][ix]<min_c) buff[iy][ix]=min_c;
      else buff[iy][ix]=data[iy][ix];
    }
  prepare_foto_block(dimx,dimy,block,buff,foto);

  free_matrix2D(buff,dimy);

  return OK;
} // end of prepare_foto_block_limites


/***************************************************************************/
int prepare_foto_log( int dimx, int dimy, double **data, char **foto ) {
/***************************************************************************/
  int ix,iy;
  double maximo,minimo,buff;

  maximo=-1e30, minimo=1e30;
  
  for(ix=0;ix<dimx;ix++)  
    for(iy=0;iy<dimy;iy++)      {
      buff=fabs(data[iy][ix]);
      if(buff>1e-30)	    {
	maximo=fMax(maximo,log(buff));
	minimo=fMin(minimo,log(buff));
      }
    }
  
  fprintf(stderr,"%f   %f\n",maximo,minimo);
  if(maximo-minimo>1e-30)    {
    maximo-=minimo;
    for(ix=0;ix<dimx;ix++)  
      for(iy=0;iy<dimy;iy++)	{
	buff=fabs(data[iy][ix]);
	if(buff>1e-30) buff=log(buff)-minimo;
	else buff=0.;
	foto[iy][ix]=(char)((int)(255*buff/maximo));
      }

  }  else    {
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	foto[iy][ix]=0x7f;    
  }

  return OK;
} // end of prepare_foto_log


/***************************************************************************/
int prepare_foto_4(int dimx, int dimy, double **data, char **foto) {
  /***************************************************************************
   * Called by write_foto_4
   ***************************************************************************/
  int ix,iy;
  double maximo,minimo,prov;
                        
  maximo = minimo = data[0][0];

  for(ix=0;ix<dimx;ix++)
    for(iy=0;iy<dimy;iy++)	{
      maximo=fMax(maximo,data[iy][ix]);
      minimo=fMin(minimo,data[iy][ix]);   
    }               
     	    
  if(maximo-minimo>1e-30)
    for(iy=0;iy<dimy/2;iy++)
      for(ix=0;ix<dimx/2;ix++)	    {
	prov=255.*(data[iy+dimy/2][ix+dimx/2]-minimo)
	  /(maximo-minimo);
	foto[iy][ix]=(char)((int)(prov));
	
	prov=255.*(data[iy+dimy/2][ix]-minimo)
	  /(maximo-minimo);
	foto[iy][ix+dimx/2]=(char)((int)(prov));
	
	prov=255.*(data[iy][ix+dimx/2]-minimo)
	  /(maximo-minimo);
	foto[iy+dimy/2][ix]=(char)((int)(prov));
	
	prov=255.*(data[iy][ix]-minimo)
	  /(maximo-minimo);
	foto[iy+dimy/2][ix+dimx/2]=(char)((int)(prov));
      }
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	foto[iy][ix]=0x7f;
  
  return OK;
} // end of prepare_foto_4            


/***************************************************************************/
int prepare_foto_log_4(int dimx, int dimy, double **data, char **foto) {
/***************************************************************************/
  int ix,iy;
  double maximo,minimo,prov,buff;

  maximo=-1e30, minimo=1e30;

  for(ix=0;ix<dimx;ix++)
    for(iy=0;iy<dimy;iy++)	{
      buff=fabs(data[iy][ix]);
      if(buff>1e-30)	{
	maximo=fMax(maximo,log(buff));
	minimo=fMin(minimo,log(buff));
      }
    }                              
  
  if(maximo-minimo>1e-30)  {
    maximo-=minimo;
    for(iy=0;iy<dimy/2;iy++)	{
      for(ix=0;ix<dimx/2;ix++)	{
	buff=fabs(data[iy+dimy/2][ix+dimx/2]);
	if(buff>1e-30) prov=(log(buff)-minimo)/maximo;
	else prov=0.;
	foto[iy][ix]=(char)((int)(255*prov));
      }
      for(ix=dimx/2;ix<dimx;ix++)	{
	buff=fabs(data[iy+dimy/2][ix-dimx/2]);
	if(buff>1e-30) prov=(log(buff)-minimo)/maximo;
	else prov=0.;
	foto[iy][ix]=(char)((int)(255*prov));
      }
    }
    
    for(iy=dimy/2;iy<dimy;iy++)      {
      for(ix=0;ix<dimx/2;ix++)	{
	buff=fabs(data[iy-dimy/2][ix+dimx/2]);
	if(buff>1e-30) prov=(log(buff)-minimo)/maximo;
	else prov=0.;
	foto[iy][ix]=(char)((int)(255*prov));
      }
      for(ix=dimx/2;ix<dimx;ix++)	    {
	buff=fabs(data[iy-dimy/2][ix-dimx/2]);
	if(buff>1e-30) prov=(log(buff)-minimo)/maximo;
	else prov=0.;
	foto[iy][ix]=(char)((int)(255*prov));
      }
    }

  }  else    {
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	foto[iy][ix]=0x7f;
  }
  
  return OK; 
} // end of prepare_foto_log_4


/***************************************************************************/
int prepare_char_block( int dimx, int dimy, double block, char **grayin,
			char **grayout ) {
  /***************************************************************************
   * Called by write_RGB_block
   *           write_video_RGB_block
   *           write_binary_block
   ***************************************************************************/
  int sizex,sizey;
  int ix,iy,ibx,iby;
  int beff,buff,prov;

  sizex = dimx*block, sizey = dimy*block;

  if(block >= 1)  {
    beff = (int)block;
    for( iy=0; iy<dimy; iy++ )      
      for( ix=0; ix<dimx; ix++ )
	for( iby=0; iby<beff; iby++ )
	  for( ibx=0; ibx<beff; ibx++ )
	    grayout[iby+beff*iy][ibx+beff*ix] =
	      grayin[iy][ix];
    
  }  else    {
    beff = (int)(1./block);
    for( iy=0; iy<sizey; iy++ )    
      for( ix=0; ix<sizex; ix++ )	    {
	buff = 0;
	for( iby=0; iby<beff; iby++ )
	  for( ibx=0; ibx<beff; ibx++ )	      {
	    prov = (int)grayin[iy*beff+ibx][ix*beff+ibx];
	    if(prov < 0) prov += 256;
	    buff += prov;
	  }
	buff = buff/(beff*beff);
	grayout[iy][ix] = (char)buff;
      }  
  }
  
  return OK;
} // end of prepare_char_block


/***************************************************************************/
int write_gray_foto( int dimx, int dimy, int dimz, int dimv, int iz, 
		     int bin, char *base, char *ext, char ***cont ) {
  /***************************************************************************/

  char name[MAXNAMELENGTH];
  int ic;

  if(dimz == 1) {
    if((dimv==1)||(dimv==3)) {
      sprintf(name,"%s.%s.gif",base,ext);
      if(dimv==1)   write_pgm( dimx, dimy, bin, name, cont[0] );
      else write_ppm( dimx, dimy, bin, name, cont[0], cont[1], cont[2] );
    } else	
      for(ic=0;ic<dimv;ic++)	{
	sprintf(name,"%s.%s-V%02d.gif",base,ext,ic);
	write_pgm( dimx, dimy, bin, name, cont[ic] );
      }
    
  } else { // dimz > 1
    if((dimv==1)||(dimv==3)) {
      sprintf(name,"%s.%s-Z%02d.gif",base,ext,iz);
      if(dimv==1)   write_pgm( dimx, dimy, bin, name, cont[0] );
      else write_ppm( dimx, dimy, bin, name, cont[0], cont[1], cont[2] );
    } else	
      for(ic=0;ic<dimv;ic++)	{
	sprintf(name,"%s.%s-Z%02d-V%02d.gif",base,ext,iz,ic);
	write_pgm( dimx, dimy, bin, name, cont[ic] );
      }
  }
  
  return OK;
} // end of write_gray_foto


/***************************************************************************/
int write_pgm( int dimx, int dimy, int bin, char *nombre_out, 
	       char **cont ) {
  /***************************************************************************/
 FILE* canal;
 char cab_raw[]=
   "P5\n# Created with write_pgm\n";
 char cab_asc[]=
   "P2\n# Created with write_pgm\n";
 int ix,iy;
 int dato;
 int l;

 if(bin)    {
   TrackNull( canal=fopen(nombre_out,"wb"), "Error opening file in 'w' mode" );
   l=strlen(cab_raw);
   fwrite(cab_raw,sizeof(char),l,canal);
   fprintf(canal,"%d %d\n%d\n",dimx,dimy,255);
   for(iy=0;iy<dimy;iy++) fwrite(cont[iy],sizeof(char),dimx,canal);

 }  else    {
   TrackNull( canal=fopen(nombre_out,"wt"), "Error opening file in 'w' mode" );
   l=strlen(cab_asc);
   fwrite(cab_asc,sizeof(char),l,canal);
   fprintf(canal,"%d %d\n%d\n",dimx,dimy,255);
   for(iy=0;iy<dimy;iy++)
     for(ix=0;ix<dimx;ix++)       {
       dato=(int) cont[iy][ix];
       if(dato<0) dato=256+dato;
       fprintf(canal,"%d\n",dato);
     }
 }
 
  fclose(canal);
  
  return OK;
} // end of write_pgm


/***************************************************************************/
int write_ppm( int dimx, int dimy, int bin, char *nombre_out, char **Red, 
	       char **Green, char **Blue ) {
/***************************************************************************/
  char cab_raw[]=
    "P6\n # Created with write_ppm\n";
  char cab_asc[]=
    "P3\n # Created with write_ppm\n";
  int ix,iy,l;
  int datoR,datoG,datoB;
  int salida;
  FILE* canal;
         
  if(bin)    {
    TrackNull( canal=fopen(nombre_out,"wb"), "Error opening file in 'w' mode" );
    l=strlen(cab_raw);
    fwrite(cab_raw,sizeof(char),l,canal);
    fprintf(canal,"%d %d\n%d\n",dimx,dimy,255);
    for(iy=0;iy<dimy;iy++)	
      for(ix=0;ix<dimx;ix++)	{
	fwrite(&Red[iy][ix],sizeof(char),1,canal);
	fwrite(&Green[iy][ix],sizeof(char),1,canal);
	fwrite(&Blue[iy][ix],sizeof(char),1,canal);
      }

  }  else    {
    TrackNull( canal=fopen(nombre_out,"wt"), "Error opening file in 'w' mode" );
    l=strlen(cab_asc);
    fwrite(cab_asc,sizeof(char),l,canal);
    fprintf(canal,"%d %d\n%d\n",dimx,dimy,255);
    
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)	{
	datoR=(int) Red[iy][ix];
	if(datoR<0) datoR=256+datoR;
	datoG=(int) Green[iy][ix];
	if(datoG<0) datoG=256+datoG;
	datoB=(int) Blue[iy][ix];
	if(datoB<0) datoB=256+datoB;
	fprintf(canal,"%d  %d  %d\n",datoR,datoG,datoB);
      }
  }
  
  fclose(canal);
  
  return OK;
} // end of write_ppm


/***************************************************************************/
int write_foto( int dimx, int dimy, int dimv, int dimz, int iz, 
		char *base, char *ext,  double ***data ) {
  /***************************************************************************/
  
  char nombre[MAXNAMELENGTH];
  int ic;

  double blout=
#ifdef _PARSE_IO_PARAMETERS_
  p_io->block_out;
#else
  BLOCKOUT;
#endif /*!_PARSE_IO_PARAMETERS_*/

  if(
#ifdef _PARSE_IO_PARAMETERS_
     p_io->flag_video ||
#else
     FLAG_VIDEO ||
#endif /*!_PARSE_IO_PARAMETERS_*/
     (dimz==1))  {
    if((dimv==1)||(dimv==3))    {
      sprintf(nombre,"%s.%s.gif",base,ext);
      if(iz==0) write_video_block(dimx,dimy,dimv,blout,0,nombre,data);
      write_video_block(dimx,dimy,dimv,blout,1,nombre,data);
      if(iz==dimz-1) write_video_block(dimx,dimy,dimv,blout,2,nombre,data);
    }    else	{
      for(ic=0;ic<dimv;ic++)	{
	sprintf(nombre,"%s.%s-V%02d.gif",base,ext,ic);
	if(iz==0) write_video_block(dimx,dimy,1,blout,0,nombre,&(data[ic]));
	write_video_block(dimx,dimy,1,blout,1,nombre,&(data[ic]));
	if(iz==dimz-1) write_video_block(dimx,dimy,1,blout,2,nombre,&(data[ic]));
      }
    }
  }  else    {
    if((dimv==1)||(dimv==3))      {
      sprintf(nombre,"%s.%s-Z%02d.gif",base,ext,iz);
      write_video_block(dimx,dimy,dimv,blout,0,nombre,data);
      write_video_block(dimx,dimy,dimv,blout,1,nombre,data);
      write_video_block(dimx,dimy,dimv,blout,2,nombre,data);
    }    else      {
      for(ic=0;ic<dimv;ic++)	{
	sprintf(nombre,"%s.%s-Z%02d-V%02d.gif",base,ext,iz,ic);
	write_foto_block(dimx,dimy,blout,nombre,data[ic]);
      }
    }
  } 
  
  return OK; 
} // end of write_foto


/***************************************************************************/
int write_mask_block( int dimx, int dimy, double block, char* nombre,
		       char **mask ) {
  /***************************************************************************
   * Called by visualise_gris_trozos
   ***************************************************************************/
  char **foto;
  int sizex,sizey;
  int ix,iy,ibx,iby;
  int beff,prov,buff;

  sizex=dimx*block, sizey=dimy*block;
  TrackNullAlloc( foto=cmatrix2D(sizey,sizex) );
  
  if(block>=1)    {
    beff=(int)block;
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	for(iby=0;iby<beff;iby++)
	  for(ibx=0;ibx<beff;ibx++)
	    foto[iby+beff*iy][ibx+beff*ix]=mask[iy][ix];
    
  } else {
    beff=(int)(1./block);
    for(iy=0;iy<sizey;iy++)	
      for(ix=0;ix<sizex;ix++)	    {
	buff=0;
	for(iby=0;iby<beff;iby++)	    {
	  for(ibx=0;ibx<beff;ibx++)		    {
	    prov=(int)mask[iy*beff+iby][ix*beff+ibx];
	    if(prov<0) prov+=256;
	    buff+=prov;
	  }
	}
	buff=buff/(beff*beff);
	foto[iy][ix]=(char)buff;
      }      
  }
  
  write_gif_rgb(sizex,sizey,nombre,foto,foto,foto);
  free_cmatrix2D(foto,sizey);
  
  return OK;
} // end of write_mask_block


/***************************************************************************/
int write_foto_block( int dimx, int dimy, double block, char* nombre,
		      double **data ) {
  /***************************************************************************
   * Called by write_foto
   *           write_image_foto, 
   *           write_foto_vec_block
   ***************************************************************************/
  
  char **foto;
  int sizex,sizey;
  int ix,iy,kx,ky;

  sizex=dimx*block, sizey=dimy*block;
  TrackNullAlloc( foto=cmatrix2D(sizey,sizex) );

  prepare_foto_block(dimx,dimy,block,data,foto);
  write_gif_rgb(sizex,sizey,nombre,foto,foto,foto);

  free_cmatrix2D(foto,sizey);

  return OK;
} // end of write_foto_block


/***************************************************************************/
int write_binary_foto( int dimx, int dimy, int dimz, int iz, 
		       char *base, char *ext, char **bin ) {
  /***************************************************************************/

  char name[MAXNAMELENGTH];

  if(dimz == 1) 
    sprintf(name,"%s.%s.gif",base,ext);
  else
    sprintf(name,"%s.%s-Z%02d.gif",base,ext,iz);
  
  write_binary_block( dimx, dimy, 
#ifdef _PARSE_IO_PARAMETERS_
		      p_io->block_out,
#else
		      BLOCKOUT,
#endif /*!_PARSE_IO_PARAMETERS_*/
		      name, bin );

  return OK;
} // end of write_binary_foto


/***************************************************************************/
int write_foto_RGB( int dimx, int dimy, int dimv, int dimz, int iz, char *base,
		    char *ext,  char ***Red, char ***Green, char ***Blue) {
  /***************************************************************************
   * Called by  write_image_reference
   ***************************************************************************/
  char nombre[MAXNAMELENGTH];
  int ic;

  double blout=
#ifdef _PARSE_IO_PARAMETERS_
    p_io->block_out;
#else
  BLOCKOUT;
#endif /*!_PARSE_IO_PARAMETERS_*/
  
  if(
#ifdef _PARSE_IO_PARAMETERS_
     p_io->flag_video ||
#else
     FLAG_VIDEO ||
#endif
     (dimz==1))  {
    if(dimv==1)    {
      sprintf(nombre,"%s.%s.gif",base,ext);
      if(iz==0) write_video_RGB_block(dimx,dimy,blout,0,nombre,
				      Red[0],Green[0],Blue[0]);
      write_video_RGB_block(dimx,dimy,blout,1,nombre,
			    Red[0],Green[0],Blue[0]);
      if(iz==dimz-1) write_video_RGB_block(dimx,dimy,blout,2,nombre,
					   Red[0],Green[0],Blue[0]);
    }    else	{
      for(ic=0;ic<dimv;ic++)	    {
	sprintf(nombre,"%s.%s-V%02d.gif",base,ext,ic);
	if(iz==0) write_video_RGB_block(dimx,dimy,blout,0,nombre,
					Red[ic],Green[ic],Blue[ic]);
	write_video_RGB_block(dimx,dimy,blout,1,nombre,
			      Red[ic],Green[ic],Blue[ic]);
	if(iz==dimz-1) write_video_RGB_block(dimx,dimy,blout,2,nombre,
					     Red[ic],Green[ic],Blue[ic]);
      }
    } 
  }  else    {
    if(dimv==1)     {
      sprintf(nombre,"%s.%s-Z%02d.gif",base,ext,iz);
      write_video_RGB_block(dimx,dimy,blout,0,nombre,
			    Red[0],Green[0],Blue[0]);
      write_video_RGB_block(dimx,dimy,blout,1,nombre,
			    Red[0],Green[0],Blue[0]);
      write_video_RGB_block(dimx,dimy,blout,2,nombre,
			    Red[0],Green[0],Blue[0]);
    }     else	{
      for(ic=0;ic<dimv;ic++)	{
	sprintf(nombre,"%s.%s-Z%02d-V%02d.gif",base,ext,iz,ic);
	write_video_RGB_block(dimx,dimy,blout,0,nombre,
			      Red[ic],Green[ic],Blue[ic]);
	write_video_RGB_block(dimx,dimy,blout,1,nombre,
			      Red[ic],Green[ic],Blue[ic]);
	write_video_RGB_block(dimx,dimy,blout,2,nombre,
			      Red[ic],Green[ic],Blue[ic]);
      }
    } 
  }
  
  return OK;
} // end of write_foto_RGB


/***************************************************************************/
void paleta( double rep, char *Red, char *Green, char *Blue) {
  /***************************************************************************
   * Called by visualise_color
   ***************************************************************************/
  double rep2;
  int gamma=0;
  int scR,scG,scB;
  int ic;
  
  rep2=5.*rep;
  ic=(int) rep2;
  if(ic==5) ic=4;
  rep2-=(double)ic;
  rep2=1.-rep2;
  
  switch(ic)    {
  case 0:
    scR=255-gamma;
    scG=255-gamma;
    scB=(255-gamma)*rep2;
    break;
  case 1:
    scR=(255-gamma)*rep2;
    scG=(255-gamma)/2;
    scB=0.;
    break;
  case 2:
    scR=(255-gamma)*(1-rep2);
    scG=0.;//(255-gamma)*rep2;
    scB=0.;
    break;
  case 3:
    scR=(255-gamma)*rep2;
    scG=0.;
    scB=(255-gamma)*(1.-rep2);
    break;
  case 4:
    scR=0.;
    scG=0.;
    scB=(255-gamma)*rep2;
    break;
  }

  *Red=(char)(gamma+scR);
  *Green=(char)(gamma+scG);
  *Blue=(char)(gamma+scB);

} // end of paleta


/***************************************************************************/
int write_image_reference( int dimx, int dimy, int dimv, int dimz, int dimt,
			   int iz, char *base, char **nombre_tipo, char *ext,
			   double ***data ) {
  /***************************************************************************/
  
  char nombre[MAXNAMELENGTH];
  char ***ajustado;
  double mmd0[2],mmd[2];
  int it;
  
  TrackNullAlloc( ajustado=cmatrix3D(dimv,dimy,dimx) );
  
  for(it=0;it<dimt;it++)    {
    extrema(dimx,dimy,data[it],&mmd[0],NULL);
    if(it)      {
      mmd[0]=fMin(mmd[0],mmd0[0]);
      mmd[1]=fMax(mmd[1],mmd0[1]);
    }    else      {
      mmd0[0]=mmd[0];
      mmd0[1]=mmd[1];
    }
  }
  
  for(it=0;it<dimt;it++)	  {
    if(it==3) op_log(dimx,dimy,data[it],data[it], NULL);
    prepare_foto_block_limites(dimx,dimy,mmd0[0],mmd0[1],/*frst,res=*/1.,
			       data[it],ajustado[0]);
    sprintf(nombre,"%s%s",base,nombre_tipo[it]);

    write_foto_RGB(dimx,dimy,dimv,dimz,iz,nombre,ext,ajustado,ajustado,ajustado); 
  }
  
  free_cmatrix3D(ajustado,dimv,dimy);

  return OK;
} // end of write_image_reference


/***************************************************************************/
int write_foto_color( int dimx, int dimy, int dimz, int iz, 
		      char *base, char *ext, double ***data) {
  /***************************************************************************/

  char name[MAXNAMELENGTH];

  if(dimz == 1)        sprintf(name,"%s.%s.ppm",base,ext);
  else                 sprintf(name,"%s.%s-Z%02d.ppm",base,ext,iz);
  
  write_foto_color_block( dimx, dimy, 3, 
#ifdef _PARSE_IO_PARAMETERS_
			  p_io->block_out,
#else
			  BLOCKOUT,
#endif /*!_PARSE_IO_PARAMETERS_*/
			  name, data);
  
  return OK;
} // end of write_foto_color


/***************************************************************************/
int write_image_foto( int leff, int dimy, char *name_in, double **datos) {
  /***************************************************************************/
  char name[90];
  int error;

  sprintf(name,"%s.gif",name_in);
#ifdef _PARSE_IO_PARAMETERS_
  error = write_foto_block(leff,dimy,p_io->block_out,name,datos);
#else
  error = write_foto_block(leff,dimy,BLOCKOUT,name,datos);
#endif /*!_PARSE_IO_PARAMETERS_*/
  
  return error;
} // end of write_image_foto



/***************************************************************************/
int visualise_color( int dimx, int dimy, char *name_in, 
		     double m, double M, double **data ) {
  /***************************************************************************/
  char name[MAXNAMELENGTH];
  char **Red,**Green,**Blue;
  double **dataeff;
  double rep;
  int dimxf,dimyf;
  int ix,iy;
  double blout;
  
#ifdef _PARSE_IO_PARAMETERS_
  blout=p_io->block_out;
#else
  blout=BLOCKOUT;
#endif /*!_PARSE_IO_PARAMETERS_*/

  dimxf=blout*dimx, dimyf=blout*dimy;
  
  TrackNullAlloc( Red=cmatrix2D(dimyf,dimxf) );
  TrackNullAlloc( Green=cmatrix2D(dimyf,dimxf) );
  TrackNullAlloc( Blue=cmatrix2D(dimyf,dimxf) );
  TrackNullAlloc( dataeff=matrix2D(dimyf,dimxf) );
  
  coarseres(dimx,dimy,1./blout,data,dataeff);
  
  for(iy=0;iy<dimyf;iy++)   
    for(ix=0;ix<dimxf;ix++)	{
      rep=(dataeff[iy][ix]-m)/(M-m);
      if(rep < 0.) rep=0.;
      if(rep > 1.) rep=1.;
      paleta(rep,&(Red[iy][ix]),&(Green[iy][ix]),&(Blue[iy][ix]));
    }
  
  // sprintf(nombre,"sing_total.%s.ppm",dext);
  sprintf(name,"%s.ppm",name_in);
  write_ppm(dimxf,dimyf,1,name,Red,Green,Blue);
  
  free_cmatrix2D(Red,dimyf);
  free_cmatrix2D(Green,dimyf);
  free_cmatrix2D(Blue,dimyf);
  free_matrix2D(dataeff,dimyf);

  return OK;
}


/***************************************************************************/
int visualise_gris( int dimx, int dimy, char *name_in, 
		     double m, double M, double **data ) {
/***************************************************************************/
  char name[MAXNAMELENGTH];
  double **aux;
  int ix,iy;
  
  TrackNullAlloc( aux=matrix2D(dimy,dimx) );
  
  for(iy=0;iy<dimy;iy++)	  {
    for(ix=0;ix<dimx;ix++)      {
      aux[iy][ix]=-data[iy][ix];
      if(aux[iy][ix]>-m) aux[iy][ix]=-m;
      if(aux[iy][ix]<-M) aux[iy][ix]=-M;
    }
  }
  // sprintf(nombre,"sing_total.%s.%s.gif",dext,nombre_wv);
  sprintf(name,"%s.gif",name_in);
  write_foto_4(dimx,dimy,name,aux);
  
  free_matrix2D(aux,dimy);
}


