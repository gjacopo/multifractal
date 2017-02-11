/* ===================================
** io_netcdf.c
** started on Tue Jan 30 09:47:18 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <io_netcdf.h>
//    NetCDF I/O routines

DataStructure ping_netcdf( char *filename)
/*
   PURPOSE:   Ping for NetCDF data files
*/
{
   DataStructure output=_DataStructure;
   size_t len;

   int *map_d,*map_v;
   int *dimids;

   int D_space;
   int status=0;
   int varnum;
   int ndims,nvars,natts,unldimid;
   int ndims0;

   int ncid,varid,dimid;
   int nmap_v;

   status = nc_open(filename,NC_NOWRITE,&ncid);
   if (status != NC_NOERR) return (output);
   nc_inq(ncid,&ndims,&nvars,&natts,&unldimid);

/* Establishing the map between dimensions and coordinate variables */

   map_d=(int *)calloc(ndims,sizeof(int));
   map_v=(int *)calloc(nvars,sizeof(int));
   nmap_v=nc_map_dv(ncid,ndims,nvars,map_d,map_v);
   if (nmap_v==0)
   {
    free(map_d);
    free(map_v);
    return(output);
   }

/*
   Counting the number of valid variables according to their dimensions
*/

   varnum=-1;
   for (varid=0;(varid<nvars)&&(varnum==-1);varid++)
   {
    dimid=map_v[varid];

    if (dimid==-1) // so this is NOT a coordinate variable
    {
     nc_inq_varndims(ncid,varid,&ndims0);
     if ((ndims0>=1)&&(ndims0<=4)) // let it be at least 1D
                                   // but no more than 4D
     {
      dimids=(int *) calloc(ndims0,sizeof(int));
      nc_inq_vardimid(ncid,varid,dimids);
      switch (ndims0)
      {
       case 1:
        D_space=1;
        varnum++;
        break;
       case 2:
        if (dimids[0]!=unldimid) D_space=2;
        else                     D_space=1;
        varnum++;
        break;
       case 3:
        D_space=2;
        varnum++; // we accept multiple components signals
        break;
       case 4:
        D_space=2;
        if (dimids[0]==unldimid) varnum++; // 2D+comps+time
        break;
      }
      if (varnum==0) // we found the requested variable, the first one
      {
       output.Reader=&simple_netcdf_reader;
       nc_inq_dimlen(ncid,dimids[ndims0-1],&len);
       output.dimx=(int)len;
       output.dimy=1;
       output.dimv=1;
       output.dimz=1;

       switch (ndims0)
       {
        case 1:
         break;
        case 2:
         nc_inq_dimlen(ncid,dimids[ndims0-2],&len);
         if (D_space==1)
         {
          output.dimz=(int)len;
          output.dimy=output.dimv=1;
         }
         else
         {
          output.dimy=(int)len;
          output.dimv=output.dimz=1;
         }
         break;
        case 3:
        case 4:
         nc_inq_dimlen(ncid,dimids[ndims0-2],&len);
         output.dimy=len;
         nc_inq_dimlen(ncid,dimids[0],&len);
         if (dimids[0]==unldimid) output.dimz=(int)len;
         nc_inq_dimlen(ncid,dimids[0],&len);
         if (dimids[0]==unldimid) output.dimz=(int)len;
         else                     output.dimv=(int)len;
         if (ndims0==4)
         {
          nc_inq_dimlen(ncid,dimids[1],&len);
          output.dimv=(int)len;
         }
         break;
       }
      }

/*    Memory release before quiting the variable identification    */

      free(dimids);
     }
    }
   }

/* Closing file and releasing memory */

   nc_close(ncid);
   free(map_v);
   free(map_d);

   return(output);
}
int simple_netcdf_reader( char *filename, DataWindow select, Matrix* data)
/*
   PURPOSE:  The one indicated by the routine name. WARNING! This routine
   blindly trusts the window selection, so be sure it fits your data
*/
{
   size_t *start,*count;
   char missattname[]="missing_value";
   char *buffer;
   nc_type xtype;
   double miss;
   double Nev,data_av;
   int *map_d,*map_v;
   int *dimids;
   int hasunldim;
   int status=0;
   int varnum,attnum,attlen;
   int ndims,nvars,natts,unldimid;
   int ndims0;
   int ncid,varid,dimid;
   int nmap_v;
   int ix;

   nc_open(filename,NC_NOWRITE,&ncid);
   nc_inq(ncid,&ndims,&nvars,&natts,&unldimid);

/* Establishing the map between dimensions and coordinate variables */

   map_d=(int *)calloc(ndims,sizeof(int));
   map_v=(int *)calloc(nvars,sizeof(int));
   nmap_v=nc_map_dv(ncid,ndims,nvars,map_d,map_v);
   if (nmap_v==0)
   {
    free(map_d);
    free(map_v);
    return(-1);
   }

/*
   Counting the number of valid variables according to their dimensions
*/

   varnum=-1;
   for (varid=0;(varid<nvars)&&(varnum==-1);varid++)
   {
    dimid=map_v[varid];
    if (dimid==-1) // so this is NOT a coordinate variable
    {
     nc_inq_varndims(ncid,varid,&ndims0);
     if ((ndims0>=1)&&(ndims0<=4)) // let it be at least 1D but no more than 4D
     {
      dimids=(int *) calloc(ndims0,sizeof(int));
      nc_inq_vardimid(ncid,varid,dimids);
      hasunldim=(dimids[0]==unldimid)?1:0;
      switch (ndims0)
      {
       case 1:
       case 2:
       case 3:
        varnum++; // we accept multiple components signals
        break;
       case 4:
        if (hasunldim) varnum++; // 2D+comps+time
        break;
      }
      if (varnum==0) // we found the requested variable, the first one
      {

/*     Allocating memory for receiving the data     */

       create_matrix(select.dimx,select.dimy,0,data);

/*     Defining the region to be written down in the matrix     */

       start=(size_t *) calloc(ndims0,sizeof(size_t));
       count=(size_t *) calloc(ndims0,sizeof(size_t));

       switch (ndims0)
       {
        case 1:
         start[0]=select.ix0;
         count[0]=select.dimx;
         break;
        case 2:
         start[1]=select.ix0;
         count[1]=select.dimx;
         if (hasunldim)
         {
          start[0]=select.iz0;
          count[0]=1;
         }
         else
         {
          start[0]=select.iy0;
          count[0]=select.dimy;
         }
         break;
        case 3:
        case 4:
         start[ndims0-1]=select.ix0;
         count[ndims0-1]=select.dimx;
         start[ndims0-2]=select.iy0;
         count[ndims0-2]=select.dimy;
         if (ndims0==4)
         {
          start[1]=select.iv0;
          count[1]=1;
         }
         if (hasunldim) start[0]=select.iz0;
         else           start[0]=select.iv0;
         count[0]=1;
         break;
       }

/*     Actual reading     */

       nc_inq_vartype(ncid,varid,&xtype); // NC_CHAR must be treated separately
       if (xtype==NC_CHAR)
       {
        buffer=(char *) calloc(select.dimx*select.dimy,sizeof(char));
        nc_get_vara_text(ncid,varid,start,count,buffer);
        for (ix=0;ix<select.dimx*select.dimy;ix++)
        {
         data->p[ix]=(double)(int)buffer[ix];
         if (data->p[ix]<0.) data->p[ix]+=256.;
        }
        free(buffer);
       }
       else nc_get_vara_double(ncid,varid,start,count,data->p);

/*     Defining the mask, if necessary     */

       status=nc_inq_attid(ncid,varid,missattname,&attnum);
       if (status==NC_NOERR)
       {
        nc_inq_attlen(ncid,varid,missattname,&attlen);
        if (attlen==1) // so we have an actual mask
        {

/*       Allocating memory for the mask       */

         data->pm=(char *) calloc(data->dimx*data->dimy,sizeof(char));

/*       Retrieving the value associated to the masked points       */

         nc_get_att_double(ncid,varid,missattname,&miss);

/*       Masking points conveniently and obtaining the average on valid points       */

         Nev=0.;
         data_av=0.;
         for (ix=0;ix<data->dimx*data->dimy;ix++)
         {
          if (data->p[ix]==miss) data->pm[ix]=MASKED;
          else
          {
           data->pm[ix]=VALID;
           Nev+=1.;
           data_av+=data->p[ix];
          }
         }
         if (Nev>0.) data_av/=Nev;

/*       Second pass to substitute the (eventually large) mask value by average       */

         for (ix=0;ix<data->dimx*data->dimy;ix++)
         {
          if (data->pm[ix]==MASKED) data->p[ix]=data_av;
         }

/*       Refitting matrix in order to be able to referentiate 2D the mask       */

         fit_matrix(data);
        }
       }

/*     Memory release before quiting the reading     */

       free(start);
       free(count);
      }

/*    Memory release before closing the loop    */

      free(dimids);
     }
    }
   }

/* Closing file and releasing memory */

   nc_close(ncid);
   free(map_v);
   free(map_d);

   return(0);
}


int nc_map_dv( int ncid, int ndims, int nvars, int *map_d, int *map_v)
{
   char dimname[NC_MAX_NAME];

   int dimid0;
   int nmap_v;
   int ndims0;
   int varid,dimid;
   int status;

   for (varid=0;varid<nvars;varid++) map_v[varid]=-1;
   nmap_v=nvars;

   for (dimid=0;dimid<ndims;dimid++)
   {
      nc_inq_dimname(ncid,dimid,dimname);
      status=nc_inq_varid(ncid,dimname,&varid);
      if (status==NC_NOERR)
      {
         nc_inq_varndims(ncid,varid,&ndims0);
         if (ndims0>1) map_d[dimid]=-1; // if it is not a coordinate variable
         else
         {
            nc_inq_vardimid(ncid,varid,&dimid0);
            if (dimid0!=dimid)
            {
               map_d[dimid]=-1; // if dimensions don't fit
            }
            else
            {
               map_d[dimid]=varid;
               map_v[varid]=dimid;
               nmap_v--;
            }
         }
      }
      else
      {
         map_d[dimid]=-1;
      }
   }

   return(nmap_v);
}


int write_simple_netcdf( char *dataname, Matrix data)
/*
   PURPOSE: To provide a simple NetCDF write interface for matrices
*/
{
   char filename[90];
   const char missattname[]="missing_value";

   size_t start[4],count[4];

   double *buffer;
   const double missv=1e30;

   int ncid;
   int ndims;
   int dimids[4];
   int varid;
   int ix;
   int lost_it;
   int attlen;

/*
   Bear in mind that the type GeoImage is rather general; in its core it is
   always a tritensor, but we will record it as that if the V coordinate
   is greater than 1, and as a matrix in other case.
*/

   sprintf(filename,"%s.nc",dataname);

   nc_create(filename,NC_WRITE,&ncid);

/* Defining the dimensions */

   ndims=1;
   nc_def_dim(ncid,"X",data.dimx,&lost_it);
   if (data.dimy>1)
   {
      nc_def_dim(ncid,"Y",data.dimy,&lost_it);
      ndims++;
      dimids[1]=0; // Fastest coordinate: X
      dimids[0]=1; // Slowest coordinate: Y
   }
   else dimids[0]=0; // We have only X coordinate

/* Defining the variable */

   varid=0;
   nc_def_var(ncid,data.Name,NC_DOUBLE,ndims,&dimids[0],&lost_it);

/* If required, stamping an attibute for the missing value */

   if (data.Mask!=NULL) nc_put_att_double(ncid,varid,missattname,NC_DOUBLE,1,&missv);

/* Stamping the origin */

   attlen=1+strlen(originname);
   nc_put_att_text(ncid,NC_GLOBAL,"origin",attlen,originname);

/* Stamping the factory seal */

   attlen=1+strlen(factoryname);
   nc_put_att_text(ncid,NC_GLOBAL,"factory_seal",attlen,factoryname);

/* End definitions */

   nc_enddef(ncid);

/* Start writing */

/* Defining the geometry */

   varid=0;
   if (ndims==2)
   {
      start[1]=0;
      count[1]=data.dimx;
      start[0]=0;
      count[0]=data.dimy;
   }
   else
   {
      start[0]=0;
      count[0]=data.dimx;
   }

/* Separating cases if data are masked or not */

   if (data.Mask==NULL) // No mask
   {
      nc_put_vara_double(ncid,varid,&start[0],&count[0],data.p);
   }
   else // we have to padd masked points with missing value previous to write
   {
      buffer=(double *) calloc(data.dimx*data.dimy,sizeof(double));
      for (ix=0;ix<data.dimx*data.dimy;ix++)
      {
         if (data.pm[ix]==MASKED) buffer[ix]=missv;
         else                     buffer[ix]=data.p[ix];
      }
      nc_put_vara_double(ncid,varid,&start[0],&count[0],buffer);
      free(buffer);
   }

/* Closing file */

   nc_close(ncid);

   return(0);
}
