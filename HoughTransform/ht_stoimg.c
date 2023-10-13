/*****************************************************************************
 * Houghtool - The Software Package for efficiency measuring and visualization
 * of the HT and it's variants for line detection
 *
 * Lappeenranta University of Technology, Department of Information Technology
 * Laboratory of Information Processing
 *
 * Authors: Petri Hirvonen, Jouni Ikonen (Jouni.Ikonen@lut.fi)
 *	    Pekka Kultanen, Heikki K{lvi{inen (Heikki.Kalviainen@lut.fi)
 *
 * File:    ht_stoimg.c
 * Purpose: store images
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_formats.h"
#include "ht_cvl.h"

/*****************************************************************************/
FILE *cvl_write_open(name,dim_x,dim_y)
/*

  Open CVL image file for writing.

  Parameters :  name - output image filename
                dim_x,dim_y - image size

*/
char *name;
int dim_x,dim_y;
{
  FILE *data=stdout;
  struct cvl_p cvl_header;

  if (strlen(name))
    if ((data=fopen(name,"wb"))==(FILE *)NULL) {
      err_prnt(name);
      return 0;
    }

  cvl_header.cv_tag = 044520;
  cvl_header.cv_lcom = 0;
  cvl_header.cv_type = 0;
  cvl_header.cv_dims = 2;
  cvl_header.cv_hpls = 1;
  cvl_header.cv_plns = 1;
  cvl_header.cv_rows = dim_y;
  cvl_header.cv_cols = dim_x;
  cvl_header.cv_bnds = 1;
  cvl_header.cv_bp = 8;
  cvl_header.cv_ebb = 0;
  cvl_header.cv_sbb = 0;
  cvl_header.cv_bb = 0;
  cvl_header.cv_sbp = 0;
  cvl_header.cv_ebp = 0;

  if (fwrite(&cvl_header,sizeof(struct cvl_p),1,data)!=1) {
    err_prnt(name);
    fclose(data);
    return 0;
  }

  return data;
}

/*****************************************************************************/
FILE *pgm_write_open(name,dim_x,dim_y,proc)
/*

  Open PGM image file for writing.

  Parameters :  name - output image filename
                dim_x,dim_y - image size

*/
char *name, *proc;
int dim_x,dim_y;
{
  FILE *data=stdout;

  if (strlen(name))
    if ((data=fopen(name,"wb"))==(FILE *)NULL) {
      err_prnt(name);
      return 0;
    }
  if (proc!=(char *)NULL && strlen(proc)) {
    if (fprintf(data,"P5\n# created by %s\n%d %d\n255\n",
                proc,dim_x,dim_y)==EOF) {
      err_prnt(name);
      fclose(data);
      return 0;
    }
  } else
    if (fprintf(data,"P5\n%d %d\n255\n",dim_x,dim_y)==EOF) {
      err_prnt(name);
      fclose(data);
      return 0;
    }

  return data;
}

/*****************************************************************************/
FILE *raw_write_open(name,dim_x,dim_y)
/*

  Open image file for writing.

  Parameters :  name - output image filename
                dim_x,dim_y - image size

*/
char *name;
int dim_x,dim_y;
{
  FILE *data=stdout;

  if (strlen(name))
    if ((data=fopen(name,"wb"))==(FILE *)NULL) {
      err_prnt(name);
      return 0;
    }

  return data;
}

/*****************************************************************************/
pic_to_disk(file_name,pic,dim_x,dim_y,format,proc)
/*

  Write image to disk.

  Parameters :  file_name - output image filename
                pic - image
                dim_x,dim_y - image size
                format - image format
                proc - program name (for PGM comment field)

*/
char *file_name, *proc;
pic_type pic[][MAX_SIZE],dim_x,dim_y,format;
{
  FILE *data;
  int i,j;
  char name[256];

  if (file_name==(char *)NULL)
    name[0]='\0';
  else
    strcpy(name,file_name);

  switch(format) {
    case CVL_PIC: if (strlen(name) && !(check_file_name(name,".cvl")))
                    strcat(name,".cvl");
                  if (!(data=cvl_write_open(name,dim_x,dim_y)))
                    return 0;
                  break;
    case PGM_PIC: if (strlen(name) && !(check_file_name(name,".pgm")))
                    strcat(name,".pgm");
                  if (!(data=pgm_write_open(name,dim_x,dim_y,proc)))
                    return 0;
                  break;
    case SKE_PIC: if (strlen(name) && !(check_file_name(name,".ske")))
                    strcat(name,".ske");
                  if (!(data=raw_write_open(name,dim_x,dim_y)))
                    return 0;
                  break;
    case BIN_PIC: if (strlen(name) && !(check_file_name(name,".bin")))
                    strcat(name,".bin");
                  if (!(data=raw_write_open(name,dim_x,dim_y)))
                    return 0;
                  break;
    default:      return 0;
  }

  if (format==SKE_PIC) { /* four bytes per pixel (try to avoid this) */
    for (i=0; i<dim_y; i++)
      if (fwrite(pic[i],sizeof(pic_type),dim_x,data)!=dim_x) {
        err_prnt(name);
        fclose(data);
        return 0;
      }
  } else
  for (i=0; i<dim_y; i++) /* one byte per pixel */
    for (j=0; j<dim_x; j++)
      if ((fprintf(data,"%c",pic[i][j])) == EOF) {
        err_prnt(name);
        fclose(data);
        return 0;
      }

  fclose(data);
}

/*****************************************************************************/
invert_pic_to_disk(file_name,pic,dim_x,dim_y,format,proc)
/*

  Write image to disk.

  Parameters :  file_name - output image filename
                pic - image
                dim_x,dim_y - image size
                format - image format
                proc - program name (for PGM comment field)

*/
char *file_name, *proc;
pic_type pic[][MAX_SIZE],dim_x,dim_y,format;
{
  FILE *data;
  int i,j;
  char name[256];
  pic_type tmp_pic[MAX_SIZE];

  if (file_name==(char *)NULL)
    name[0]='\0';
  else
    strcpy(name,file_name);

  switch(format) {
    case CVL_PIC: if (strlen(name) && !(check_file_name(name,".cvl")))
                    strcat(name,".cvl");
                  if (!(data=cvl_write_open(name,dim_x,dim_y)))
                    return 0;
                  break;
    case PGM_PIC: if (strlen(name) && !(check_file_name(name,".pgm")))
                    strcat(name,".pgm");
                  if (!(data=pgm_write_open(name,dim_x,dim_y,proc)))
                    return 0;
                  break;
    case SKE_PIC: if (strlen(name) && !(check_file_name(name,".ske")))
                    strcat(name,".ske");
                  if (!(data=raw_write_open(name,dim_x,dim_y)))
                    return 0;
                  break;
    case BIN_PIC: if (strlen(name) && !(check_file_name(name,".bin")))
                    strcat(name,".bin");
                  if (!(data=raw_write_open(name,dim_x,dim_y)))
                    return 0;
                  break;
    default:      return 0;
  }

  if (format==SKE_PIC) { /* four bytes per pixel (try to avoid this) */
    for (i=0; i<dim_y; i++){
	for (j=0; j<dim_x; j++) tmp_pic[j]=pic[i][j]?0:255;
      if (fwrite(tmp_pic,sizeof(pic_type),dim_x,data)!=dim_x) {
        err_prnt(name);
        fclose(data);
        return 0;
      }
    }
  } else
  for (i=0; i<dim_y; i++) /* one byte per pixel */
    for (j=0; j<dim_x; j++)
      if ((fprintf(data,"%c",pic[i][j]?(char)0:(char)255)) == EOF) {
        err_prnt(name);
        fclose(data);
        return 0;
      }

  fclose(data);
}

