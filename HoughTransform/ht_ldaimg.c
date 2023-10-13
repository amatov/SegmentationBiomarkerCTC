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
 * File:    ht_ldaimg.c
 * Purpose: image loading
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_formats.h"
#include "ht_visilog.h"
#include "ht_cvl.h"

int cols, rows;
char s[256];

/*****************************************************************************/
long conv_long(item)

/*

  Converts MicroSoft-type long integer representation to
  SUN-type long integer and returns it.

  Parameters : item - value to be converted

*/
long *item;
{
  short tmp;
  union {
    short shorts[2];
    unsigned char bytes[4];
    long longs;
  }  temp;

  temp.longs = *item;
  tmp = temp.shorts[1];
  temp.shorts[1] = temp.shorts[0];
  temp.shorts[0] = tmp;

  tmp = temp.bytes[1];
  temp.bytes[1] = temp.bytes[0];
  temp.bytes[0] = tmp;

  tmp = temp.bytes[2];
  temp.bytes[2] = temp.bytes[3];
  temp.bytes[3] = tmp;

  *item=temp.longs;

  return temp.longs;
}

/*****************************************************************************/
check_file_name(fname, ident)
/*

  Check whether fname includes string ident.

  Parameters :  fname - file name string
                ident - string

*/
char *fname, *ident;
{
  int i, check_chars=strlen(fname)-strlen(ident)+1;

  for (i=0; i<check_chars; i++, fname++) {
    if (!strncmp(fname, ident, strlen(ident)))
      return 1;
  }
  return 0;
}

/*****************************************************************************/
image_format(fname)
/*

  Try to find out the input image format.

  Parameters :  fname - file name string

*/
char *fname;
{
  FILE *fd=stdin;
  struct cvl_p cvl_header;
  struct desc_im vis_header;
  char pgm_header_line[256], magic[3];

  if (strlen(fname))
    if( (fd = fopen(fname, "r"))==NULL) {
      err_prnt(fname);
      return NULL;
    }

  /* CVL ?? */
  if (fread(&cvl_header, sizeof(struct cvl_p),1,fd) != 1) {
    err_prnt("reading header");
    if (strlen(fname))
      fclose(fd);
    return NULL;
  }

  if (cvl_header.cv_tag == 044520) {
    if (strlen(fname))
      fclose(fd);
    else
      rewind(fd);
    return CVL_PIC;
  }
  rewind(fd);

  /* PGM ?? */
  while (fgets(pgm_header_line,255,fd) && pgm_header_line[0]=='#');
  sscanf(pgm_header_line,"%2s",magic);

  if (!strncmp(magic,"P5",2)) {
    if (strlen(fname))
      fclose(fd);
    else
      rewind(fd);
    return PGM_PIC;
  }
  rewind(fd);

  /* VIS ?? */
  if (fread(&vis_header,sizeof(struct desc_im),1,fd)!=1) {
    err_prnt("reading header");
    if (strlen(fname))
      fclose(fd);
    return NULL;
  }

  conv_long(&vis_header.i_magic);
  if (vis_header.i_magic == I_MAGICD) {
    if (strlen(fname))
      fclose(fd);
    else
      rewind(fd);
    return VIS_PIC;
  }
  if (strlen(fname))
    fclose(fd);
  else
    rewind(fd);

  /* SKE ?? */
    if (strlen(fname))
      if (check_file_name(fname,".ske"))
        return SKE_PIC;

  /* BIN ?? */
    if (strlen(fname))
      if (check_file_name(fname,".bin") || check_file_name(fname,".raw"))
        return BIN_PIC;

 /* PNT  if a cordinate file */
    if (strlen(fname))
      if (check_file_name(fname,".cor"))
        return COR_PIC;

  sprintf(s,"%s: Unrecognized file format\n",strlen(fname)?fname:"'<stdin>'");
  err_prnt(s);

  return NULL;
}

/*****************************************************************************/
FILE *cvl_open(filename,c,r)
/*

  Open CVL input image file.

  Parameters :  filename - file name string
                c,r - image size

*/
char *filename; int *c,*r;
{
  FILE *fd=stdin;
  char *comment;
  struct cvl_p cvl_header;

  if (strlen(filename))
    if( (fd = fopen(filename, "r"))==NULL) {
      err_prnt(filename);
      return NULL;
    }

  if (fread(&cvl_header, sizeof(struct cvl_p),1,fd) != 1) {
    err_prnt("reading header");
    return NULL;
  }
  if (cvl_header.cv_tag != 044520) {
    sprintf(s,"%s: Not cvl image (Tag: %lo)\n",strlen(filename)?
            filename:"'<stdin>'",cvl_header.cv_tag);
    err_prnt(s);
    return NULL;
  }
  rows = cvl_header.cv_rows;
  cols = cvl_header.cv_cols;
  *c = cols; *r = rows;
  if (cvl_header.cv_lcom == 0) return fd;
  if ((comment=(char *)malloc(cvl_header.cv_lcom+1))==NULL) {
    err_prnt("malloc");
    return NULL;
  }
  if (fread(comment,1,cvl_header.cv_lcom,fd) != cvl_header.cv_lcom) {
    err_prnt("reading comment field");
    return NULL;
  }
  comment[cvl_header.cv_lcom-1]=0;
  /* printf("Title: %s",comment); */
  free(comment);
  return fd;
}

/*****************************************************************************/
FILE *pgm_open(filename,c,r)
/*

  Open PGM raw input image file.

  Parameters :  filename - file name string
                c,r - image size

*/
char *filename; int *c,*r;
{
  FILE *fd=stdin;
  int max_val=-1;
  char pgm_line[256],magic[3];

  magic[0]='\0';
  *c=*r=-1;

  if (strlen(filename))
    if( (fd = fopen(filename, "r"))==NULL) {
      err_prnt(filename);
      return NULL;
    }

  while (fgets(pgm_line,255,fd) && pgm_line[0]=='#');
  sscanf(pgm_line,"%2s%5d%5d%5d",magic,c,r,&max_val);

  if (strncmp(magic,"P5",2)) {
    sprintf(s,"%s: wrong magic number: %s\n",strlen(filename)?
            filename:"'<stdin>'",magic);
    err_prnt(s);
    return 0;
  }

  if ((*c)==-1) {
    while (fgets(pgm_line,255,fd) && pgm_line[0]=='#');
    sscanf(pgm_line,"%5d%5d%5d",c,r,&max_val);
  }

  if ((*r)==-1) {
    while (fgets(pgm_line,255,fd) && pgm_line[0]=='#');
    sscanf(pgm_line,"%5d%5d",r,&max_val);
  }

  if (max_val==-1) {
    while (fgets(pgm_line,255,fd) && pgm_line[0]=='#');
    sscanf(pgm_line,"%5d",&max_val);
  }
 
  if (max_val<0 || max_val>255) {
    err_prnt(strlen(filename)?filename:"'<stdin>'");
    return 0;
  }

  cols = *c;
  rows = *r;

  return fd;
}

/*****************************************************************************/
FILE *vis_open(filename,c,r)
/*

  Open VISILOG input image file.

  Parameters :  filename - file name string
                c,r - image size

*/
char *filename; int *c,*r;
{
  FILE *fd=stdin;
  struct desc_im header;
  int i=0,j=0;

  if (strlen(filename))
    if( (fd = fopen(filename, "r"))==NULL) {
      err_prnt(filename);
      return NULL;
    }

  if (fread(&header, sizeof(struct desc_im),1,fd) != 1) {
    err_prnt("reading header");
    return NULL;
  }

  conv_long(&header.i_magic);
  conv_long(&header.i_dimx);
  conv_long(&header.i_dimy);
  conv_long(&header.i_code);
  conv_long(&header.hdr_sys);
  conv_long(&header.hdr_data);

  if (header.i_magic != I_MAGICD) {
    err_prnt(strlen(filename)?filename:"'<stdin>'");
    return NULL;
  }

  *c=header.i_dimx;
  *r=header.i_dimy;

  if (header.i_code!=I_CHR) {
    err_prnt(strlen(filename)?filename:"'<stdin>'");
    return NULL;
  }

  cols = *c;
  rows = *r;

  return fd;
}

/*****************************************************************************/
FILE *raw_open(filename,c,r)
/*

  Open raw input image file.

  Parameters :  filename - file name string
                c,r - image size

*/
char *filename;
int *c,*r;
{
  FILE *fd=stdin;

  if (strlen(filename))
    if( (fd = fopen(filename, "r"))==NULL) {
      err_prnt(filename);
      return NULL;
    }

  cols = *c = MAX_SIZE; /* Could something else be concluded here ?? */
  rows = *r = MAX_SIZE;

  return fd;
}
/*****************************************************************************/
FILE *cor_open(filename,c,r)
/*

  Open cordinate input file.

  Parameters :  filename - file name string
                c,r - image size

*/
char *filename;
int *c,*r;
{
  FILE *fd=stdin;

  if (strlen(filename))
    if( (fd = fopen(filename, "r"))==NULL) {
      err_prnt(filename);
      return NULL;
    }

  cols = *c = MAX_SIZE; /* Could something else be concluded here ?? */
  rows = *r = MAX_SIZE;

  return fd;
}

/*****************************************************************************/
char *read_image(fd)
/*

  Read input image (one byte per pixel).

  Parameters :  fd - file discriptor

*/
FILE *fd;
{
  int Nread,Nrequired;
  char *Imagebuff;

  Nrequired = cols * rows ;
  Imagebuff = (char *)malloc(Nrequired);
  if (!Imagebuff) {
    err_prnt("mallocing Imagebuff");
    return 0;
  }

  if (fread(Imagebuff,1,Nrequired,fd) != Nrequired) {
    err_prnt("reading image");
    return 0;
  }
  fclose(fd);
  return Imagebuff;
}

/*****************************************************************************/
int *read_ske_image(fd)
/*

  Read input image (four bytes per pixel).

  Parameters :  fd - file discriptor

*/
FILE *fd;
{
  int Nread,Nrequired,*Imagebuff;

  Nrequired = cols * rows ;

  if (!(Imagebuff = (int *)malloc(Nrequired*sizeof(int)))) {
    err_prnt("mallocing Imagebuff");
    return 0;
  }
  if (fread(Imagebuff,sizeof(int),Nrequired,fd) != Nrequired) {
    err_prnt("reading image");
    return 0;
  }
  fclose(fd);
  return Imagebuff;
}

/*****************************************************************************/
char *read_cor_image(fd)
/*

  Read input image (four bytes per pixel).

  Parameters :  fd - file discriptor

*/
FILE *fd;
{
  char *Imagebuff;
  int Nread,Nrequired,i;
  int x1, y1, cnt=0;
  char read_buff[250];
  Nrequired = cols * rows ;

  if (!(Imagebuff = (char *)malloc(Nrequired*sizeof(unsigned)))) {
    err_prnt("mallocing Imagebuff");
    return 0;
  }
	for(i=0;i<Nrequired;i++) Imagebuff[i]=(unsigned)0;

	while(fgets(read_buff, 50, fd)) {
		for(i=0;i<strlen(read_buff);i++) {
			if(read_buff[i]==' ' || read_buff[i]=='\n') {
				cnt++;
				read_buff[i]='\0';
				if(cnt==1) sscanf(read_buff,"%d", &x1);
				else {
					sscanf(read_buff,"%d", &y1);
					cnt=0;
					Imagebuff[y1*cols+x1]=
						(unsigned)255;
				}
				strcpy(read_buff,&read_buff[i+1]);
				i=0;	
			}	
		}	
	}


  fclose(fd);
  return Imagebuff;
}
 
/*****************************************************************************/
pic_from_disk(filename,pic,shrink,dim_x,dim_y)
/*

  Read image from disk into an array.

  Parameters :  filename - input image file name
                pic - image array
                shrink - shrink image verically by 1/3 (on/off)
                dim_x,dim_y - image size

*/
char *filename;
pic_type pic[][MAX_SIZE];
int shrink, *dim_x, *dim_y;
{
  FILE *image;
  int i,j,c,r,*im_int,*p_int,format;
  char *im,*p,name[256];

  if (filename==(char *)NULL)
    name[0]='\0';
  else
    strcpy(name,filename);

  if (strlen(name))
    switch(format=image_format(name)) {
      case CVL_PIC: if ((image = cvl_open(name,&c,&r))==NULL)
                      return 0;
                    break;
      case PGM_PIC: if ((image = pgm_open(name,&c,&r))==NULL)
                      return 0;
                    break;
      case VIS_PIC: if ((image = vis_open(name,&c,&r))==NULL)
                      return 0;
                    break;
      case SKE_PIC:
      case BIN_PIC: if ((image = raw_open(name,&c,&r))==NULL)
                      return 0;
                    break;
      case COR_PIC: if ((image = cor_open(name,&c,&r))==NULL)
                      return 0;
                    break;

      default:      return 0;
    }
  else
    if ((image = pgm_open(name,&c,&r))==NULL)
      return 0;

  if (c>MAX_SIZE || r>MAX_SIZE) {
    sprintf(s,"image too is wide, define BIG_IMAGES in Makefile ");
    strcat(s,"and recompile program");
    err_prnt(s);
    fclose(image);
    return 0;
  }

  for (i=0; i<r; i++)
    for (j=0; j<c; j++)
      pic[i][j]=(pic_type)0;
   
  if (format==SKE_PIC) { /* four bytes per pixel */
    if (!(im_int = read_ske_image(image))) return 0;
    p_int=im_int;

    if (shrink)
      for (i=0; i<r; i++)
        for(j=0; j<c; j++) {
          if (!pic[(int)(2*i/3)][j])
            pic[(int)(2*i/3)][j] = *p_int;
          p_int++;
        }
    else
      for (i=0; i<r; i++)
        for(j=0; j<c; j++) {
          pic[i][j] = *p_int;
          p_int++;
        }

    free(im_int);

  } else { /* one byte per pixel */
    if(format==COR_PIC) {
	if (!(im = read_cor_image(image))) return 0; 
    }
    else 
	if (!(im = read_image(image))) return 0;
    p=im;

    if (shrink)
      for(i=0;i<r;i++) 
        for(j=0;j<c;j++) {
          if (!pic[(int)(2*i/3)][j])
            pic[(int)(2*i/3)][j] = (unsigned char)(*p);
          p++;
        }
    else
      for(i=0;i<r;i++) 
        for(j=0;j<c;j++) {
          pic[i][j] = (unsigned char)(*p);
          p++;
        }

    free(im);
  }

  *dim_x=c; *dim_y=r;
  return 1;
}

/*****************************************************************************/
make_object_from_image(name,pic,shrink,dim_x,dim_y,noise)
/*

  Read image from disk into an array and add noise to it.

  Parameters :  name - input image file name
                pic - image array
                shrink - shrink image verically by 1/3 (on/off)
                dim_x,dim_y - image size
                noise - additive noise percentage

*/
char *name;
pic_type pic[][MAX_SIZE];
int shrink, *dim_x, *dim_y;
double noise;
{
  int i, j;

  if (!pic_from_disk(name,pic,shrink,dim_x,dim_y))
    return 0;

  for(i=0;i<(*dim_y);i++) 
    for(j=0;j<(*dim_x);j++)
      if (pic[i][j])
        pic[i][j] = (pic_type)OBJECT_PIX_VAL;
  if (noise>0.0 && noise<100.0)
    add_noise_to_pic(pic,noise,*dim_x,*dim_y,1,OBJECT_PIX_VAL);

  return 1;
}

/*****************************************************************************/
real_params_from_disk(filename,params)
/*

  Read parameter reference file into a array.

  Parameters :  filename - paramter file name
                params - parameter array

*/
char *filename;
double params[][2];
{
  FILE *fd;
  int i, num_of_params=-1, size_of_params;
  char str[256], magic[11];

  magic[0]='\0';
  params[0][0]=params[1][0]=params[1][1]=-1;

  if( (fd = fopen(filename, "r"))==NULL) {
    err_prnt(filename);
    return 0;
  }

  while (fgets(str,255,fd) && str[0]=='#'); /*reject comment lines */
  sscanf(str,"%11s%5d%lf%lf",magic,&num_of_params,&params[1][0],&params[1][1]);

  if (strncmp(magic,"LINE_PARAMS",11)) {
    err_prnt(filename);
    return 0;
  }

  if (num_of_params==-1) {
    while (fgets(str,255,fd) && str[0]=='#'); /* official comment field */
    sscanf(str,"%5d%lf%lf",&num_of_params,&params[1][0],&params[1][1]);
  }
  params[0][0]=(double)num_of_params;

  if (params[0][0]<0 || params[0][0]>MAX_LINES) {
    err_prnt(filename);
    return 0;
  }

  if (params[1][0]==-1) {
    while (fgets(str,255,fd) && str[0]=='#');
    sscanf(str,"%lf%lf",&params[1][0],&params[1][1]);
  }

  if (params[1][1]==-1) {
    while (fgets(str,255,fd) && str[0]=='#');
    sscanf(str,"%lf",&params[1][1]);
  }

  if (params[1][0]<0 || params[1][0]>100 ||
      params[1][1]<0 || params[1][1]>100) {
    err_prnt(filename);
    return 0;
  }

  size_of_params=(int)(params[0][0])+2;
  for (i=2; i<size_of_params; i++)
    if (fscanf(fd,"%lf%lf",&params[i][0],&params[i][1])!=2) {
      err_prnt(filename);
      return i-1;
    }

  return size_of_params-2;
}
