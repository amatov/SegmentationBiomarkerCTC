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
 * File:    ht_imgutils.c
 * Purpose: image processing subroutines
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"

/*****************************************************************************/
init_pic(pic, dim_x, dim_y)
/*

  Zero the image.

  Parameters :  pic - image
                dim_x,dim_y - image size

*/
pic_type pic[][MAX_SIZE];
int dim_x, dim_y;
{
  int i,j;

  for (i=0; i<dim_y; i++)
    for (j=0; j<dim_x; j++)
      pic[i][j]=(pic_type)0;
}

/*****************************************************************************/
copy_pic(pic, pic2, dim_x, dim_y)
/*

  Copy an image to an other.

  Parameters :  pic - image
                dim_x,dim_y - image size

*/
pic_type pic[][MAX_SIZE], pic2[][MAX_SIZE];
int dim_x, dim_y;
{
  int i,j;

  for (i=0; i<dim_y; i++)
    for (j=0; j<dim_x; j++)
      pic[i][j]=pic2[i][j];
}

/*****************************************************************************/
add_noise_to_pic(pic,Noise,dim_x,dim_y,amplitude,coef)
/*

  Add noise to image.

  Parameters :  pic - image
                Noise - noise percentage
                dim_x,dim_y - image size
                amplitude - noise amplitude
                coef - coefficient by witch the noise is multiplied

*/
pic_type pic[][MAX_SIZE],amplitude;
double Noise;
int dim_x,dim_y,coef;
{
  int i,j,x,y;
  long k,count=0L,num_of_noise_points;

  init_random_generator();
/*
  for(i=0; i<dim_y; i++)
    for(j=0; j<dim_x; j++)
      if (pic[i][j]==(pic_type)OBJECT_PIX_VAL)
        count++;

  num_of_noise_points=(long)((double)count/(1.0-Noise/100.0));
*/
  num_of_noise_points=(long)(dim_x*dim_y*Noise/100.0);

  for (k=0L; k<num_of_noise_points; k++) {
    x=(int)rnd(dim_x-1);
    y=(int)rnd(dim_y-1);
    pic[x][y]=(pic_type)(coef*rnd(amplitude));
  }
}

/*****************************************************************************/
order_12(x,y)
int *x,*y;
/*

  Add noise to image.

  Parameters :  pic - image
                Noise - noise percentage
                dim_x,dim_y - image size
                amplitude - noise amplitude
                coef - coefficient by witch the noise is multiplied

*/
{
  int tmp;

  if (*x > *y) {
    tmp = *x;
    *x = *y;
    *y = tmp;
  }
}

/*****************************************************************************/
put_line_to_pic(pic,dim_x,dim_y,s_x,s_y,e_x,e_y,val)
pic_type pic[][MAX_SIZE],val;
int dim_x,dim_y,s_x,s_y,e_x,e_y;
/*

  Draw a line to image space.

  Parameters :  pic - image
                dim_x,dim_y - image size
                s_x,s_y,e_x,e_y - line endpoints
                val - drawing gray level

*/
{
  int x, y;
  double a, b;

  if (s_x<0 || s_x>=dim_x || s_y<0 || s_y>=dim_y ||
      e_x<0 || e_x>=dim_x || e_y<0 || e_y>=dim_y)
    return; /* line is (partly) out of image */

  if (s_x==e_x) {
    x=s_x;
    if (s_y==e_y) {
      pic[s_y][s_x] = (pic_type)val;
      return;
    }
    order_12(&s_y,&e_y);
    y=s_y;
    do
      pic[y][x] = (pic_type)val;
    while (y++ < e_y);
  } else
    if (s_y==e_y) {
      y=s_y;
      if (s_x == e_x) {
        pic[s_y][s_x] = (pic_type)val;
        return;
      }
      order_12(&s_x,&e_x);
      x=s_x;
      do
        pic[y][x] = (pic_type)val;
      while (x++ < e_x);
    } else {
      a=(double)(e_y-s_y)/(double)(e_x-s_x);
      b=(double)s_y-a*s_x;
      if (fabs(a) <= 1.0) {
        order_12(&s_x,&e_x);
        x=s_x;
        do {
          y = nint(a*x+b);
          pic[y][x] = (pic_type)val;
        } while (x++ < e_x); 
      } else {
        order_12(&s_y,&e_y);
        y=s_y;
        do {
          x = nint((y - b)/a);
          pic[y][x] = (pic_type)val;
        } while (y++ < e_y);
      }
    }
}

/*****************************************************************************/
store_all_line_segments(sto_pic,dim_x,dim_y,start,val, params_out)
/*

  Draw lines to output image space.

  Parameters :  sto_pic - image
                dim_x,dim_y - image size
                start - line structure
                val - drawing gray level
		params_out - save file for parameters

*/
pic_type sto_pic[][MAX_SIZE],val;
int dim_x,dim_y;
line *start;
FILE *params_out;
{
  line *p=start;
  int i, j;

  while (p!=(line *)NULL) {
    put_line_to_pic(sto_pic,dim_x,dim_y,p->start_x,p->start_y,
                                        p->end_x,p->end_y,val);
    p=p->next;
  }
}

/*****************************************************************************/
store_all_maxs_line(gpic,dim_x,dim_y,maxs,n,val,ParametersOut,round,start)
/*

  Draw lines to image space.

  Parameters :  sto_pic - image
                dim_x,dim_y - image size
                maxs - list of accumulator maxima
                n - size of maximum list
                val - drawing gray level
		params_out - save file for parameters
		ParametersOut - output parameter file
		round - rounding value
		start - line structure

*/
pic_type gpic[][MAX_SIZE],val;
double maxs[][2],round;
int dim_x,dim_y,n;
char ParametersOut[256];
line *start;
{
  int i,s_x,s_y,e_x,e_y;
  double rho,theta;
  FILE *params_out=NULL;
  
  for (i=0; i<n; i++) {
    if (maxs[i][0]>0.0)
      if (isinf(maxs[i][0])) {
        s_x=e_x=nint(maxs[i][1]);
        s_y=0; e_y=dim_y-1;
      } else {
        s_x=(-maxs[i][1]/maxs[i][0])<=0?0:nint(-maxs[i][1]/maxs[i][0]);
        s_x=s_x>=dim_x?dim_x-1:s_x;
        s_y=maxs[i][1]<=0?0:nint(maxs[i][1]);
        s_y=s_y>=dim_y?dim_y-1:s_y;

        e_x=((dim_y-1-maxs[i][1])/maxs[i][0])<=0?0:
             nint((dim_y-1-maxs[i][1])/maxs[i][0]);
        e_x=e_x>=dim_x?dim_x-1:e_x;
        e_y=(maxs[i][0]*(dim_x-1)+maxs[i][1])<=0?0:
             nint(maxs[i][0]*(dim_x-1)+maxs[i][1]);
        e_y=e_y>=dim_y?dim_y-1:e_y;
      }
    else
      if (maxs[i][0]==0.0) {
        s_x=0; e_x=dim_x-1;
        s_y=e_y=nint(maxs[i][1]);
      } else
        if (isinf(maxs[i][0])) {
          s_x=e_x=nint(maxs[i][1]);
          s_y=0; e_y=dim_y-1;
        } else {
          s_x=(-maxs[i][1]/maxs[i][0])<=0?0:nint(-maxs[i][1]/maxs[i][0]);
          s_x=s_x>=dim_x?dim_x-1:s_x;
          s_y=(maxs[i][0]*(dim_x-1)+maxs[i][1])<=0?0:
               nint(maxs[i][0]*(dim_x-1)+maxs[i][1]);
          s_y=s_y>=dim_y?dim_y-1:s_y;

          e_x=((dim_y-1-maxs[i][1])/maxs[i][0])<=0?0:
               nint((dim_y-1-maxs[i][1])/maxs[i][0]);
          e_x=e_x>=dim_x?dim_x-1:e_x;
          e_y=maxs[i][1]<=0?0:nint(maxs[i][1]);
          e_y=e_y>=dim_y?dim_y-1:e_y;
        }
/*
    printf("a: %lf b: %lf sx: %d sy: %d ex: %d ey: %d\n",
           maxs[i][0],maxs[i][1],s_x,s_y,e_x,e_y);
*/  

    put_line_to_pic(gpic,dim_x,dim_y,s_x,s_y,e_x,e_y,val);
  }
/* Save parameters using a, b or rho, theta presentation */ 
 
     if (strlen(ParametersOut))
      if ((params_out=fopen(ParametersOut,"w"))==(FILE *)NULL) {
        err_prnt(ParametersOut);
      } else {
	write_params(start,params_out);
	if(check_file_name(ParametersOut,".rho")) {
	  fprintf(params_out,"rho\t\ttheta\n");
	  for (i=0; i<n; i++) {
	    theta=rint(round*atan2(1.0,-maxs[i][0]))/round;
      	    rho=rint(round*(maxs[i][1]/(sqrt(pow(maxs[i][0],2.0)+1))))/round;
	    fprintf(params_out,"%-12f %-12f\n",rho,theta);
	  }
	} else {
	  fprintf(params_out,"a\t\tb\n");
	  for (i=0; i<n; i++)
	    fprintf(params_out,"%-12f %-12f\n",maxs[i][0],maxs[i][1]);
	}
        fclose(params_out);
      }
}

/*****************************************************************************/
store_all_maxs_rho_theta(gpic,dim_x,dim_y,maxs_double,n,val,ParametersOut,start)
/*

  Draw lines to image space.

  Parameters :  gpic - image
                dim_x,dim_y - image size
                maxs_double - list of accumulator maxima
                n - size of maximum list
                val - drawing gray level
		ParametersOut - output parameter file
		start - line structure

*/
pic_type gpic[][MAX_SIZE],val;
int dim_x,dim_y,n;
double maxs_double[][2];
char ParametersOut[256];
line *start;
{
  int i,x,y;
  double b,k,rho,theta;
  FILE *params_out=NULL;

  for (i=0; i<n; i++)
    if (fabs(sin(maxs_double[i][1]))<0.001) {
      x=(int)maxs_double[i][0];
      if (x>=0 && x<dim_x) 
        for (y=0; y<dim_y; y++)
          gpic[y][x] = (pic_type)val;
    } else {
      k = -cos(maxs_double[i][1])/sin(maxs_double[i][1]);
      b = maxs_double[i][0]/sin(maxs_double[i][1]);

      if (fabs(k)<=1.0)
        for (x=0; x<dim_x; x++) {
          y=(int)(k*x + b);
          if (y>=0 && y<dim_y)
            gpic[y][x] = (pic_type)val;
         }
      else
        for (y=0; y<dim_y; y++) {
          x=(int)((y - b)/k);
          if (x>=0 && x<dim_x) 
            gpic[y][x] = (pic_type)val;
        }
	
    }

/* Save parameters using a, b or rho, theta presentation */  
 
    if (strlen(ParametersOut))
      if ((params_out=fopen(ParametersOut,"w"))==(FILE *)NULL) {
        err_prnt(ParametersOut);
      } else {
  	write_params(start,params_out);
	if(!check_file_name(ParametersOut,".rho")) { /* It WAS a bug */
          fprintf(params_out,"a\t\tb\n");	     /* in an old version here.*/
	  for (i=0; i<n; i++) {
	    theta=atan2(1.0,-maxs_double[i][0]);
      	    rho=maxs_double[i][1]/(sqrt(pow(maxs_double[i][0],2.0)+1));
	    fprintf(params_out,"%-12f %-12f\n",rho,theta);
	  }
	} else {
	  fprintf(params_out,"rho\t\ttheta\n");
	  for (i=0; i<n; i++) {
	    fprintf(params_out,"%-12f %-12f\n", 
		maxs_double[i][0],maxs_double[i][1]);
	  }
	}
        fclose(params_out);
      }
}

/*****************************************************************************/
write_params(start,params_out)
/*

  Write start and endpoints of lines to a file.

  Parameters :	start - line structure
		params_out - save file for parameters

*/
line *start;
FILE *params_out;
{
  line *p=start;
  int i, j;

  fprintf(params_out,"Found lines\n(X1, Y1)  -  (X2, Y2)\n");

  while (p!=(line *)NULL) {
    fprintf(params_out, "(%3d, %3d) - (%3d, %3d)\n", 
	p->start_x,p->start_y,p->end_x,p->end_y);
    p=p->next;
  }
}
