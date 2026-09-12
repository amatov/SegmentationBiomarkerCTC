/*****************************************************************************
 * Standard Hough Transform (SHT) and Randomized Hough Transform (RHT)
 * Programs for testing the efficiency of the HT and it's variants,
 * beta version
 *
 * Lappeenranta University of Technology, Department of Information Technology
 * Laboratory of Information Processing
 *
 * Authors: Petri Hirvonen, Jouni Ikonen (Jouni.Ikonen@lut.fi)
 *	    Pekka Kultanen, Heikki K{lvi{inen (Heikki.Kalviainen@lut.fi)
 *
 * File: sht_PSfuncs.c (functions for writing PS files (SHT))
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ht_hough.h"

int BoundingBox[]= {25,25,241,241};  /* 3.000 in x 3.000 in */

FILE *psfile;

float find_maxv(),find_minv();
double rint(),fabs(),sqrt(),sin(),cos(),atan2();

move_abs_2PS(x,y)
float x,y;
  {
    fprintf(psfile,"N\n");
    fprintf(psfile,"%6.2f %6.2f\n",x,y);
    fprintf(psfile,"M\n");
  }

line_abs_2PS(x,y)
float x,y;
  {
    fprintf(psfile,"%6.2f %6.2f\n",x,y);
    fprintf(psfile,"L S\n");
  }

put_pixelPS(i,j,val)
int i,j;
int val;
  {
    fprintf(psfile,"N %6.2f %6.2f P\n",
            j/256.0*(BoundingBox[2]-BoundingBox[0])+BoundingBox[0],
            (256.0-i)/256.0*(BoundingBox[3]-BoundingBox[1])+BoundingBox[1]);
  }

show_line_c2PS(x1,y1,dim_x,dim_y,i_dim,j_dim)
int x1,y1;
int dim_x,dim_y,i_dim,j_dim;
  {
    int x,y;
    double b,k,rho,theta,rho_max;

    rho_max=sqrt(2.0)*sqrt((double)dim_x*dim_x+(double)dim_y*dim_y);

    theta=(double)((y1-(double)j_dim/2.0)*(PI/2.0)/((double)j_dim/2.0));
    rho=(double)(x1*2.0*rho_max/i_dim-rho_max);
    if (fabs(sin(theta))<0.001)
      for (y=0; y<255; y++)
      {
        x=(int)rho;

        if (x>=0 && x<=256) 
          put_pixelPS(x,y,255);
      }
    else
      {
      k = -cos(theta)/sin(theta);
      b = rho/sin(theta);

      if (fabs(k)<=1.0)
        for (x=0; x<255; x++)
          {
            y=(int)(k*x + b);
            if (y>=0 && y<=255)
              put_pixelPS(x,y,255);
          }
      else
        for (y=0; y<255; y++)
          {
            x=(int)((y - b)/k);
            if (x>=0 && x<=255) 
              put_pixelPS(x,y,255);
          }
      }
  }

VisiToPSwithLines_c2(width,height,depth,pic,start,maxs,n)
int width,height,depth;
int pic[256][256];
line *start;
int maxs[][2];
int n;
  {
    int i,j;
    line *p=start;

    if ((psfile=fopen("greyline.ps","w"))==(FILE *)NULL)
      err_prnt("Cannot create PS-file.\n");

    BoundingBox[0]=25;
    BoundingBox[1]=25;
    BoundingBox[2]=205;
    BoundingBox[3]=205;

    fprintf(psfile,"%%%%BoundingBox: %d %d %d %d\n",
		   BoundingBox[0],BoundingBox[1],BoundingBox[2],BoundingBox[3]);
 
    fprintf(psfile,"\n%%begin(plot)\n\n");

    fprintf(psfile,"\n%%macros for newpath, moveto, and {lineto} repeat stroke\n");
    fprintf(psfile,"/N {newpath} def /M {moveto}  def /L {lineto} def /S {stroke} def\n\n");
    fprintf(psfile,"/MR {rmoveto}  def /LR {rlineto} def /P {moveto currentpoint lineto stroke} def\n\n");
    fprintf(psfile,"\n.5 setlinewidth   1 setlinecap   1 setlinejoin\n\n");

    fprintf(psfile,"\ngsave\n");

    fprintf(psfile,"/picstr %d string def\n",width/2);
    fprintf(psfile,"25 25 translate\n");

    fprintf(psfile,"%d %d scale\n",180,180);

    fprintf(psfile,"%d %d %d\n",width/2,height/2,depth);
    fprintf(psfile,"[%d 0 0 %d 0 %d]\n",width/2,-height/2,height/2);
    fprintf(psfile,"{currentfile picstr readhexstring pop}\n");
    fprintf(psfile,"image\n");

    for (i=0; i<width; i+=2)
      for (j=0; j<height; j+=2)
        {
        if (!(j % 30)) fprintf(psfile,"\n");
        fprintf(psfile,"%02X",pic[i][j]);
        }

    fprintf(psfile,"\n\ngrestore\n");

/*    fprintf(psfile,"\n0 setgray\n");

    while (p!=(line *)NULL)
      {
      move_abs_2PS((p->start_y/256.0*(BoundingBox[2]-BoundingBox[0])+BoundingBox[0]),
		   ((255.0-p->start_x)/256.0*(BoundingBox[3]-BoundingBox[1])+BoundingBox[1]));
      line_abs_2PS((p->end_y/256.0*(BoundingBox[2]-BoundingBox[0])+BoundingBox[0]),
		   ((255.0-p->end_x)/256.0*(BoundingBox[3]-BoundingBox[1])+BoundingBox[1]));
      p=p->next;
      } */
    fprintf(psfile,"\n0 setgray\n");
    for (i=0; i<n; i++) show_line_c2PS(maxs[i][0],maxs[i][1],
                                       256,256,256,256);
/* BoundingBox:n piirto debuggausta varten */

    fprintf(psfile,"\nN\n");
    fprintf(psfile,"%d %d\n",BoundingBox[0],BoundingBox[1]);
    fprintf(psfile,"M\n");
    fprintf(psfile,"%d %d\n",BoundingBox[2],BoundingBox[1]);
    fprintf(psfile,"L\n");
    fprintf(psfile,"%d %d\n",BoundingBox[2],BoundingBox[3]);
    fprintf(psfile,"L\n");
    fprintf(psfile,"%d %d\n",BoundingBox[0],BoundingBox[3]);
    fprintf(psfile,"L\n");
    fprintf(psfile,"%d %d\n",BoundingBox[0],BoundingBox[1]);
    fprintf(psfile,"L S\n");

    fprintf(psfile,"\n\n%%end(plot)\n");

    fprintf(psfile,"\n/#copies 1 def\n");
    fprintf(psfile,"showpage\n");

    fclose(psfile);

    return(0);
  }
