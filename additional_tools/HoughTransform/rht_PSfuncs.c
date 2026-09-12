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
 * File:    rht_PSfuncs.c
 * Purpose: writing PS files of dynamic accumulator
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"
#define LARGE 0
#define SMALL 1

int BoundingBox[]= {25,25,241,241};  /* 3.000 in x 3.000 in */

FILE *psfile;

float find_maxv(),find_minv();

/*****************************************************************************/
DrawAccumulatorSpacePS(AccuSpace,header,x_label,y_label,xmin,xmax,ymin,ymax,
                       accumax)
/*

  Make a PostScript file of the contents of the accumulator.

  Parameters :  AccuSpace - accumulator
                header - heading text string
                x_label,y_label - axis labels
                xmin,xmax,ymin,ymax - current parameter ranges
                accumax - maximum value of the accumulator

*/
InfMat *AccuSpace;
char *header,*x_label,*y_label;
double xmin,xmax,ymin,ymax;
int accumax;
{
  if (!init_PS())
    return 0;
  headingPS(header);
  axis_x_abs_2PS(x_label,(float)xmin,(float)xmax);
  axis_y_abs_2PS(y_label,(float)ymin,(float)ymax);
  draw_dataPS(AccuSpace,(float)xmin,(float)xmax,
                        (float)ymin,(float)ymax,accumax);
  close_PS();
}

/*****************************************************************************/
DrawLineDataPS(x,y,n,header,x_label,y_label)
/*

  Make a PostScript file of the contents of the accumulator during the detection
  process.

  Parameters :  x,y - data vectors
                n - size of the data vectors
                header - heading text string
                x_label,y_label - axis labels

*/
float x[],y[];
int n;
char *header,*x_label,*y_label;
{
  float xmin,xmax,ymin,ymax,max;

  xmin=find_minv(x,n);
  xmax=find_maxv(x,n);
  ymin=find_minv(y,n);
  ymax=find_maxv(y,n);

  BoundingBox[0]=25;
  BoundingBox[1]=25;
  BoundingBox[2]=265;
  BoundingBox[3]=200;

  if (!init_PS())
    return 0;
  headingPS(header);
  axis_x_abs_2PS(x_label,xmin,xmax);
  axis_y_abs_2PS(y_label,ymin,ymax);
  draw_line_dataPS(x,y,n,xmin,xmax,ymin,ymax);
  close_PS();
}

/*****************************************************************************/
draw_dataPS(AccuSpace,min_val_x,max_val_x,min_val_y,max_val_y,accumax)
/*

  Make a PostScript file of the contents of the accumulator.

  Parameters :  AccuSpace - accumulator
                min_val_x,max_val_x,min_val_y,max_val_y - current parameter ranges
                accumax - maximum value of the accumulator

*/
InfMat *AccuSpace;
float max_val_x,min_val_x,max_val_y,min_val_y;
int accumax;
{
  InfIndex *p=AccuSpace;
  Accu *s;
  float orig_x,orig_y,x,y,x_len,y_len,scale_x,scale_y;

  orig_x=(BoundingBox[2]-BoundingBox[0])*0.1+BoundingBox[0];
  orig_y=(BoundingBox[3]-BoundingBox[1])*0.1+BoundingBox[1];

  x_len=(BoundingBox[2]-BoundingBox[0])*0.9-(BoundingBox[2]-BoundingBox[0])*0.1;
  y_len=(BoundingBox[3]-BoundingBox[1])*0.9-(BoundingBox[3]-BoundingBox[1])*0.1;

  scale_x=x_len/(max_val_x-min_val_x);
  scale_y=y_len/(max_val_y-min_val_y);

  if (p->NextIndex==(InfIndex *)NULL) {
    err_prnt("Impossible to draw NULL or infinite accumulator.\n");
    return;
  }

  p=p->NextIndex;

  s=p->Accu; 

  do {
    if ( (p->X_coord>min_val_x) && (p->X_coord<max_val_x) ) {
      s=p->Accu;
      do {
        if ((s->Y_coord>min_val_y) && (s->Y_coord<max_val_y)) {
          x=(float)((p->X_coord-min_val_x)*scale_x+orig_x);
          y=(float)((s->Y_coord-min_val_y)*scale_y+orig_y);
          move_abs_2PS(x,y); 
          set_line_indexPS(255/*(int)rint(127.0*s->Accumulator/accumax)+125*/); 
          line_abs_2PS(x,y);
        }
        s=s->Up;
      } while (s);
    }
    p=p->NextIndex;
  } while (p);
}

/*****************************************************************************/
draw_line_dataPS(x,y,n,min_val_x,max_val_x,min_val_y,max_val_y)
/*

  Make a PostScript file of the contents of the accumulator during the detection
  process.

  Parameters :  x,y - data vectors
                n - size of the data vectors
                min_val_x,max_val_x,min_val_y,max_val_y - parameter ranges

*/
float x[],y[];
int n;
float max_val_x,min_val_x,max_val_y,min_val_y;
{
  float umin,umax,vmin,vmax,orig_x,orig_y,
        x_len,y_len,scale_x,scale_y;
  int i;

  orig_x=(BoundingBox[2]-BoundingBox[0])*0.1+BoundingBox[0];
  orig_y=(BoundingBox[3]-BoundingBox[1])*0.1+BoundingBox[1];

  x_len=(BoundingBox[2]-BoundingBox[0])*0.9-(BoundingBox[2]-BoundingBox[0])*0.1;
  y_len=(BoundingBox[3]-BoundingBox[1])*0.9-(BoundingBox[3]-BoundingBox[1])*0.1;
     
  scale_x=x_len/(max_val_x-min_val_x);
  scale_y=y_len/(max_val_y-min_val_y);

  set_line_indexPS(255); 
    
  for(i=1;i<n; i++) {
    move_abs_2PS((x[i-1]-min_val_x)*scale_x+orig_x,
                 (y[i-1]-min_val_y)*scale_y+orig_y);
    line_abs_2PS((x[i]-min_val_x)*scale_x+orig_x,
                 (y[i]-min_val_y)*scale_y+orig_y);
  }
}

/*****************************************************************************/
headingPS(name)
/*

  Print heading string to the PostScript file.

  Parameters :  name - heading string

*/
char *name;
{
  float label_x,label_y;

  label_x=(BoundingBox[2]-BoundingBox[0])*0.5+BoundingBox[0];
  label_y=(BoundingBox[3]-BoundingBox[1])*0.91+BoundingBox[1];

  textPS(label_x,label_y,name,LARGE);
}

/*****************************************************************************/
axis_x_abs_2PS(name,min_val,max_val)
/*

  Print x-axis label string and draw sticks.

  Parameters :  name - label string
                min_val,max_val - parameter ranges

*/
char *name;
float min_val, max_val;
{
  float orig_x,orig_y,end_x,label_x,label_y,
        label_inc,inc,stick_y,stick_height,label_shift;
  char number[20];
  int i;

  orig_x=(BoundingBox[2]-BoundingBox[0])*0.1+BoundingBox[0];
  orig_y=(BoundingBox[3]-BoundingBox[1])*0.1+BoundingBox[1];

  end_x=(BoundingBox[2]-BoundingBox[0])*0.9+BoundingBox[0];

  move_abs_2PS(orig_x,orig_y);
  line_abs_2PS(end_x,orig_y);

  label_x=(BoundingBox[2]-BoundingBox[0])*0.75+BoundingBox[0];
  label_y=(BoundingBox[3]-BoundingBox[1])*0.01+BoundingBox[1];

  x_axis_textPS(label_x,label_y,name,LARGE);

  label_x=orig_x;
  label_y=(BoundingBox[3]-BoundingBox[1])*0.05+BoundingBox[1];

  stick_y=(BoundingBox[3]-BoundingBox[1])*0.09+BoundingBox[1];
  stick_height=(BoundingBox[3]-BoundingBox[1])*0.02;

  inc=(max_val-min_val)/10.0;
  label_inc=(end_x-orig_x)/10.0;
  label_shift=label_inc*0.2;

  sprintf(number,"%-6.2f",min_val);
  x_axis_textPS(label_x/*-label_shift*/,label_y,number,SMALL);

  move_abs_2PS(label_x,stick_y);
  line_rel_2PS(0.0,stick_height);

  for(i=0; i<10; i+=2) {
    label_x+=(2*label_inc);
    sprintf(number,"%-6.2f",min_val+=(2*inc));
    x_axis_textPS(label_x/*-label_shift*/,label_y,number,SMALL);
    move_abs_2PS(label_x,stick_y);
    line_rel_2PS(0.0,stick_height);
  } 
}

/*****************************************************************************/
axis_y_abs_2PS(name,min_val,max_val)
/*

  Print y-axis label string and draw sticks.

  Parameters :  name - label string
                min_val,max_val - parameter ranges

*/
char *name;
float min_val, max_val;
{
  float orig_x,orig_y,end_y,label_x,label_y,
        label_inc,inc,stick_x,stick_height,label_shift;
  char number[20];
  int i;

  orig_x=(BoundingBox[2]-BoundingBox[0])*0.1+BoundingBox[0];
  orig_y=(BoundingBox[3]-BoundingBox[1])*0.1+BoundingBox[1];

  end_y=(BoundingBox[3]-BoundingBox[1])*0.9+BoundingBox[1];

  move_abs_2PS(orig_x,orig_y);
  line_abs_2PS(orig_x,end_y);

  label_x=(BoundingBox[2]-BoundingBox[0])*0.04+BoundingBox[0];
  label_y=(BoundingBox[3]-BoundingBox[1])*0.75+BoundingBox[1];

  y_axis_textPS(label_x,label_y,name,LARGE);

  label_y=orig_y;
  label_x=(BoundingBox[2]-BoundingBox[0])*0.07+BoundingBox[0];

  stick_x=(BoundingBox[2]-BoundingBox[0])*0.09+BoundingBox[0];
  stick_height=(BoundingBox[3]-BoundingBox[1])*0.02;

  inc=(max_val-min_val)/10.0;
  label_inc=(end_y-orig_y)/10.0;
  label_shift=label_inc*0.2;

  sprintf(number,"%-6.2f",min_val);
  y_axis_textPS(label_x,label_y /*-label_shift */,number,SMALL);

  move_abs_2PS(stick_x,label_y);
  line_rel_2PS(stick_height,0.0);

  for(i=0; i<10; i+=2) {
    label_y+=(2*label_inc);
    sprintf(number,"%-6.2f",min_val+=(2*inc));
    y_axis_textPS(label_x,label_y/*-label_shift*/,number,SMALL);

    move_abs_2PS(stick_x,label_y);
    line_rel_2PS(stick_height,0.0);
  } 
}


/*****************************************************************************/
init_PS()
/*

  Initialize PostScript file.

  Parameters :

*/
{
  static int n=0;
  char name[20];

  n++; 
  sprintf(name,"hough%d.ps",n);
  if ((psfile=fopen(name,"w"))==(FILE *)NULL) {
    err_prnt("Cannot create PS-file.\n");
    return FALSE;
  }

  fprintf(psfile,"%%!PS-Adobe-1.0\n");
  fprintf(psfile,"%%%%Creator: Houghtool\n");
  fprintf(psfile,"%%%%Title: Houghtool Graph\n");
  fprintf(psfile,"%%%%BoundingBox: %d %d %d %d\n",
	  BoundingBox[0],BoundingBox[1],BoundingBox[2],BoundingBox[3]);
  fprintf(psfile,"%%%%EndComments\n\n");

  fprintf(psfile,"erasepage\n\n");

  fprintf(psfile,"%%begin(plot)\n\n");

  fprintf(psfile,"%% procedure to get translation in inches\n");
  fprintf(psfile,"/unscale {scale_pic div 72 mul} def\n\n");

  fprintf(psfile,"%% set scaling, rotation and translation parameters\n");
  fprintf(psfile,"/rot_pic 0.0 def\n");
  fprintf(psfile,"/scale_pic 1.0 def\n"); /*.5706 */
  fprintf(psfile,"/tranx_pic 0.0 def    /trany_pic 0.0 def\n");
  fprintf(psfile,"/font_spec {/Times-Roman findfont} def\n\n");

  fprintf(psfile,"%% control stuff for printing two pictures per page\n");
  fprintf(psfile,"%%/move_up {0.0 unscale 5.15 unscale translate} def\n");
  fprintf(psfile,"%%/move_down {0.0 unscale -5.15 unscale translate} def\n");
  fprintf(psfile,"/switch 1 def\n");
  fprintf(psfile,"%%/print_pic {switch 0 eq \n");
  fprintf(psfile,"%%      {move_down /switch 1 def} \n");
  fprintf(psfile,"%%      {copypage erasepage move_up /switch 0 def}\n"); 
  fprintf(psfile,"%%      ifelse } def\n\n");

  fprintf(psfile,"%% last page print control\n");
  fprintf(psfile,
          "/print_last_page {switch 1 eq {copypage /switch 0 def} if} def\n\n");

  fprintf(psfile,"%% rotate, scale, and translate picture\n");
  fprintf(psfile,"%%rot_pic rotate\n");
  fprintf(psfile,"%%scale_pic scale_pic scale\n");
  fprintf(psfile,"%%tranx_pic unscale trany_pic unscale translate\n\n");

  fprintf(psfile,"%% define font and font sizes for standard text and ");
  fprintf(psfile,"numbers and exponential\n%% numbers.\n");
  fprintf(psfile,"/largefont font_spec 1.000000e+001 scalefont def\n");
  fprintf(psfile,"/smallfont font_spec 0.700000e+001 scalefont def\n\n");

  fprintf(psfile,"%% line type setup\n");
  fprintf(psfile,"/solid   { []        0 setdash } def\n");
  fprintf(psfile,"/dotted  { [0 4]     0 setdash } def\n");
  fprintf(psfile,"/dashed  { [4]       0 setdash } def\n");
  fprintf(psfile,"/dotdash { [0 4 3 4] 0 setdash } def\n\n");

  fprintf(psfile,"/centre {stringwidth pop 2 div neg 0 rmoveto} def\n");
  fprintf(psfile,"/right  {stringwidth pop neg 0 rmoveto} def\n");
  fprintf(psfile,"/centre_rot90 \n");
  fprintf(psfile,"{translate 90 rotate 0 0 moveto centre show 0 0 moveto ");
  fprintf(psfile,"-90 rotate translate} def\n");
  fprintf(psfile,"%%line width, line cap, and joint spec\n");
  fprintf(psfile,".5 setlinewidth   1 setlinecap   1 setlinejoin\n\n");

  fprintf(psfile,"%%macros for newpath, moveto, and {lineto} repeat stroke\n");
  fprintf(psfile,"/N {newpath} def /M {moveto}  def ");
  fprintf(psfile,"/L {lineto} def /S {stroke} def\n");
  fprintf(psfile,"/MR {rmoveto}  def /LR {rlineto} def ");
  fprintf(psfile,"/P {moveto currentpoint lineto stroke} def\n\n");
  /* Draw BoundingBox for debugging */
  fprintf(psfile,"N\n");
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
  return TRUE;
}

/*****************************************************************************/
close_PS()
/*

  Close PostScript file.

  Parameters :

*/
{
  fprintf(psfile,"\n%%end(plot)\n\n");
  fprintf(psfile,"%% print page\n");
  fprintf(psfile,"%%print_pic\n\n");
  fprintf(psfile,"%% print last page if neccessary\n");
  fprintf(psfile,"print_last_page\n\n");
  fclose(psfile);
}

/*****************************************************************************/
move_rel_2PS(x,y)
/*

  Move to new relative location.

  Parameters : x,y - new location

*/
float x,y;
{
  fprintf(psfile,"N\n");
  fprintf(psfile,"%6.2f %6.2f\n",x,y);
  fprintf(psfile,"MR\n");
}

/*****************************************************************************/
line_rel_2PS(x,y)
/*

  Draw line to new relative location.

  Parameters : x,y - new location

*/
float x,y;
{
  fprintf(psfile,"%6.2f %6.2f\n",x,y);
  fprintf(psfile,"LR S\n");
}

/*****************************************************************************/
set_line_indexPS(n)
/*

  Set new line index (gray level).

  Parameters : n - index

*/
int n;
{
  fprintf(psfile,"%g setgray\n",(255.0-n)/255.0);
}

/*****************************************************************************/
textPS(x,y,text,size)
/*

  Print text to PS file.

  Parameters : x,y - text location
               text - text string
               size - font size

*/
float x,y;
char *text;
int size;
{
  if (size==LARGE)
    fprintf(psfile,"largefont setfont\n");
  else
    fprintf(psfile,"smallfont setfont\n");
  fprintf(psfile,"%6.2f %6.2f M\n",x,y);
  fprintf(psfile,"(%s)\n",text);
  fprintf(psfile,"centre\n");
  fprintf(psfile,"(%s) show\n",text);
}

/*****************************************************************************/
x_axis_textPS(x,y,text,size)
/*

  Print x-axis text to PS file.

  Parameters : x,y - text location
               text - text string
               size - font size

*/
float x,y;
char *text;
int size;
{
  if (size==LARGE)
    fprintf(psfile,"largefont setfont\n");
  else
    fprintf(psfile,"smallfont setfont\n");
  fprintf(psfile,"%6.2f %6.2f M\n",x,y);
  fprintf(psfile,"(%s)\n",text);
  fprintf(psfile,"centre\n");
  fprintf(psfile,"(%s) show\n",text);
}

/*****************************************************************************/
y_axis_textPS(x,y,text,size)
/*

  Print y-axis text to PS file.

  Parameters : x,y - text location
               text - text string
               size - font size

*/
float x,y;
char *text;
int size;
{
  if (size==LARGE)
    fprintf(psfile,"largefont setfont\n");
  else
    fprintf(psfile,"smallfont setfont\n");
  fprintf(psfile,"%6.2f %6.2f\n",-x,-y);
  fprintf(psfile,"(%s)\n",text);
  fprintf(psfile,"(%s)\n",text);
  fprintf(psfile,"%6.2f %6.2f\n",x,y);
  fprintf(psfile,"centre_rot90\n");
}

/*****************************************************************************/
move_abs_2PS(x,y)
/*

  Move to new absolute location.

  Parameters : x,y - new location

*/
float x,y;
{
  fprintf(psfile,"N\n");
  fprintf(psfile,"%6.2f %6.2f\n",x,y);
  fprintf(psfile,"M\n");
}

/*****************************************************************************/
line_abs_2PS(x,y)
/*

  Draw line to new absolute location.

  Parameters : x,y - new location

*/
float x,y;
{
  fprintf(psfile,"%6.2f %6.2f\n",x,y);
  fprintf(psfile,"L S\n");
}

/*****************************************************************************/
put_pixelPS(i,j,val)
/*

  Draw a point to PS file.

  Parameters : i,j, - point location
               val - point  gray level (not used)

*/
int i,j;
int val;
{
  fprintf(psfile,"N %6.2f %6.2f P\n",
          j/256.0*(BoundingBox[2]-BoundingBox[0])+BoundingBox[0],
          (256.0-i)/256.0*(BoundingBox[3]-BoundingBox[1])+BoundingBox[1]);
}

/*****************************************************************************/
float find_minv(x,n)
/*

  Find minimum value of the vector.

  Parameters : x - vector
               n - vector size

*/
float x[];
int n;
{
  int i;
  float fmin=1.0e10;

  for(i=0; i<n; i++)
    if (x[i]<fmin) fmin=x[i];

  return fmin;
}

/*****************************************************************************/
float find_maxv(x,n)
/*

  Find maximum value of the vector.

  Parameters : x - vector
               n - vector size

*/
float x[];
int n;
{
  int i;
  float fmax= -1.0e10;

  for(i=0; i<n; i++)
    if (x[i]>fmax) fmax=x[i];

  return fmax;
}
