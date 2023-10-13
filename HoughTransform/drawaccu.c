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
 * File:    drawaccu.c
 * Purpose: dynamic accmulator visualization functions, used with the XHoughtool
 *          only
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"
#include "xhoughtool.h"

float find_minv(), find_maxv();

/*****************************************************************************/
DrawAccumulatorSpace_graph(graph_quad,AccuSpace,header,x_label,y_label,
                           xmin,xmax,ymin,ymax,accumax)
quadrant graph_quad; /* drawing area */
InfMat *AccuSpace; /* dynamic accumulator */
char *header,*x_label,*y_label; /* labels */
double xmin,xmax,ymin,ymax; /* current accumulator space ranges */
int accumax; /* current accumulator maximum */
{
  ClearQuadrant_graph(graph_quad);

  heading(graph_quad,header);
  axis_x_abs_2(graph_quad,x_label,(float)xmin,(float)xmax); /* axis labels */
  axis_y_abs_2(graph_quad,y_label,(float)ymin,(float)ymax);
  draw_data(graph_quad,AccuSpace,
            (float)xmin,(float)xmax,(float)ymin,(float)ymax,accumax);
  FlushDisp();
}

/*****************************************************************************/
DrawLineData_graph(graph_quad,x,y,n,header,x_label,y_label)
quadrant graph_quad; /* drawing area */
float x[],y[]; /* data vectors */
int n; /* data vectors lenght */
char *header,*x_label,*y_label; /* labels */
{
  float xmin,xmax,ymin,ymax,max;

  ClearQuadrant_graph(graph_quad);

  xmin=find_minv(x,n);
  xmax=find_maxv(x,n);
  ymin=find_minv(y,n);
  ymax=find_maxv(y,n);

  if (ymin==ymax)
    ymax++;

  heading(graph_quad,header);
  axis_x_abs_2(graph_quad,x_label,xmin,xmax);
  axis_y_abs_2(graph_quad,y_label,ymin,ymax);
  draw_line_data(graph_quad,x,y,n,xmin,xmax,ymin,ymax);
  FlushDisp();
}

/*****************************************************************************/
DrawAccumulatorHistogram_graph(graph_quad,hist,size,header,x_label,y_label)
quadrant graph_quad; /* drawing area */
long hist[]; /* histogram */ 
int size; /* histogram lenght */
char *header,*x_label,*y_label; /* labels */
{
  int i;
  long maxval=0L;

  ClearQuadrant_graph(graph_quad);

  for(i=0; i<size; i++)
    if (hist[i]>maxval)
      maxval=hist[i];

  heading(graph_quad,header);
  axis_x_abs_2_hist(graph_quad,x_label,size);
  axis_y_abs_2(graph_quad,y_label,0.0,(float)maxval);
  draw_histogram(graph_quad,hist,size,maxval);
  FlushDisp();
}

/*****************************************************************************/
DrawHistoValues_graph(graph_quad,hist,wtime,AccuSpace,a_range,size,values,nhisto,
                      diff,header,x_label,y_label)
quadrant graph_quad; /* drawing area */
long hist[]; /* histogram */
InfMat3D *AccuSpace; /* dynamic accumulator */
int wtime,a_range,size,values,nhisto,diff; /* histogram specifications */
char *header,*x_label,*y_label; /* labels */
{
  int i,i2,histovalues=values;
  float minprob=100.0,maxprob=0.0,prob;

  ClearQuadrant_graph(graph_quad);

  if (nhisto)  
    for(i=size-1;(i>size-histovalues-1) && (i>size-nhisto-1) && (i>-1);i--) {
      i2=i-diff;
      if (hist[i2] && i2>=0) {
        prob=100.*(i+1)/(float)wtime;
        if (prob<minprob)
          minprob=prob;
        if (prob>maxprob)
          maxprob=prob;
      } else
        histovalues++;
    }
  else 
    for(i=size-1;(i>size-histovalues-1) && (i>-1);i--) {
      i2=i-diff;
      if (hist[i2] && i2>=0) {
        prob=100.*(i+1)/(float)wtime;
        if (prob<minprob)
          minprob=prob;
        if (prob>maxprob)
          maxprob=prob;
      } else
        histovalues++;
    }

  heading(graph_quad,header);
  axis_x_abs_2(graph_quad,x_label,0.0,(float)a_range);
  axis_y_abs_2(graph_quad,y_label,minprob,maxprob);
  draw_histovalues(graph_quad,hist,wtime,AccuSpace,size,values,nhisto,diff,
                   0.0,(float)a_range,minprob,maxprob);
  FlushDisp();
}

/*****************************************************************************/
draw_data(graph_quad,AccuSpace,min_val_x,max_val_x,min_val_y,max_val_y,accumax)
quadrant graph_quad; /* drawing area */
InfMat *AccuSpace; /* dynamic accumulator */
float max_val_x,min_val_x,max_val_y,min_val_y; /* current accumulator ranges */
int accumax; /* current accumulator maximum */
{
  InfIndex *p=AccuSpace;
  Accu *s;
  float umin,umax,vmin,vmax,orig_x,orig_y,
        x,y,x_len,y_len,scale_x,scale_y;
  int i;

  umin=0.0; /*(float)graph_quad.ulc_x;*/
  vmin=5.0; /*(float)graph_quad.ulc_y;*/
  umax=255.0; /*(float)graph_quad.lrc_x;*/
  vmax=200.0; /*(float)graph_quad.lrc_y;*/

  orig_x=umin+35.0;
  orig_y=vmax-15.0;

  x_len=umax-umin-orig_x-15.0;
  y_len=vmax-vmin-15.0-15.0;

  scale_x=x_len/(max_val_x-min_val_x);
  scale_y=y_len/(max_val_y-min_val_y);

  if (p->NextIndex==(InfIndex *)NULL)
    return;

  p=p->NextIndex;

  s=p->Accu; 

  do { /* finite values only */
    if ( (p->X_coord>min_val_x) && (p->X_coord<max_val_x) ) {
      s=p->Accu;
      do {
        if ((s->Y_coord>min_val_y) && (s->Y_coord<max_val_y)) {
          x=(float)((p->X_coord-min_val_x)*scale_x+orig_x);
          y=(float)(-((s->Y_coord-min_val_y)*scale_y)+orig_y);

          PutPixel_graph(graph_quad,nint(x),nint(y));
        }
        s=s->Up;
      } while (s);
    }
    p=p->NextIndex;
  } while (p);
}

/*****************************************************************************/
draw_line_data(graph_quad,x,y,n,min_val_x,max_val_x,min_val_y,max_val_y)
quadrant graph_quad; /* drawing area */
float x[],y[]; /* data vectors */
int n; /* data vectors lenght */
float max_val_x,min_val_x,max_val_y,min_val_y; /* current accumulator ranges */
{
  float umin,umax,vmin,vmax,orig_x,orig_y,
        x_len,y_len,scale_x,scale_y;
  int i;

  umin=0.0; /*(float)graph_quad.ulc_x;*/
  vmin=5.0; /*(float)graph_quad.ulc_y;*/
  umax=255.0; /*(float)graph_quad.lrc_x;*/
  vmax=200.0; /*(float)graph_quad.lrc_y;*/

  orig_x=umin+35.0;             /* (umax-umin)*0.1+umin; */
  orig_y=vmax-15.0;             /* (vmax-vmin)*0.1+vmin; */

  x_len=umax-umin-orig_x-15.0;  /* (umax-umin)*0.9-(umax-umin)*0.1; */
  y_len=vmax-vmin-15.0-15.0;    /* (vmax-vmin)*0.9-(vmax-vmin)*0.1; */

  scale_x=x_len/(max_val_x-min_val_x);
  scale_y=y_len/(max_val_y-min_val_y);

  for(i=1; i<n; i++)
  PutVector_graph(graph_quad, nint((x[i-1]-min_val_x)*scale_x+orig_x),
                              nint(-((y[i-1]-min_val_y)*scale_y)+orig_y),
                              nint((x[i]-min_val_x)*scale_x+orig_x),
                              nint(-((y[i]-min_val_y)*scale_y)+orig_y));
}

/*****************************************************************************/
draw_histogram(graph_quad,hist,size,maxval)
quadrant graph_quad; /* drawing area */
long hist[]; /* histogram */
int size; /* histogram size */
long maxval; /* current accumulator maximum */
{
  float umin,umax,vmin,vmax,orig_x,orig_y,x_len,y_len,scale_x,scale_y;
  int i;

  umin=0.0; /*(float)graph_quad.ulc_x;*/
  vmin=5.0; /*(float)graph_quad.ulc_y;*/
  umax=255.0; /*(float)graph_quad.lrc_x;*/
  vmax=200.0; /*(float)graph_quad.lrc_y;*/

  orig_x=umin+35.0;
  orig_y=vmax-15.0;

  x_len=umax-umin-orig_x-15.0;
  y_len=vmax-vmin-15.0-15.0;

  scale_x=x_len/(float)size;
  scale_y=y_len/((float)maxval-0.0);

  for (i=0; i<size; i++) {
    PutVector_graph(graph_quad, nint(orig_x+i*scale_x),
                                nint(orig_y-(float)hist[i]*scale_y),
                                nint(orig_x+i*scale_x),
                                nint(orig_y));
    PutVector_graph(graph_quad, nint(orig_x+i*scale_x),
                                nint(orig_y-(float)hist[i]*scale_y),
                                nint(orig_x+(i+1)*scale_x),
                                nint(orig_y-(float)hist[i]*scale_y));
    PutVector_graph(graph_quad, nint(orig_x+(i+1)*scale_x),
                                nint(orig_y-(float)hist[i]*scale_y),
                                nint(orig_x+(i+1)*scale_x),
                                nint(orig_y));
  }

}

/*****************************************************************************/
draw_histovalues(graph_quad,hist,wtime,AccuSpace,size,values,nhisto,diff,
                 min_val_x,max_val_x,min_val_y,max_val_y)
quadrant graph_quad; /* drawing area */
long hist[]; /* histogram */
InfMat3D *AccuSpace; /* 3D dynamic accumulator */
int wtime,size,values,nhisto,diff; /* histogram specifications */
float min_val_x,max_val_x,min_val_y,max_val_y; /* current accumulator ranges */
{
  int i,i2;

  if (nhisto)  
    for(i=size-1;(i>size-values-1) && (i>size-nhisto-1) && (i>-1);i--) {
      i2=i-diff;
      if (hist[i2] && i2>=0)
        draw_histoval_accu(graph_quad,100.*(i+1)/(float)wtime,AccuSpace,i+1,
                           min_val_x,max_val_x,min_val_y,max_val_y);
      else
        values++;
    }
  else 
    for(i=size-1;(i>size-values-1) && (i>-1);i--) {
      i2=i-diff;
      if (hist[i2] && i2>=0)
        draw_histoval_accu(graph_quad,100.*(i+1)/(float)wtime,AccuSpace,i+1,
                           min_val_x,max_val_x,min_val_y,max_val_y);
      else
        values++;
    }
}

/*****************************************************************************/
draw_histoval_accu(graph_quad,prob,AccuSpace,nhits,min_val_x,max_val_x,
                   min_val_y,max_val_y)
quadrant graph_quad; /* drawing area */
InfMat3D *AccuSpace; /* dynamic accumulator */
int nhits; /* accumulator score */
float prob; /* ? */
float min_val_x,max_val_x,min_val_y,max_val_y; /* current accumulator ranges */
{
  InfIndex3DX *p=AccuSpace;
  InfIndex3DY *s;
  InfIndex3DZ *r;
  float x,y,umin,umax,vmin,vmax,orig_x,orig_y,x_len,y_len,scale_x,scale_y;

  umin=0.0; /*(float)graph_quad.ulc_x;*/
  vmin=5.0; /*(float)graph_quad.ulc_y;*/
  umax=255.0; /*(float)graph_quad.lrc_x;*/
  vmax=200.0; /*(float)graph_quad.lrc_y;*/

  orig_x=umin+35.0;
  orig_y=vmax-15.0;

  x_len=umax-umin-orig_x-15.0;
  y_len=vmax-vmin-15.0-15.0;

  scale_x=x_len/(max_val_x-min_val_x);
  scale_y=y_len/(max_val_y-min_val_y);

  if (p->NextIndex==(InfIndex3DX *)NULL) return 0;
    p=p->NextIndex;
    do {
      s=p->YIndex;
      do {
        r=s->ZIndex;
        do {
          if (nhits == r->Accumulator) {
            x=(float)((r->Z-min_val_x)*scale_x+orig_x);
            y=(float)(-((prob-min_val_y)*scale_y)+orig_y);
            PutPixel_graph(graph_quad,nint(x),nint(y));
          }
          r=r->NextIndex;
        } while (r);
        s=s->NextIndex;
      } while (s);
      p=p->NextIndex;
    } while (p);
}

/*****************************************************************************/
heading(graph_quad,name)
quadrant graph_quad; /* drawing area */
char *name; /* label */
{
  float umin,umax,vmin,vmax,label_x,label_y;

  umin=0.0; /*(float)graph_quad.ulc_x;*/
  vmin=5.0; /*(float)graph_quad.ulc_y;*/
  umax=255.0; /*(float)graph_quad.lrc_x;*/
  vmax=200.0; /*(float)graph_quad.lrc_y;*/
  label_x=umin+100.0;
  label_y=vmin+12.0;

  PutText_graph(graph_quad,nint(label_x),nint(label_y),name);
}

/*****************************************************************************/
axis_x_abs_2(graph_quad,name,min_val,max_val)
quadrant graph_quad; /* drawing area */
char *name; /* label */
float min_val, max_val; /* current accumulator ranges */
{
  float orig_x,orig_y,umin,umax,vmin,vmax,end_x,label_x,label_y,
        label_inc,inc,stick_height,label_shift;
  char number[20];
  int i;

  umin=0.0; /*(float)graph_quad.ulc_x;*/
  vmin=5.0; /*(float)graph_quad.ulc_y;*/
  umax=255.0; /*(float)graph_quad.lrc_x;*/
  vmax=200.0; /*(float)graph_quad.lrc_y;*/

  orig_x=umin+35.0;
  orig_y=vmax-15.0;

  end_x=umax-umin-15.0;
  PutVector_graph(graph_quad, nint(orig_x), nint(orig_y),
	                       nint(end_x), nint(orig_y));
  label_x=umax-40.0;
  label_y=vmax;

  PutText_graph_small(graph_quad, nint(label_x), nint(label_y), name);

  label_x=orig_x;
  label_y=vmax-7.0;

  stick_height=2.0;

  inc=(max_val-min_val)/10.0;
  label_inc=(end_x-orig_x)/10.0;
  label_shift=label_inc*0.2;

  sprintf(number,"%-6.2f",min_val);

  PutText_graph_small(graph_quad, nint((label_x-label_shift)),
				  nint(label_y), number);

  PutVector_graph(graph_quad, nint(label_x), nint(orig_y-stick_height),
	                      nint(label_x), nint(orig_y+stick_height));
  for(i=0; i<10; i+=2) {
    label_x+=(2*label_inc);

    sprintf(number,"%-6.2f",min_val+=(2*inc));

    PutText_graph_small(graph_quad, nint(label_x-label_shift),
                                    nint(label_y), number);

    PutVector_graph(graph_quad,nint(label_x),nint(orig_y-stick_height),
                               nint(label_x),nint(orig_y+stick_height));
  }
}

/*****************************************************************************/
axis_y_abs_2(graph_quad,name,min_val,max_val)
quadrant graph_quad; /* drawing area */
char *name; /* label */
float min_val, max_val; /* current accumulator ranges */
{
  float orig_x,orig_y,umin,umax,vmin,vmax,end_y,label_x,label_y,
        label_inc,inc,stick_height,label_shift;
  char number[20];
  int i;

  umin=0.0; /*(float)graph_quad.ulc_x;*/
  vmin=5.0; /*(float)graph_quad.ulc_y;*/
  umax=255.0; /*(float)graph_quad.lrc_x;*/
  vmax=200.0; /*(float)graph_quad.lrc_y;*/

  orig_x=umin+35.0;
  orig_y=vmax-15.0;

  end_y=vmin+15.0;

  PutVector_graph(graph_quad, nint(orig_x), nint(orig_y),
                              nint(orig_x), nint(end_y));
  label_x=umin+1.0;
  label_y=vmin+8.0;

  PutText_graph_small(graph_quad, nint(label_x), nint(label_y), name);

  label_y=orig_y;
  label_x=umin+1.0;

  stick_height=2.0;

  inc=(max_val-min_val)/10.0;
  label_inc=(end_y-orig_y)/10.0;
  label_shift=5.0;

  sprintf(number,"%-6.2f",min_val);

  PutText_graph_small(graph_quad, nint(label_x), nint(label_y+2.0), number);
  PutVector_graph(graph_quad, nint(orig_x-stick_height), nint(label_y),
                              nint(orig_x+stick_height), nint(label_y));
  for(i=0; i<10; i+=2) {
    label_y+=(2*label_inc);

     sprintf(number,"%-6.2f",min_val+=(2*inc));

     PutText_graph_small(graph_quad, nint(label_x),nint(label_y+2.0),number);

     PutVector_graph(graph_quad,nint(orig_x-stick_height),nint(label_y),
	                        nint(orig_x+stick_height),nint(label_y));
  } 
}

/*****************************************************************************/
axis_x_abs_2_hist(graph_quad,name,size)
quadrant graph_quad; /* drawing area */
char *name; /* label */
int size; /* histogram size */
  {
  float orig_x,orig_y,umin,umax,vmin,vmax,end_x,label_x,label_y,label_inc;
  char number[20];
  int i;

  umin=0.0; /*(float)graph_quad.ulc_x;*/
  vmin=5.0; /*(float)graph_quad.ulc_y;*/
  umax=255.0; /*(float)graph_quad.lrc_x;*/
  vmax=200.0; /*(float)graph_quad.lrc_y;*/

  orig_x=umin+35.0;
  orig_y=vmax-15.0;

  end_x=umax-umin-15.0;

  PutVector_graph(graph_quad, nint(orig_x), nint(orig_y),
                              nint(end_x), nint(orig_y));
  label_x=umax-30.0-(float)(2.0*strlen(name));
  label_y=vmax;

  PutText_graph_small(graph_quad, nint(label_x), nint(label_y), name);

  label_y=vmax-7.0;
  label_inc=(end_x-orig_x)/(float)size;
  label_x=orig_x-label_inc/2.0-1.0;

  for(i=0; i<size; i++) {
    label_x+=(label_inc);
    if (i == 9)
      label_x-=2.0;
    sprintf(number,"%d",i+1);

    PutText_graph_small(graph_quad, nint(label_x),
					   nint(label_y), number);
  }
}

/*****************************************************************************/
DrawAccumulatorSpace_aht(quad,accu,dim_m,dim_c,tmp_pic,dim_x,dim_y,
                         small_m,big_m,small_c,big_c,
                         low_m,high_m,low_c,high_c)
quadrant quad; /* drawing area */
pic_type accu[][ACC_MAX_SIZE], tmp_pic[][MAX_SIZE]; /* accumulator */
int dim_m,dim_c,dim_x,dim_y; /*accumulator & image sizes*/
double small_m,big_m,small_c,big_c,low_m,high_m,low_c,high_c;     /* current */
{                                                      /* accumulator ranges */

  int i, j, acc_low_m=0, acc_high_m=dim_x,
            acc_low_c=0, acc_high_c=dim_y;
  pic_type max_pix_val=0;
  double acc_m, acc_c, pix_scale;

  ClearQuadrant(quad);

  init_pic(tmp_pic,dim_x,dim_y);

  for (i=0; i<dim_c; i++) /* finding global max value of the accumulator */
    for (j=0; j<dim_m; j++)
      if (accu[i][j]>max_pix_val)
        max_pix_val=accu[i][j];

  if (max_pix_val==(pic_type)0)
    return;

  pix_scale=255./(double)max_pix_val; /* scale gray values by the max value */

  if (small_m>low_m)
    acc_low_m=(int)((small_m-low_m)/(high_m-low_m)*dim_x);
  if (big_m>low_m && big_m<high_m)
    acc_high_m=(int)((((big_m-low_m)+((big_m-small_m)/(double)dim_m))/
                      (high_m-low_m))*dim_x);
  if (small_c>low_c)
    acc_low_c=(int)((small_c-low_c)/(high_c-low_c)*dim_y);
  if (big_c>low_c && big_c<high_c)
    acc_high_c=(int)((((big_c-low_c)+((big_c-small_c)/(double)dim_c))/
                      (high_c-low_c))*dim_y);

  if (acc_high_m>dim_x)
    acc_high_m=dim_x;
  if (acc_high_c>dim_y)
    acc_high_c=dim_y;

  acc_m=(double)dim_m/(double)(acc_high_m-acc_low_m);
  acc_c=(double)dim_c/(double)(acc_high_c-acc_low_c);

  for (i=acc_low_c; i<acc_high_c; i++) /* make tmp_pic to represent the accu */
    for (j=acc_low_m; j<acc_high_m; j++)
      tmp_pic[i][j]=((pic_type)(accu[(int)((i-acc_low_m)*acc_m)]
                                    [(int)((j-acc_low_c)*acc_c)]*pix_scale));

  PutAccu_gray(quad,tmp_pic,dim_x,dim_y); /* show the accu */
}
