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
 * File:    cht_fndmaxline.c
 * Purpose: finding local maxima from the accumulator space
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_graphmacros.h"

#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[];
#endif
extern int StopFlag;

extern coordinates edge_points[1000];
int peaks[MAX_LINES][3];

/*****************************************************************************/

FindMaxs_rho_theta_one(start,acc_space,n,threshold,dim_theta,dim_rho,edge,
                   dim_x,dim_y,maxs,min_seg_length,max_seg_width,
                   max_seg_gap,false_alarms,mask_size)

/*
  Find n local maximas from Hough accumulator space and form corresponding
  linesegments and store results to global data structure line. Then form
  graph representation from linesegments and beautify formed graph. Graph
  is stored to given data structure line.

  Parameters : start - pointer to data structure line
               acc_space - accumulator space pic_type matrix,
                           size [][ACC_MAX_SIZE]
               n - number of maximas to be found
               threshold - threshold for minimum value in acc_space
               dim_theta,dim_rho - size of the accumulator space
               edge - binary image containing edges pic_type matrix,
                      size [][MAX_SIZE]
               dim_x,dim_y - size of the edge image
               maxs - arrays for storing maximas
               min_seg_length - minimum length for single linesegment
               max_seg_width -  maximum width for single linesegment
               max_seg_width -  maximum gap between pixels within
                                linesegment
	       false_alarms -   number of false alarms
	       mask_size -   size of the local maxima search window

  Remarks :
            
*/
  
line **start;
pic_type acc_space[][ACC_MAX_SIZE],threshold,edge[][MAX_SIZE];
int n,dim_theta,dim_rho,dim_x,dim_y,min_seg_length,max_seg_width,
    *false_alarms,mask_size;
double maxs[][2];
{
  int i,j,k,num_of_real_maxs=0,NumOfFalseMaxs=0,num_of_peaks=0,min_val=0;
  unsigned int val=MAXINT;
  double a,b,rho,theta,rho_min=-(double)dim_x,
         rho_max=sqrt((double)(dim_x*dim_x)+(double)(dim_y*dim_y));

  for (i=0; i<dim_rho && !StopFlag; i++) /* sort accumulator peaks */
    for (j=0; j<dim_theta && !StopFlag; j++) {
      if (acc_space[i][j]>=threshold && acc_space[i][j]>min_val &&
          detect_accu_peak(acc_space,dim_rho,dim_theta,i,j,mask_size))
        update_peak_list(peaks,n,i,j,acc_space[i][j],&num_of_peaks,&min_val);
      NotifyDispatch_MACRO();
    }

  StopFlag=0;

  for (k=0; k<num_of_peaks && val>(min_seg_length*2) && NumOfFalseMaxs<10 &&
            !StopFlag; k++) {

    i=peaks[k][0];
    j=peaks[k][1];
    acc_space[i][j]=0;

    theta=(double)j*PI/(double)dim_theta; /* calc parameters */
    rho=(double)i*(rho_max-rho_min)/(double)(dim_rho-1)+rho_min;
    if (!fabs(sin(theta))) { /* (rho, theta)  --> (a, b) */
      a=INF_DOUBLE_VAL;
      b=rint(rho);
    } else {
      a=-cos(theta)/sin(theta);
      if (isinf(a))
        b=rint(rho);
      else
        b=rho/sin(theta);
    }
 
    /* find endpoits and store lines */
    FollowAndMarkLineSegments_a_b(start,a,b,edge,dim_x,dim_y,min_seg_length,
                                  max_seg_width,max_seg_gap);


    if (RemoveEdgePoints(*start,edge,dim_x,dim_y)) {
      maxs[num_of_real_maxs][0]=rho;
      maxs[num_of_real_maxs][1]=theta;
      num_of_real_maxs++;
      NumOfFalseMaxs=0;
      for (i=val=0; i<dim_y; i++)
        for (j=0; j<dim_x; j++)
          if (edge[i][j]==(pic_type)255)
            val++;
    } else {
      (*false_alarms)++;
      NumOfFalseMaxs++;
    }

    PutAccu_gray_MACRO(quads[3],acc_space,dim_theta,dim_rho);

    NotifyDispatch_MACRO();
  }

  return num_of_real_maxs;
}

/*****************************************************************************/
detect_accu_peak(acc_space,dim_rho,dim_theta,i,j,s)
/*

  Check if there is a peak in the masked area of the accumulator array.

  Parameters : acc_space - accumulator space pic_type matrix,
                           size [][ACC_MAX_SIZE]
               dim_theta,dim_rho - size of the accumulator space
               i,j - accumulator space coordinates
               s - peak detection mask size

*/
pic_type acc_space[][ACC_MAX_SIZE];
int dim_rho,dim_theta,i,j,s;
{
  int min_i=(i-s)<0?0:i-s, max_i=(i+1+s)>dim_rho?dim_rho:i+1+s,
      min_j=(j-s)<0?0:j-s, max_j=(j+1+s)>dim_theta?dim_theta:j+1+s,
      x,y,i_peak,j_peak,val=0;

  for (x=min_i; x<max_i; x++) /* peak detection mask */
    for (y=min_j; y<max_j; y++)
      if (acc_space[x][y]>val) {
        val=acc_space[x][y];
        i_peak=x;
        j_peak=y;
      }

  if (i==i_peak && j==j_peak)
    return TRUE; /* there is a peak, indeed */
  return FALSE;
}

/*****************************************************************************/
update_peak_list(peaks,n,i,j,val,num_of_peaks,min_val)
/*

  Add a new peak to the list and sort the list.

  Parameters :	peaks - accumulator peak list
		n - peak list size
		i,j - accumulator space coordinates
		val - peak score
		dim_x, dim_y - edge image size
		num_of_peaks - number of 'old' peaks
		min_val - lowest peak score

*/
int peaks[][3],n,i,j,val,*num_of_peaks,*min_val;
{
  int x=0,y;

  if ((*num_of_peaks)==0) {
    peaks[0][0]=i;
    peaks[0][1]=j;
    peaks[0][2]=val;
    (*num_of_peaks)=1;
    return;
  }

  while (x<(*num_of_peaks) && peaks[x][2]>=val) x++;

  if (x>=n) /* peak was not high enough? */
    return;

  for (y=((*num_of_peaks)+1)>n?n:(*num_of_peaks)+1; y>x-1; y--) {
    peaks[y][0]=peaks[y-1][0];
    peaks[y][1]=peaks[y-1][1];
    peaks[y][2]=peaks[y-1][2];
  }
  peaks[x][0]=i;
  peaks[x][1]=j;
  peaks[x][2]=val;
  if ((*num_of_peaks)<n)
    (*num_of_peaks)++;
  if ((*num_of_peaks)==n)
    (*min_val)=peaks[n-1][2];
}
