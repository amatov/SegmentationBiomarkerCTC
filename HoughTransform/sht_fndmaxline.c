/*****************************************************************************
 * Houghtool - The Software Package for efficiency measuring and visualization
 * of the HT and it's variants for line detection
 *
 * Lappeenranta University of Technology, Department of Information Technology
 * Laboratory of Information Processing, Lappeenranta, Finland
 *
 * Authors: Petri Hirvonen, Jouni Ikonen (Jouni.Ikonen@lut.fi)
 *	    Pekka Kultanen, Heikki Kalviainen (Heikki.Kalviainen@lut.fi)
 *
 * File:    	sht_fndmaxline.c
 * Purpose: 	finding local maxima from the accumulator space & verifying lines
 *          	found
 * Date:    	Jun 1, 1993
 * Last change: Mar 8, 1996
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_graphmacros.h"

#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[];
#endif
extern int StopFlag;

extern pic_type hough_kernel[3][3];
extern coordinates edge_points[1000];

/*****************************************************************************/
FindMaxs_rho_theta(start,acc_space,n,threshold,dim_theta,dim_rho,edge,edge1,
                   dim_x,dim_y,maxs,min_seg_length,max_seg_width,
                   max_seg_gap,false_alarms)

/*
  Find n local maxima from Hough accumulator space and form corresponding
  linesegments and store results to global data structure line. Then form
  graph representation from linesegments and beautify formed graph. Graph
  is stored to given data structure line.

  Parameters : start - pointer to data structure line
               acc_space - accumulator space pic_type matrix,
                           size [][ACC_MAX_SIZE]
               n - number of maxima
               threshold - threshold for minimum value in acc_space
               dim_theta,dim_rho - size of the accumulator space
               edge - binary image containing edges pic_type matrix,
                      size [][MAX_SIZE]
               edge1 - binary image containing edges pic_type matrix
               dim_x,dim_y - size of the edge image
               maxs - arrays for storing maxima
               min_seg_length - minimum length for single linesegment
               max_seg_width -  maximum width for single linesegment
               max_seg_width -  maximum gap between pixels within
                                linesegment
	       false_alarms -   number of false alarms

*/
  
line **start;
pic_type acc_space[][ACC_MAX_SIZE],threshold,edge[][MAX_SIZE],edge1[][MAX_SIZE];
int n,dim_theta,dim_rho,dim_x,dim_y,min_seg_length,max_seg_width,
    *false_alarms;
double maxs[][2];
{
  int k,i,j,num_of_real_maxs=0,NumOfFalseMaxs=0;
  unsigned int val=MAXINT;
  double a,b,rho,theta,rho_min=-(double)dim_x,
         rho_max=sqrt((double)(dim_x*dim_x)+(double)(dim_y*dim_y));
/*
         rho_max=sqrt(2.0)*sqrt((double)dim_x*dim_x+(double)dim_y*dim_y);
*/
  for (k=0; k<n && val>min_seg_length && NumOfFalseMaxs<10 && !StopFlag; ) {

    find_max_value(acc_space,&val,&i,&j,dim_theta,dim_rho);

    if (val<threshold)
      return num_of_real_maxs;
/*
    theta=(double)((j-(double)dim_theta/2.0)*(PI/2.0)/((double)dim_theta/2.0));
*/
    theta=(double)j*PI/(double)dim_theta;
/*
    rho=(double)(i*2.0*rho_max/dim_rho-rho_max);
*/
    rho=(double)i*(rho_max-rho_min)/(double)(dim_rho-1)+rho_min;
    if (!fabs(sin(theta))) { /* (rho, theta) -> (a, b) */
      a=INF_DOUBLE_VAL;
      b=rint(rho);
    } else {
      a=-cos(theta)/sin(theta);
      if (isinf(a))
        b=rint(rho);
      else
        b=rho/sin(theta);
    }

    FollowAndMarkLineSegments_a_b(start,a,b,edge1,dim_x,dim_y,min_seg_length,
                                  max_seg_width,max_seg_gap);


    if (remove_edge_and_subtract_acc_space(*start,acc_space,dim_theta,dim_rho,
                                           edge,edge1,dim_x,dim_y)) {
      maxs[num_of_real_maxs][0]=rho;
      maxs[num_of_real_maxs][1]=theta;
      k++;
      num_of_real_maxs++;
      NumOfFalseMaxs=0;
      for (i=val=0; i<dim_y; i++)
        for (j=0; j<dim_x; j++)
          if (edge1[i][j]==(pic_type)OBJECT_PIX_VAL)
            val++;
    } else {
      acc_space[i][j]=0; /* removing only the local (global) maximum */
      (*false_alarms)++;
      NumOfFalseMaxs++;
    }

    PutAccu_gray_MACRO(quads[3],acc_space,dim_theta,dim_rho);

    NotifyDispatch_MACRO();
  }

  return num_of_real_maxs;
}

/*****************************************************************************/
remove_edge_and_subtract_acc_space(start,acc_space,dim_the,dim_rho,edge,
                                   edge1,dim_x,dim_y)
/*

  Remove edge points which lie under found line, from edge image
  and subtract corresponding Hough accumulator space points. 

  Parameters :  start     - beginning of the line data structure
                acc_space - Hough accumulator space, 
                            pic_type matrix size[][ACC_MAX_SIZE]
                dim_theta,dim_rho - size of the accumulator space
                edge,edge1 - binary images containing edge points,
                            pic_type matrix size[][ACC_MAX_SIZE]
                dim_x,dim_y - size of the edge images

  Remarks : No heuristic rules are used for controlling missubtraction
            from Hough accumulator space.

*/

line *start;
pic_type acc_space[][ACC_MAX_SIZE], edge[][MAX_SIZE], edge1[][MAX_SIZE];
int dim_the, dim_rho, dim_x, dim_y;
{
  line *p=start;
  coordinates *coord_pointer;
  double theta,rho,rho_min=-(double)dim_x,
         rho_max=sqrt((double)(dim_x*dim_x)+(double)(dim_y*dim_y));
/*
         rho_max=sqrt(2.0)*sqrt((double)dim_x*dim_x+(double)dim_y*dim_y);
*/

  int l,k,removed_points=0;

  while((p!=(line *)NULL) && !StopFlag) {
    if (! p->removed) {
      coord_pointer=p->edge_points;
      do {
        if (coord_pointer->x>=0 && coord_pointer->x<dim_x &&
            coord_pointer->y>=0 && coord_pointer->y<dim_y) {
          /* substract accumulator space */
          if (edge[coord_pointer->y][coord_pointer->x]=(pic_type)OBJECT_PIX_VAL)
            for (k=0; (k<dim_the) && !StopFlag; k++) {
/*
              theta=(double)(-PI/2.0+k*PI/dim_the);
*/
              theta=(double)k*PI/(double)dim_the;
              rho=(double)coord_pointer->x*cos(theta)+
                (double)coord_pointer->y*sin(theta);
/*
              l=(int)rint(rho*dim_rho/rho_max/2.0+dim_rho/2.0);
*/
              l=nint(((rho-rho_min)/(rho_max-rho_min))*(double)(dim_rho-1));
              sub_acc_space(hough_kernel,acc_space,dim_the,dim_rho,l,k);  
            }
          /* remove edges */
          edge1[coord_pointer->y][coord_pointer->x]=(pic_type)0;
          PutPixel_gray_MACRO(quads[0],coord_pointer->x,coord_pointer->y,0);
          removed_points++;
        }
      } while ( follow_segment(&coord_pointer) ); /* follow line segment */
    }
    p->removed=TRUE;
    p=p->next;
  }

  return removed_points;
}
