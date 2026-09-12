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
 * File:    drht_blockpixels.c
 * Purpose: blocking image pixels
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"

#define MAX_D 20

coordinates edge_points[5000];

/*****************************************************************************/
FollowAndGetBlockedPixels_Line(a,b,edge,dim_x,dim_y,d,max_var,vec1,count)
/*
 
  Follow the line defined by Hough accumulator space maxima location and
  get all blocked "white" pixels (edge segments).

  Parameters : a, b - Hough accumulator space coordinates 
               edge - binary image containing  edge points, 
                      pic_type matrix size [][MAX_SIZE]
               dim_x, dim_y - edge image size
               d - block width
               max_var - maximum variation of slope (a) in degrees
               vec1 - table in which pixel coordinates are stored
               count - number of pixels stored

*/
double a,b,max_var;
pic_type edge[][MAX_SIZE];
int dim_x,dim_y,d;
pic_vec_type vec1[];
long *count;
{
  int x,y,num_of_points=0,more_d=0,at_least_d;
  double delta_angle,angle,x1=0.0,y1=0.0;

  at_least_d=d;
  delta_angle=((PI/180.0)*max_var)/2.0;

  if (fabs(a)<=1.0) { /* low slope (x-dominant line) */
    angle=a>=0.0?atan(a):-(atan(a));
    if (fabs(a)<0.001) {
      y=nint(b);
      if (y>=0 && y<dim_y)
        for (x=0; x<dim_x; x++) {
          d=at_least_d; /* setting block width */
          x1=(double)x;
          more_d=abs(nint(ceil((x1*(double)(tan(angle)-
                     tan(angle+delta_angle))))));
          if (more_d>d)
            d=more_d;
          if (more_d>MAX_D)
            d=MAX_D;
          if ( some_match_y_dir(edge,dim_y,x,y,d) ) 
            gather_edge_points_y_dir(edge,dim_y,x,y,edge_points,&num_of_points,d);
        } 
    } else  
      for (x=0; x<dim_x; x++) {
        d=at_least_d; /* setting block width */
        y=nint(a*x + b);
        x1=(double)x;
        more_d=abs(nint(ceil((x1*(double)(tan(angle)-
                   tan(angle+delta_angle))))));
       if (more_d>d)
         d=more_d;
       if (more_d>MAX_D)
         d=MAX_D;
       if (y>=0 && y<dim_y)           
         if ( some_match_y_dir(edge,dim_y,x,y,d) ) 
           gather_edge_points_y_dir(edge,dim_y,x,y,edge_points,&num_of_points,d);
     }   
  } else {/* high slope (y-dominant line) */
    angle=a>=0.0?(PI/2.0)-atan(a):(PI/2.0)+atan(a);
    if (isinf(a)) {
      x=nint(b);
      if (x>=0 && x<dim_x)
        for (y=0; y<dim_y; y++) {
          d=at_least_d; /* setting block width */
          y1=(double)y;
          more_d=abs(nint(ceil((y1*(double)(tan(angle)-
                     tan(angle+delta_angle))))));
          if (more_d>d)
            d=more_d;
          if (more_d>MAX_D)
            d=MAX_D;
          if ( some_match_x_dir(edge,dim_x,x,y,d) )     
            gather_edge_points_x_dir(edge,dim_x,x,y,edge_points,&num_of_points,d);
        }
    } else
      for (y=0; y<dim_y; y++) {
        d=at_least_d; /* setting block width */
        x=nint((y - b)/a);
        y1=(double)y;         
        more_d=abs(nint(ceil((y1*(double)(tan(angle)-
                   tan(angle+delta_angle))))));
        if (more_d>d)
          d=more_d;
        if (more_d>MAX_D)
          d=MAX_D;
        if (x>=0 && x<dim_x)
          if ( some_match_x_dir(edge,dim_x,x,y,d) )     
            gather_edge_points_x_dir(edge,dim_x,x,y,edge_points,&num_of_points,d);
      }
  }

  for ( *count=0; *count<num_of_points; (*count)++) {
    vec1[*count].i=edge_points[*count].y;
    vec1[*count].j=edge_points[*count].x;
  }
}
