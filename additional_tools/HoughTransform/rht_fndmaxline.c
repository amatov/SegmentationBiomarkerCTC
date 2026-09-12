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
 * File: rht_fndmaxline.c: (functions for finding local maxima from the
 *                          dynamic accumulator space)
 *****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"
#include "ht_graphmacros.h"
#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[];
#endif

coordinates edge_points[1000];

FollowAndMarkLineSegments_a_b(start,a,b,edge,dim_x,dim_y,min_seg_length,
                              max_seg_width,max_seg_gap)
/*

  Follow line defined by Hough accumulator space maxima location and
  form line segments according to given minimum segment length and 
  given deviation from absolute line.

  Parameters :	start - data structure for line segments
		a, b - Hough accumulator space coordinates
		edge - binary image containing  edge points, 
		       pic_type matrix size [][MAX_SIZE]
		dim_x, dim_y - edge image size
		min_seg_length - minimum legth for line segment
		max_seg_width  - maximum width for line segment
		max_seg_gap    - maximum gap between pixels within
                                 line segment

  Remarks : ...

*/

line **start;
double a, b;
pic_type edge[][MAX_SIZE];
int dim_x, dim_y, min_seg_length, max_seg_width, max_seg_gap;
{
  int x,y,find_flag=FALSE,length,gap,start_x,start_y,end_x,end_y,num_of_points;

  if (fabs(a)<=1.0)
    if (fabs(a)<0.001) {
      y=nint(b);
      if (y>=0 && y<dim_y)
        for (x=0; x<dim_x; x++) {
          if ( some_match_y_dir(edge,dim_y,x,y,max_seg_width) ) {
            if (find_flag)  
              gather_edge_points_y_dir(edge,dim_y,x,y,edge_points,
                                       &num_of_points,max_seg_width);
            else {
              find_flag=TRUE;
              num_of_points=0;
              gather_edge_points_y_dir(edge,dim_y,x,y,edge_points,
                                       &num_of_points,max_seg_width);
            }
            gap=0;
          } else
            if (find_flag) {
              if (gap<max_seg_gap)
                gap++;
              else {
                set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                                    edge_points,num_of_points);
                if (length>=min_seg_length)
                  add_line_segment(start,start_x,start_y,end_x,end_y,
                                   edge_points,num_of_points);
                find_flag=FALSE;
              }
            }
          PutPixel_MACRO(quads[1],x,y); /* PutPixel() or empty */
          NotifyDispatch_MACRO();
        }
    } else
      for (x=0; x<dim_x; x++) {
        y=nint(a*x + b);
        if (y>=0 && y<dim_y) {
          if ( some_match_y_dir(edge,dim_y,x,y,max_seg_width) ) {
            if (find_flag) 
              gather_edge_points_y_dir(edge,dim_y,x,y,edge_points,
                                       &num_of_points,max_seg_width);
            else {
              find_flag=TRUE;
              num_of_points=0;
              gather_edge_points_y_dir(edge,dim_y,x,y,edge_points,
                                       &num_of_points,max_seg_width);
            }
            gap=0;
          } else
            if (find_flag) {
              if (gap<max_seg_gap)
                gap++;
              else {
                set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                                    edge_points,num_of_points);
                if (length>=min_seg_length)
                  add_line_segment(start,start_x,start_y,end_x,end_y,
                                   edge_points,num_of_points);
                find_flag=FALSE;
              }
            }
          PutPixel_MACRO(quads[1],x,y);
          NotifyDispatch_MACRO();
        }
      }
  else
    if (isinf(a)) {
      x=nint(b);
      if (x>=0 && x<dim_x)
        for (y=0; y<dim_y; y++) {
          if ( some_match_x_dir(edge,dim_x,x,y,max_seg_width) ) {
            if (find_flag) 
              gather_edge_points_x_dir(edge,dim_x,x,y,edge_points,
                                       &num_of_points,max_seg_width);
            else {
              find_flag=TRUE;
              num_of_points=0;
              gather_edge_points_x_dir(edge,dim_x,x,y,edge_points,
                                       &num_of_points,max_seg_width);
            }
            gap=0;
          } else
            if (find_flag) {
              if (gap<max_seg_gap)
                gap++;
              else {
                set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                                    edge_points,num_of_points);
                if (length>=min_seg_length)
                  add_line_segment(start,start_x,start_y,end_x,end_y,
                                   edge_points,num_of_points);
                find_flag=FALSE;
              }
            }
          PutPixel_MACRO(quads[1],x,y);
          NotifyDispatch_MACRO();
        }
    } else
      for (y=0; y<dim_y; y++) {
        x=nint((y - b)/a);
        if (x>=0 && x<dim_x) {
          if ( some_match_x_dir(edge,dim_x,x,y,max_seg_width) ) {
            if (find_flag) 
              gather_edge_points_x_dir(edge,dim_x,x,y,edge_points,
                                       &num_of_points,max_seg_width);
            else {
              find_flag=TRUE;
              num_of_points=0; 
              gather_edge_points_x_dir(edge,dim_x,x,y,edge_points,
                                       &num_of_points,max_seg_width);
            }
            gap=0;
          } else
            if (find_flag) {
              if (gap<max_seg_gap)
                gap++;
              else {
                set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                                    edge_points,num_of_points);
                if (length>=min_seg_length)
                  add_line_segment(start,start_x,start_y,end_x,end_y,
                                   edge_points,num_of_points);
                find_flag=FALSE;
              }
            }
          PutPixel_MACRO(quads[1],x,y);
          NotifyDispatch_MACRO();
        }
      }
  FlushDisp_MACRO();
  NotifyDispatch_MACRO();
}

RemoveEdgePoints(start,pic,dim_x,dim_y,vec,n)
line *start;
pic_type pic[][MAX_SIZE];
int dim_x,dim_y;
pic_vec_type vec[];
long *n;
{
  line *p=start;
  coordinates *coord_pointer;
  long n1=*n;

  while(p!=(line *)NULL) {
    if (! p->removed) {
      coord_pointer=p->edge_points;
      if (coord_pointer->x>=0 && coord_pointer->x<dim_x &&
          coord_pointer->y>=0 && coord_pointer->y<dim_y)
        do {
          /* poista ko. piste reunakuvasta */
          pic[coord_pointer->y][coord_pointer->x]=(pic_type)0;
          PutPixel_gray_MACRO(quads[0],coord_pointer->x,coord_pointer->y,0);
        } while ( follow_segment(&coord_pointer) ); /* follow line segment */

      NotifyDispatch_MACRO();
      FlushDisp_MACRO();
    }
    p->removed=TRUE;
    p=p->next;
  }

  convert_pic_vec(pic,dim_x,dim_y,vec,n);

  return (int)((*n)-n1);
}
