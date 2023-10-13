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
 * File:    ht_fndmaxline.c
 * Purpose: verifying lines found & storing discriptions of them
 * Date:    Jun 1, 1993
 *****************************************************************************/
#include "ht_hough.h"
#include "ht_graphmacros.h"
#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[];
#endif

coordinates edge_points[1000];

/*****************************************************************************/
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
		max_seg_width  - line scanning width
		max_seg_gap    - maximum gap between pixels within
                                 line segment

*/

line **start;
double a, b;
pic_type edge[][MAX_SIZE];
int dim_x, dim_y, min_seg_length, max_seg_width, max_seg_gap;
{
  int x,y,find_flag=FALSE,length,gap,start_x,start_y,end_x,end_y,num_of_points;

  if (fabs(a)<=1.0) /* low slope (x-dominant line) */
    if (fabs(a)<0.001) {
      y=nint(b);
      if (y>=0 && y<dim_y) /* scan the image along the line */
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
                if (num_of_points>=min_seg_length) {
                  set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                                    edge_points,num_of_points);
                  if (length>=min_seg_length)
                    add_line_segment(start,start_x,start_y,end_x,end_y,
                                     edge_points,num_of_points);
                }
                find_flag=FALSE;
              }
            }
          PutPixel_MACRO(quads[1],x,y);
          NotifyDispatch_MACRO();
        }
    } else
      for (x=0; x<dim_x; x++) { /* scan the image along the line */
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
                if (num_of_points>=min_seg_length) {
                  set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                                      edge_points,num_of_points);
                  if (length>=min_seg_length)
                    add_line_segment(start,start_x,start_y,end_x,end_y,
                                     edge_points,num_of_points);
                }
                find_flag=FALSE;
              }
            }
          PutPixel_MACRO(quads[1],x,y);
          NotifyDispatch_MACRO();
        }
      }
  else /* high slope (y-dominant line) */
    if (isinf(a)) {
      x=nint(b);
      if (x>=0 && x<dim_x) /* scan the image along the line */
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
                if (num_of_points>=min_seg_length) {
                  set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                                      edge_points,num_of_points);
                  if (length>=min_seg_length)
                    add_line_segment(start,start_x,start_y,end_x,end_y,
                                     edge_points,num_of_points);
                }
                find_flag=FALSE;
              }
            }
          PutPixel_MACRO(quads[1],x,y);
          NotifyDispatch_MACRO();
        }
    } else
      for (y=0; y<dim_y; y++) { /* scan the image along the line */
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
                if (num_of_points>=min_seg_length) {
                  set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                                      edge_points,num_of_points);
                  if (length>=min_seg_length)
                    add_line_segment(start,start_x,start_y,end_x,end_y,
                                     edge_points,num_of_points);
                }
                find_flag=FALSE;
              }
            }
          PutPixel_MACRO(quads[1],x,y);
          NotifyDispatch_MACRO();
        }
      }

  /* line segment leads to image border (??) */
  if (find_flag && num_of_points>=min_seg_length) {
    set_end_points_line(&start_x,&start_y,&end_x,&end_y,&length,
                        edge_points,num_of_points);
    if (length>=min_seg_length)
      add_line_segment(start,start_x,start_y,end_x,end_y,edge_points,
                       num_of_points);
  }
  FlushDisp_MACRO();
  NotifyDispatch_MACRO();
}

/*****************************************************************************/
some_match_y_dir(edge,dim_y,x,y,length)
/*

  Looks for points in binary image in y-direction in given location.
  Comparison is made length number of points up and down to detect
  points.
  Returns TRUE if found FALSE if not found. Value OBJECT_PIX_VAL corresponds to
  object and other values to background.

  Parameters : edge - binary image to match, pic_type matrix size [][MAX_SIZE]
               dim_y - image size in y-direction
               x, y - location where matching is done
               length - number of points to look for edge points

*/

pic_type edge[][MAX_SIZE];
int dim_y,x,y,length;
{
  register i;
  int low_y=y-length, high_y=y+length+1;

  if (low_y<0) low_y=0;
  if (high_y>dim_y) high_y=dim_y;

  for (i=low_y; i<high_y; i++)
    if (edge[i][x]==(pic_type)OBJECT_PIX_VAL)
      return TRUE;

  return FALSE;
}

/*****************************************************************************/
some_match_x_dir(edge,dim_x,x,y,length)

/*

  Looks for points in binary image in x-direction in given location.
  Comparison is made length number of points to left and to right to detect
  points.
  Returns TRUE if found FALSE if not found. Value OBJECT_PIX_VAL corresponds to
  object and other values to background.

  Parameters : edge - binary image to match, pic_type matrix size [][MAX_SIZE]
               dim_x - image size in x-direction
               x, y - location where matching is done
               length - number of points to look for edge points

*/

pic_type edge[][MAX_SIZE];
int dim_x,x,y,length;
{
  register j;
  int low_x=x-length, high_x=x+length+1;

  if (low_x<0) low_x=0;
  if (high_x>dim_x) high_x=dim_x;

  for (j=low_x; j<high_x; j++)
    if (edge[y][j]==(pic_type)OBJECT_PIX_VAL)
      return TRUE;

  return FALSE;
}

/*****************************************************************************/
gather_edge_points_y_dir(edge,dim_y,x,y,edge_points,num_of_points,length)

/*

  Looks for points in binary image in y-direction in given location.
  Comparison is made length number of points up and down to points.
  Found edge points are stored to line data structure in coordinates
  structure and number of edge points is returned via num_of_points.
  Value OBJECT_PIX_VAL corresponds to object and other values to background.

  Parameters : edge - binary image to match, pic_type matrix size [][MAX_SIZE]
               dim_y - image size in y-direction
               x, y - location where matching is done
               edge_points - pointer to data structure containing 
                             coordinates of found line segment
               num_of_points - number of found edge points (input/output)
               length - number of points to look for edge points

*/

pic_type edge[][MAX_SIZE];
int dim_y,x,y,*num_of_points,length;
coordinates *edge_points;
{
  register i;
  int low_y=y-length, high_y=y+length+1;

  if (low_y<0) low_y=0;
  if (high_y>dim_y) high_y=dim_y;

  for (i=low_y; i<high_y; i++) {
    PutPixel_gray_MACRO(quads[1],x,i,MID_GRAY_LEVEL);
    if (edge[i][x]==(pic_type)OBJECT_PIX_VAL) {
      edge_points[*num_of_points].x=x;
      edge_points[(*num_of_points)++].y=i;
      PutPixel_MACRO(quads[1],x,i);
    }
  }
}

/*****************************************************************************/
gather_edge_points_x_dir(edge,dim_x,x,y,edge_points,num_of_points,length)

/*

  Looks for points in binary image in x-direction in given location.
  Comparison is made length number of points to left and to right to points.
  Found edge points are stored to line data structure in coordinates
  structure and number of edge points is returned via num_of_points.
  Value OBJECT_PIX_VAL corresponds to object and other values to background.

  Parameters : edge - binary image to match, pic_type matrix size [][MAX_SIZE]
               dim_x - image size in x-direction
               x, y - location where matching is done
               edge_points - pointer to data structure containing 
                             coordinates of found line segment
               num_of_points - number of found edge points (input/output)
               length - number of points to look for edge points

*/

pic_type edge[][MAX_SIZE];
int dim_x,x,y,*num_of_points,length;
coordinates *edge_points;
{
  register j;
  int low_x=x-length, high_x=x+length+1;

  if (low_x<0) low_x=0;
  if (high_x>dim_x) high_x=dim_x;

  for (j=low_x; j<high_x; j++) {
    PutPixel_gray_MACRO(quads[1],j,y,MID_GRAY_LEVEL);
    if (edge[y][j]==(pic_type)OBJECT_PIX_VAL) {
      edge_points[*num_of_points].x=j;
      edge_points[(*num_of_points)++].y=y;
      PutPixel_MACRO(quads[1],j,y);
    }
  }
}

/*****************************************************************************/
set_end_points_line(start_x,start_y,end_x,end_y,length,edge_points,
                    num_of_points)
/*

  Set starting and ending points of the line segment. Also set the length of 
  the line segment.

  Parameters : start_x, start_y - starting coordinates of the line segment
               end_x, end_y - ending coordinates of the line segment
               length - length of the line segment
               edge_points -  pointer to data structure containing coordinates 
                              of found circle segment
               num_of_points - number of found edge points

*/

int *start_x, *start_y, *end_x, *end_y, *length, num_of_points;
coordinates *edge_points;
{
  (*start_x) = edge_points[0].x;
  (*start_y) = edge_points[0].y;
  (*end_x) = edge_points[num_of_points-1].x;
  (*end_y) = edge_points[num_of_points-1].y;
  (*length) = nint(sqrt( /* euclidian length of the line segment */
                 (double)((edge_points[0].x-edge_points[num_of_points-1].x)*
                          (edge_points[0].x-edge_points[num_of_points-1].x))+
                 (double)((edge_points[0].y-edge_points[num_of_points-1].y)*
                          (edge_points[0].y-edge_points[num_of_points-1].y))));
}

/*****************************************************************************/
add_line_segment(start,start_x,start_y,end_x,end_y,edge_points,num_of_edge_points)

/*

  Adds line segment to linked list of line data structures.

  Parameters : start_x,start_y - start coordinates of the line segment
               end_x,end_y - end coordinates of the line segment
               egde_points - line segment coordinates
               num_of_edge_points - number of edge points in line segment

*/

line **start;
int start_x, start_y, end_x, end_y;
coordinates *edge_points;
int num_of_edge_points;
{
  static count=0;
  int i;
  line *line_segment;

  if ((line_segment=(line *)malloc(sizeof(line)))==(line *)NULL) {
    err_prnt("allocating space for line segment");
    exit(-1);
  }

  add_segment(start,line_segment);

  line_segment->name = (++count);
  line_segment->start_label= -1;
  line_segment->start_x=start_x;
  line_segment->start_y=start_y;
  line_segment->end_label= -1;
  line_segment->end_x=end_x;
  line_segment->end_y=end_y;
  line_segment->removed=FALSE;

  if ((line_segment->edge_points=(coordinates *)malloc((num_of_edge_points+1)*
                                                         sizeof(coordinates)))
	                          == (coordinates *)NULL) {
    err_prnt("allocating space for coordinates of a line segment");
    exit(-1);
  }

  for (i=0; i<num_of_edge_points; i++) {
    line_segment->edge_points[i].x=edge_points[i].x;
    line_segment->edge_points[i].y=edge_points[i].y;
  }
  line_segment->edge_points[num_of_edge_points].x= -1;
  line_segment->edge_points[num_of_edge_points].y= -1;

  return TRUE;
}

/*****************************************************************************/
add_segment(start,line_segment)

/*
  Add line segment data structure to linked list to position
  start.

  Parameters : start - pointer to pointer to line segment data structure
					   (input/output)
               line_segment - pointer to linesegment data structure to
			                  be added to linked list

*/

line **start, *line_segment;
{
  line_segment->next= *start;
  *start=line_segment;
}

/*****************************************************************************/
follow_segment(coord_pointer)

/*
  Follow found linesegment and increment pointer in linesegments
  coordinate vector.

  Parameters : coord_pointer - pointer to linesegment coord_pointer
                               data structure (input/output)

*/

coordinates **coord_pointer;
{
  (*coord_pointer)++;
  if ((**coord_pointer).x < 0)
    return FALSE;
  else
    return TRUE;
}


/*****************************************************************************/
DeleteAllLineSegments(start)

/*
  Delete line segment structure and free memory. 

  Parameters : start - pointer to linked list which contains
                       line segment data structures

*/

line *start;
{
  line *p=start ,*tmp_p;

  tmp_p=p;
  if (p!=(line *)NULL) {
    p->name=p->start_x=p->start_y=p->start_label=p->end_x=p->end_y=p->end_label=0;
    p->removed=1;
    free(p->edge_points);
    p->edge_points=(coordinates *)NULL;

    p=p->next;
    tmp_p->next=(line *)NULL;

    while (p!=(line *)NULL) {
      tmp_p=p->next;
      free(p->edge_points);
      free(p);
      p=tmp_p;
    }
  }
}

/*****************************************************************************/
find_max_value(acc_space,val,x,y,dim_the,dim_rho)

/*

  Finds maximum value from Hough accumulator space.

  Parameters : acc_space - Hough accumulator space,
                           pic_type matrix, size [][ACC_MAX_SIZE]
               val - found maximum value (output)
               x,y - location in acc_space matrix which contains maxima
               dim_the, dim_rho - size of the acc_space

*/
pic_type acc_space[][ACC_MAX_SIZE], *val;
int *x,*y,dim_the,dim_rho;
{
  int i,j;

  *val=0;

  for (i=0; i<dim_rho; i++)
    for (j=0; j<dim_the; j++)
      if (acc_space[i][j]>(*val)) {
        *val=acc_space[i][j];
        *x=i;
        *y=j;
      }
}

/*****************************************************************************/
RemoveEdgePoints(start,pic,dim_x,dim_y)
/*

  Remove edge points from the image corresponding to those of stored in the
  line structure.

  Parameters : start - pointer to linked list which contains
                       line segment data structures
               pic - edge image, pic_type matrix, size [][ACC_MAX_SIZE]
               dim_x, dim_y - size of the edge image

*/
line *start;
pic_type pic[][MAX_SIZE];
int dim_x,dim_y;
{
  line *p=start;
  coordinates *coord_pointer;
  int removed_points=0;

  while(p!=(line *)NULL) {
    if (! p->removed) {
      coord_pointer=p->edge_points;
      do {
        if (coord_pointer->x>=0 && coord_pointer->x<dim_x &&
            coord_pointer->y>=0 && coord_pointer->y<dim_y) {

          /* remove point from the edge image */
          pic[coord_pointer->y][coord_pointer->x]=(pic_type)0;
          PutPixel_gray_MACRO(quads[0],coord_pointer->x,coord_pointer->y,0);
          removed_points++;
        }
      } while ( follow_segment(&coord_pointer) ); /* follow line segment */

      NotifyDispatch_MACRO();
      FlushDisp_MACRO();
    }
    p->removed=TRUE;
    p=p->next;
  }

  return removed_points;
}
