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
 * File:    graphics.c
 * Purpose: drawing graphics, used with XHoughtool only
 * Date:    Jun 1, 1993
 *****************************************************************************/
 
#include "ht_hough.h"
#include "xhoughtool.h"

extern GC gc;
extern Xv_Font font, small_font;
extern XFontStruct *cur_font;
extern Display *dpy;
extern XID xid_canvas, xid_graph_canvas;
extern unsigned long *pixel_table, bg_level, fg_level;
extern inversed_graylevels;

/*****************************************************************************/
FlushDisp()
/*** flushes display buffer onto the display ***/
{
  XFlush(dpy);
}

/*****************************************************************************/
change_font(fontid)
/*** changes font between normal and small ***/
int fontid;
{
  if (fontid==XHT_SMALL_FONT)
    cur_font = (XFontStruct *)xv_get(small_font, FONT_INFO);
  else
    cur_font = (XFontStruct *)xv_get(font, FONT_INFO);
  XSetFont(dpy, gc, cur_font->fid);
}

/*****************************************************************************/
PutPic(quad,pic,dim_x,dim_y)
/*** displays binary image on a window area ***/
quadrant quad; /* area */
pic_type pic[][MAX_SIZE]; /* image */
int dim_x, dim_y; /* image size */
{
  int i, j, prev_gray_index=127, new_gray_index;

  ClearQuadrant(quad);

  for (i=0; (i<dim_y && i<QUAD_Y_SIZE); i++)
    for (j=0; (j<dim_x && j<QUAD_X_SIZE); j++)
      if (pic[i][j])
        XDrawPoint(dpy, xid_canvas, gc, quad.ulc_x+j, quad.ulc_y+i);

  XFlush(dpy);
}

/*****************************************************************************/
PutPic_gray(quad,pic,dim_x,dim_y)
/*** displays gray level image on a window area ***/
quadrant quad; /* area */
pic_type pic[][MAX_SIZE]; /* image */
int dim_x, dim_y; /* image size */
{
  int i, j, prev_gray_index=127, new_gray_index;

  if (inversed_graylevels) /* white background */
    prev_gray_index=0;

  ClearQuadrant(quad);

  if (inversed_graylevels) {
    prev_gray_index=0;
    for (i=0; (i<dim_y && i<QUAD_Y_SIZE); i++)
      for (j=0; (j<dim_x && j<QUAD_X_SIZE); j++)
        if (pic[i][j]) {
          if ((new_gray_index=127-(pic[i][j]/2)) != prev_gray_index) {
            /* set new drawing graylevel if necessary */
            XSetForeground(dpy,gc,pixel_table[new_gray_index]);
            prev_gray_index=new_gray_index;
          }
          XDrawPoint(dpy, xid_canvas, gc, quad.ulc_x+j, quad.ulc_y+i);
        }
  } else
    for (i=0; (i<dim_y && i<QUAD_Y_SIZE); i++)
      for (j=0; (j<dim_x && j<QUAD_X_SIZE); j++)
        if (pic[i][j]) {
          if ((new_gray_index=(pic[i][j]/2)) != prev_gray_index) {
            /* set new drawing graylevel if necessary */
            XSetForeground(dpy,gc,pixel_table[new_gray_index]);
            prev_gray_index=new_gray_index;
          }
          XDrawPoint(dpy, xid_canvas, gc, quad.ulc_x+j, quad.ulc_y+i);
        }

  /* reset default foreground graylevel */
  XSetForeground(dpy,gc,fg_level);

  XFlush(dpy);
}

/*****************************************************************************/
PutAccu_gray(quad,accu,dim_x,dim_y)
/*** displays static accumulator on a window area ***/
quadrant quad; /* area */
pic_type accu[][MAX_SIZE]; /* accumulator array */
int dim_x, dim_y; /* array size */
{
  int i, j, prev_gray_index=127, new_gray_index;
  pic_type max_pix_val=0;
  double pix_scale;

  ClearQuadrant(quad);

  for (i=0; i<dim_y; i++) /* finding global max of the accumulator */
    for (j=0; j<dim_x; j++)
      if (accu[i][j]>max_pix_val)
        max_pix_val=accu[i][j];

  if (max_pix_val==(pic_type)0)
    return;

  pix_scale=255.0/(double)max_pix_val; /* scale gray values by the global maxima */

  for (i=0; (i<dim_y && i<QUAD_Y_SIZE); i++)
    for (j=0; (j<dim_x && j<QUAD_X_SIZE); j++)
      if (accu[i][j]) {
        if ((new_gray_index=(nint(accu[i][j]*pix_scale)/2)) != prev_gray_index) {
          /* set new drawing graylevel if necessary */
          XSetForeground(dpy,gc,pixel_table[new_gray_index]);
          prev_gray_index=new_gray_index;
        }
        XDrawPoint(dpy, xid_canvas, gc, quad.ulc_x+j, quad.ulc_y+i);
      }

  /* reset default graylevel (white) */
  XSetForeground(dpy,gc,fg_level);

  XFlush(dpy);
}

/*****************************************************************************/
PutPixel(quad,x,y)
/*** displays one point on a window area ***/
quadrant quad; /* area */
int x,y; /* point coordinates */
{
  if (x>=0 && x<QUAD_X_SIZE && y>=0 && y<QUAD_X_SIZE)
    XDrawPoint(dpy, xid_canvas, gc, quad.ulc_x+x, quad.ulc_y+y);
}

/*****************************************************************************/
PutPixel_gray(quad,x,y,val)
/*** displays a gray level point ***/
quadrant quad; /* area */
int x,y,val; /* point coordinates, gray level */
{
  if (val==MID_GRAY_LEVEL) /* not white and not black */
    val=inversed_graylevels?(pic_type)70:(pic_type)150;
  if (val!=(pic_type)OBJECT_PIX_VAL) {
    XSetForeground(dpy,gc,pixel_table[val/2]);
    if (x>=0 && x<QUAD_X_SIZE && y>=0 && y<QUAD_X_SIZE)
      XDrawPoint(dpy, xid_canvas, gc, quad.ulc_x+x, quad.ulc_y+y);
    XSetForeground(dpy,gc,fg_level);
  } else
    if (x>=0 && x<QUAD_X_SIZE && y>=0 && y<QUAD_X_SIZE)
      XDrawPoint(dpy, xid_canvas, gc, quad.ulc_x+x, quad.ulc_y+y);
}

/*****************************************************************************/
PutPixel_graph(graph_quad,x,y)
/*** displays a point on a graph area ***/
quadrant graph_quad; /* graph area */
int x,y; /* point coordinates */
{
  XDrawPoint(dpy, xid_graph_canvas, gc, graph_quad.ulc_x+x, graph_quad.ulc_y+y);
}

/*****************************************************************************/
PutVector(quad,x1,y1,x2,y2)
/*** displays a vector on a window area ***/
quadrant quad; /* area */
int x1,y1,x2,y2; /* end point coordinates */
{
  XDrawLine(dpy, xid_canvas, gc,
            quad.ulc_x+x1, quad.ulc_y+y1, quad.ulc_x+x2, quad.ulc_y+y2);
  XFlush(dpy);
}

/*****************************************************************************/
PutVector_graph(graph_quad,x1,y1,x2,y2)
/*** displays a vector on a graph area ***/
quadrant graph_quad; /* graph area */
int x1,y1,x2,y2; /* end point coordinates */
{
  XDrawLine(dpy, xid_graph_canvas, gc,
            graph_quad.ulc_x+x1, graph_quad.ulc_y+y1,
            graph_quad.ulc_x+x2, graph_quad.ulc_y+y2);
  XFlush(dpy);
}

/*****************************************************************************/
PutText(quad,x,y,str)
/*** displays a text string on a window area ***/
quadrant quad; /* area */
int x,y; /* starting coordinates */
char *str; /* text string */
{
  XDrawString(dpy, xid_canvas, gc, quad.ulc_x+x, quad.ulc_y+y, str, strlen(str));
  XFlush(dpy);
}

/*****************************************************************************/
PutLongText(quad,x,y,val,draw)
/*** displays a long value as a text string on a window area ***/
quadrant quad; /* area */
int x,y,draw; /* starting coordinates, condition */
long val; /* value to be displayed */
{
  char str[30];
  static int prev_x=0, prev_y=0;
  static long prev_val=0L;

  if (draw) { /* condition */
    if (prev_x==x && prev_y==y) { /* `unprint' the old text */
      XSetForeground(dpy,gc,bg_level);
      sprintf(str,"%ld",prev_val);
      PutText(quad,x,y,str);
      XSetForeground(dpy,gc,fg_level);
    }
    sprintf(str,"%ld",val);
    PutText(quad,x,y,str);
    prev_x=x; prev_y=y;
    prev_val=val;
  }
}

/*****************************************************************************/
PutText_graph(graph_quad,x,y,str)
/*** displays a text string on a graph area ***/
quadrant graph_quad; /* graph area */
int x,y; /* starting coordinates */
char *str; /* text string */
{
  XDrawString(dpy, xid_graph_canvas, gc,
              graph_quad.ulc_x+x, graph_quad.ulc_y+y, str, strlen(str));
  XFlush(dpy);
}

/*****************************************************************************/
PutText_graph_small(graph_quad,x,y,str)
/*** displays a text string on a graph area using small font ***/
quadrant graph_quad; /* graph area */
int x,y; /* starting coordinates */
char *str; /* text string */
{
  change_font(XHT_SMALL_FONT);
  XDrawString(dpy, xid_graph_canvas, gc,
              graph_quad.ulc_x+x, graph_quad.ulc_y+y, str, strlen(str));
  XFlush(dpy);
  change_font(XHT_FONT);
}

/*****************************************************************************/
ClearQuadrant(quad)
/*** glears a window area ***/
quadrant quad; /* area */
{
  XClearArea(dpy,xid_canvas,quad.ulc_x,quad.ulc_y,QUAD_X_SIZE,QUAD_Y_SIZE,0);
  XFlush(dpy);
}

/*****************************************************************************/
ClearQuadrant_graph(graph_quad)
/*** glears a graph area ***/
quadrant graph_quad; /* graph area */
{
  XClearArea(dpy,xid_graph_canvas,graph_quad.ulc_x,graph_quad.ulc_y,256,200,0);
  XFlush(dpy);
}

/*****************************************************************************/
ShowImage(quad,filename,shrink)
/*** loads and shows an image on the screen (128 graylevels used only) ***/
quadrant quad; /* area */
char *filename; /* filename of the image to be displayed */
int shrink; /* shrink on/off */
{
  int i, j, dim_x, dim_y, prev_gray_index=127, new_gray_index;

  pic_type pic[MAX_SIZE][MAX_SIZE];

  if (!pic_from_disk(filename,pic,shrink,&dim_x,&dim_y))
    return 0;

  PutPic_gray(quad,pic,dim_x,dim_y);

  return 1;
}

/*****************************************************************************/
DrawLines(quad,start)
/*** display lines on a window area ***/
quadrant quad; /* area */
line *start; /* line structure */
{
  line *p=start;

  while (p!=(line *)NULL) {
    PutVector(quad,p->start_x,p->start_y,p->end_x,p->end_y);
    p=p->next;
  }

  XFlush(dpy);
}

/*****************************************************************************/
DrawCircles(quad,circle_start)
/*** display circles on a window area ***/
quadrant quad; /* area */
circle *circle_start; /* circle structure */
{
  circle *p = circle_start;

  while (p!=(circle *)NULL) {
    draw_circle(quad,p->start_x,p->start_y,p->end_x,p->end_y,
                p->centre_x,p->centre_y,p->radius);
    p=p->next;
  }

  XFlush(dpy);
}

/*****************************************************************************/
draw_circle(quad,s_x,s_y,e_x,e_y,a,b,radius)
/*** display a circle on a window area ***/
quadrant quad; /* area */
int s_x,s_y,e_x,e_y; /* end point coordinates */
double a,b,radius; /* circle center points and radius */
{
  double t,s_t,e_t,t_add;
  int x,y;

  if ((s_t=atan2(s_y-b,s_x-a))<0.0) /* set end angles */
    s_t+=2.0*PI;
  if ((e_t=atan2(e_y-b,e_x-a))<0.0)
    e_t+=2.0*PI;
  if (e_t<=s_t)
    e_t+=2.0*PI;
  t_add=PI/(7.0*radius);
  t=s_t+t_add;

  while (t>s_t && t<e_t) {
    x=(int)rint(radius*cos(t)+a);
    y=(int)rint(radius*sin(t)+b);
    if (x>0 && x<QUAD_X_SIZE && y>0 && y<QUAD_Y_SIZE)
      PutPixel(quad,x,y);
    t+=t_add;
  }
}
