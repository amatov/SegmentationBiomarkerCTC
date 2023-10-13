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
 * File:    ht_hough.h
 * Purpose: header include for hough algorithms
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <values.h>

#ifndef TRUE
#define TRUE			1
#endif
#ifndef FALSE
#define FALSE			0
#endif
#ifndef NULL
#define NULL			0
#endif

#ifdef 	HUGE
#define	INF_DOUBLE_VAL		HUGE
#else
#ifdef 	HUGE_VAL
#define	INF_DOUBLE_VAL		HUGE_VAL
#else
#define	INF_DOUBLE_VAL		(infinity())
#endif
#endif

#define MAX_PIX_VAL		255
#define OBJECT_PIX_VAL		255

#ifdef 	BIG_IMAGES
#define MAX_SIZE		512 /* for 512x512 images */
#else
#define MAX_SIZE		256 /* for 256x256 images */
#endif

#define MAX_FIT_WIN_SIZE	MAX_SIZE/2
#define ACC_MAX_SIZE		MAX_SIZE
#define MAX_VEC_COORD		MAX_SIZE*MAX_SIZE

#define MAX_TIME		100000
#define MAX_LINES		1000

#define OBJECT_GRAY_LEVEL	255

#define PI			3.14159265

typedef int pic_type;

typedef struct picvec
{
  int i,j,label;
} pic_vec_type;

typedef struct coord
{
  int x,y;
} coordinates;

typedef struct coord_t
{
  double X,Y;
  int Value;
} Coord;

typedef struct coord3D_t
{
  double X,Y,Z;
  int Value;
} Coord3D;

typedef struct lineseg
{
  int name;
  int start_x, start_y;
  int start_label;
  int end_x,end_y;
  int end_label;
  int removed;
  coordinates *edge_points;
  struct lineseg *next;
} line;

typedef struct circleseg
{
  int name;
  int start_x, start_y;
  int start_label;
  int end_x, end_y;
  int end_label;
  double centre_x, centre_y;
  double radius;
  int removed;
  coordinates *edge_points;
  struct circleseg *next;
} circle;
