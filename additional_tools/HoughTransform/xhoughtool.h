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
 * File:    xhoughtool.h
 * Purpose: header include for the XHoughtool
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include <xview/xview.h>
#include <xview/canvas.h>
#include <xview/font.h>
#include <xview/cms.h>
#include <xview/xv_xrect.h>

#define MAX_COLS 1024 /* Maximum size of the drawing canvases */
#define MAX_ROWS 1024

#define XHT_FONT		1
#define XHT_SMALL_FONT		2

#define MID_GRAY_LEVEL		999

#define NO_ACTIVE_METHOD	0
#define SHT_LINE_METHOD		1
#define RHT_LINE_METHOD		2
#define DRHT_LINE_METHOD	3
#define WRHT_LINE_METHOD	4
#define RWRHT_LINE_METHOD	5
#define CFHT_LINE_METHOD	6
#define PROBHT_LINE_METHOD	7
#define AHT_LINE_METHOD		8
#define CHT_LINE_METHOD		9
#define DCHT_LINE_METHOD	10

#define QUAD_X_SIZE		256
#define QUAD_Y_SIZE		256

typedef struct quad_t {
  int ulc_x, ulc_y, lrc_x, lrc_y;
} quadrant;
