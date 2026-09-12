#include <stdio.h>
#include <xview/xview.h>
#include <xview/canvas.h>
#include <xview/panel.h>
#include <xview/tty.h>
#include <xview/openmenu.h>
#include <xview/seln.h>
#include <string.h>

typedef struct quad_t 
  {
    int ulc_x, ulc_y, lrc_x, lrc_y;
  } quadrant;

#define XHT_FONT       1
#define XHT_SMALL_FONT 2
