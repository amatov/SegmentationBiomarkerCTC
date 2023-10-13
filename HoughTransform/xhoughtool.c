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
 * File:    	xhoughtool.c
 * Purpose: 	graphical user interface for Hough Transform algorithms
 * Date:    	Jun 1, 1993
 * Last change: Dec 10, 1995
 *****************************************************************************/

#include "ht_hough.h"
#include "xhoughtool.h"
#include "ht_formats.h"

#include <malloc.h>
#include <string.h>
#include <strings.h>
#include <sys/times.h>

#include <X11/X.h>
#include <X11/Xlib.h>

#include <xview/panel.h>
#include <xview/tty.h>
#include <xview/openmenu.h>
#include <xview/seln.h>
#include <xview/scrollbar.h>
#include <xview/notice.h>
#include <xview/notify.h>
#include <xview/svrimage.h>
#include <xview/icon.h>

#include <xview/sel_attrs.h>
#include <xview/textsw.h>
#include <xview/win_input.h>
#include <xview/win_event.h>

extern int errno;
extern char *sys_errlist[];
extern int sys_nerr;


char proc_name[256];

char 	**definition;
char 	**topic;

Display *dpy;
XID xid_canvas, xid_graph_canvas;

Notify_value do_null_processing();
struct itimerval timer;

GC gc, norm_gc, inv_gc;
XGCValues gcvalues;
Cms cms, norm_cms, inv_cms;
Xv_singlecolor xv_color [129]; /* 128 graylevels in use */
unsigned long *pixel_table, *norm_pixel_table, *inv_pixel_table, bg_level,
              fg_level;
int inversed_graylevels=FALSE, grid_graylevel=0;
Textsw 	textsw2;
Xv_Font font, small_font;
XFontStruct *cur_font;
char  font_name[]="-adobe-times-medium-r-normal--10-100-75-75-p-54-iso8859-1",
small_font_name[]="-adobe-times-medium-r-normal--8-80-75-75-p-44-iso8859-1";

short closed_bits[]= { /* XView type monochrome icon */
#include "xhoughtool.icon"
};

Frame frame, popup_frame, popup_help_frame;

Canvas canvas, graph_canvas;

Panel control_panel, popup_panel,popup_help_panel;

Panel_item
	fname_item, message_item, /* control panel items */
	/* panel items for all popups */
	File_item, Format_item, PS_item, Shrink_item, NumOfMaxs_item,
	Threshold_item, MinSegLen_item, MaxSegWidth_item, MaxSegGap_item,
	Noise_item, Save_item, Ask_item, Param_item,
	/* panel items for SHT dialog popup */
	Dim_rho_item, Dim_theta_item,
	/* panel items for RHT dialog popup */
	MinDist_item, MaxDist_item, Accuracy_item,
	/* panel items for DRHT dialog popup */
	Accuracy1_item, Threshold1_item, BlockWidth_item, MaxVar_item,
	/* panel items for WRHT dialog popup */
	ThreshForPoints_item, TolForFitting_item, WinSize_item, DontRmAccu_item,
	ConnectiveCheck_item, Sectoring_item,
	/* panel items for RWRHT dialog popup */
	MinWinSize_item, MaxWinSize_item, AccumulationsLimit_item,
	/* panel items for CFHT dialog popup */
	Msize_item, Epsilon_item, Gamma_item, DontRemovePoints_item,
	/* panel items for ProbHT dialog popup */
	SampleLevel_item,
	/* panel items for AHT dialog popup */
	AccuRangeSel_item, Dim_m_item, Dim_c_item, BinarizeLevel_item,
	Accur_m_item, Accur_c_item, 
	/* panel items for CHT dialog popup */
	Segments_item, Overlap_item, MaskSize_item,
	Help_item,Help_item_list;

Menu props_menu, show_menu, clear_menu, non_prob_menu, prob_menu;

quadrant quads[] = {
	{  0,   0, 255, 255}, {257,   0, 512, 255}, {514,   0, 769, 255},
	{  0, 257, 255, 512}, {257, 257, 512, 512}, {514, 257, 769, 512} };

quadrant active_quad = {  0,   0, 255, 255};

quadrant graph_quads[] = {
	{  0,   0, 254, 200}, {  0, 204, 254, 404}, {  0, 408, 254, 608}};
/*
quadrant quads[] = {
	{  0,   0, 255, 255}, {256,   0, 511, 255}, {512,   0, 767, 255},
	{  0, 256, 255, 511}, {256, 256, 511, 511}, {512, 256, 767, 511} };

quadrant active_quad = {  0,   0, 255, 255};

quadrant graph_quads[] = {
	{  0,   0, 256, 200}, {  0, 205, 256, 400}, {  0, 409, 256, 600}};
*/
void repaint_proc(), panel_repaint_proc(), StopTraining(), prepare_for_running(),
     show_results(), adjust_speed(), change_colormap(), change_grid_gl(),
     warn_prnt(),help_topic_selection();

clear_sel_proc(), clear_graph_proc(), clear_all_proc(), pause_proc(),
start_hough(), show_norm_proc(), show_shrink_proc(), quit_proc(),
select_active_quad_proc(), cycle_event_proc(),

select_sht_line(), select_rht_line(), select_drht_line(),
select_wrht_line(), select_rwrht_line(), select_cfht_line(),
select_probht_line(), select_aht_line(), select_cht_line(),
select_dcht_line(),
select_help(),
done_sht_line_proc(), done_rht_line_proc(), done_drht_line_proc(),
done_wrht_line_proc(), done_rwrht_line_proc(), done_cfht_line_proc(),
done_probht_line_proc(), done_aht_line_proc(), done_cht_line_proc(),
done_dcht_line_proc(),
close_popup_proc(), cancel_popup_proc(), save_image(), save_parameters(),
no_active_method_proc(); /* for procedures not yet implemented */

(*hough_proc)();     /* points to a done_xxxxx_proc() */

int StopFlag=0, StoppedForTraining=0;

/*** Parameter definitions for methods ***/
char Filename[256],  Filename2[256];
int PostScript=0, Shrink=0, NumOfMaxs=999, MinSegLen=10, MaxSegWidth=2,
    MaxSegGap=5, OutputFormat=1;
double real_params[2][2], Noise=0.0;
struct sht_private_parameters { /* parameter structures and default values */
  pic_type Threshold;
  int Dim_rho, Dim_theta;
} sht_param = {10,ACC_MAX_SIZE, ACC_MAX_SIZE};
struct rht_private_parameters {
  int MinDist, MaxDist, Threshold;
  double Accuracy;
} rht_param = {5, 30, 2, 0.01};
struct drht_private_parameters {
  int MinDist, MaxDist, Threshold, Threshold1, BlockWidth;
  double Accuracy, Accuracy1, MaxVar;
} drht_param = {5, 30, 2, 4, 3, 0.02, 0.005, 3.0};
struct wrht_private_parameters {
  int Threshold, ThreshForPoints, WinSize, DontRmAccu, ConnectiveCheck,
      Sectoring;
  double Accuracy, TolForFitting;
} wrht_param = {1, 8, 20, TRUE, TRUE, TRUE, 0.01, 0.25};
struct rwrht_private_parameters {
  int MinDist, MaxDist, Threshold, MinWinSize, MaxWinSize, AccumulationsLimit;
  double Accuracy;
} rwrht_param = {5, 30, 2, 20, 40, FALSE, 0.01};
struct cfht_private_parameters {
  int Msize, ThreshForPoints, Threshold, DontRemovePoints;
  double TolForFitting, Epsilon, Gamma;
} cfht_param = {6, 6, 2, FALSE, 0.2, 1.0, 0.5};
struct probht_private_parameters {
  pic_type Threshold;
  int Dim_rho, Dim_theta;
  double SampleLevel;
} probht_param = {2, ACC_MAX_SIZE, ACC_MAX_SIZE, 20.0};
struct aht_private_parameters {
  int AccuRangeSel, Dim_m, Dim_c;
  double BinarizeLevel, Accur_m, Accur_c;
} aht_param = {1, 9, 9, 0.9, 0.01, 1.0};
struct cht_private_parameters {
  int Threshold, Segments, Overlap, Dim_rho, Dim_theta, MaskSize;
} cht_param = {2, 16, 5, ACC_MAX_SIZE, ACC_MAX_SIZE, 3};
struct dcht_private_parameters {
  int Threshold, Dim_theta;
} dcht_param = {6, ACC_MAX_SIZE};

pic_type pic[MAX_SIZE][MAX_SIZE], pic1[MAX_SIZE][MAX_SIZE],
         gpic[MAX_SIZE][MAX_SIZE], rpic[MAX_SIZE][MAX_SIZE]; /* static pics */
int dim_x, dim_y; /* size of the pic */
int  SavePic=0, Ask=1; /* Parameters for save menu */
char ParametersOut[256]=NULL;

void show_usage(), print_options();

/*****************************************************************************/
main(argc, argv)
int argc;
char *argv[];
{
  extern char *optarg;
  extern int optind;
  char *cptr=rindex(argv[0],'/');
  Server_image closed_window;
  Icon icon;
  int i,c;

  strcpy(proc_name,cptr?++cptr:argv[0]);

  if (argc>5) {
    fprintf(stderr,"Usage: %s [-f Edge_pic] [-I] [-G]\n",proc_name);
    exit(1);
  }

  for (i=1; i<argc; i++)
    if (*argv[i]!='-') {
      if (strcmp(argv[i-1], "-f")) {
        strcpy(Filename, argv[i]);
        strcpy(argv[i], "-n"); /* bluffing getopt() */
      }
      break;
    }

  while ((c = getopt(argc, argv, "f:IGhnu")) != -1)
    switch (c) {
      case 'f': strcpy(Filename, optarg); 
                break;
      case 'I': inversed_graylevels=TRUE; 
                break;
      case 'G': grid_graylevel=255;
      case 'h': print_options();
                exit(0);
      case 'n': break; /* bluff */
      case 'u': 
      default : show_usage();
                exit(0);
    }

  xv_init(XV_INIT_ARGC_PTR_ARGV, &argc, argv, NULL);

  frame = xv_create(NULL, FRAME, XV_LABEL, "XHoughtool 1.1", NULL);
  closed_window=(Server_image)xv_create(NULL, SERVER_IMAGE,
	XV_WIDTH,		64,
	XV_HEIGHT,		64,
	SERVER_IMAGE_BITS,	closed_bits,
	NULL);
  icon=(Icon)xv_create(frame, ICON,
	ICON_IMAGE,		closed_window,
	XV_X,			64,
	XV_Y,			64,
	NULL);
  xv_set(frame, FRAME_ICON, icon, NULL);
  set_up_luts();
  init_graph_canvas();
  init_control_panel();
  init_display_canvas();
  
  window_fit(frame);

  dpy = (Display *)xv_get(frame, XV_DISPLAY); /* connection to display */

  font = xv_find(frame, FONT, FONT_NAME, font_name, NULL); /* set fonts */
  if(!font) {
    ErrMesg(frame,"Sorry, Font not found... Using default font.");
    font = xv_get(frame, XV_FONT);
  }
  small_font = xv_find(graph_canvas, FONT, FONT_NAME, small_font_name, NULL);
  if(!small_font) {
    ErrMesg(frame,"Sorry, Small font not found... Trying with other fonts.");
    small_font = xv_find(frame, FONT, FONT_NAME, font_name, NULL);
    if (!small_font)
      small_font = xv_get(frame, XV_FONT);
  }
  cur_font = (XFontStruct *)xv_get(font, FONT_INFO);
  gcvalues.font = cur_font->fid;
  gcvalues.graphics_exposures = False;
  fg_level = gcvalues.foreground = WhitePixel(dpy, DefaultScreen(dpy));
  bg_level = gcvalues.background = BlackPixel(dpy, DefaultScreen(dpy));
  gc = norm_gc = XCreateGC(dpy, RootWindow(dpy, DefaultScreen(dpy)),
                 GCForeground | GCBackground | GCFont | GCGraphicsExposures,
                 &gcvalues);
  gcvalues.foreground = BlackPixel(dpy, DefaultScreen(dpy));
  gcvalues.background = WhitePixel(dpy, DefaultScreen(dpy));
  inv_gc = XCreateGC(dpy, RootWindow(dpy, DefaultScreen(dpy)),
                 GCForeground | GCBackground | GCFont | GCGraphicsExposures,
                 &gcvalues);

  if (inversed_graylevels) {
    fg_level = BlackPixel(dpy, DefaultScreen(dpy));
    bg_level = WhitePixel(dpy, DefaultScreen(dpy));
    gc=inv_gc;
  }

  hough_proc = no_active_method_proc; /* current running method */

	init_help();
	
  xv_main_loop(frame);
  while(StoppedForTraining) {   
    /*  start the work loop, and do notify_dispatch there every now and then */
    prepare_for_running();
    (*hough_proc)(); /* run the selected method */
    show_results();
    StoppedForTraining = 0;
    XFlush(dpy);
    notify_start();
  }

  exit(0);
}

/*****************************************************************************/
set_up_luts()
/*** Setting up color look up tables for canvases ***/
{
  int i;

  for (i=0; i<128; i++)
    xv_color[i].red=xv_color[i].green=xv_color[i].blue=i*2+1;

  cms = norm_cms = (Cms)xv_create(NULL, CMS,
	CMS_SIZE,		128, /* asking for 128 graylevels */
	CMS_NAME,		"images",
	CMS_COLORS,		xv_color,
	NULL);

  pixel_table=norm_pixel_table=(unsigned long *)xv_get(norm_cms, CMS_INDEX_TABLE);

  for (i=0; i<128; i++)
    xv_color[i].red=xv_color[i].green=xv_color[i].blue=255-(i*2);

  inv_cms = (Cms)xv_create(NULL, CMS,
	CMS_SIZE,		128, /* asking for 128 graylevels */
	CMS_NAME,		"images",
	CMS_COLORS,		xv_color,
	NULL);

  inv_pixel_table=(unsigned long *)xv_get(inv_cms, CMS_INDEX_TABLE);

  if (inversed_graylevels) {
    cms=inv_cms;
    pixel_table=inv_pixel_table;
  }
}

/*****************************************************************************/
init_graph_canvas()
/*** Left hand side drawing canvas initialization ***/
{
  graph_canvas = xv_create(frame, CANVAS,
  	CANVAS_AUTO_SHRINK,	FALSE,
	WIN_X,			0,
	WIN_Y,			0,
	XV_WIDTH,		256,
	XV_HEIGHT,		613,
	WIN_DYNAMIC_VISUAL,	TRUE,
	WIN_CMS,		cms,
	CANVAS_WIDTH,		MAX_COLS, /* 1024 */
	CANVAS_HEIGHT,		MAX_ROWS, /* 1024 */
        CANVAS_X_PAINT_WINDOW,	TRUE,
	CANVAS_REPAINT_PROC,	repaint_proc,
	NULL);

  xid_graph_canvas = (XID)xv_get(canvas_paint_window(graph_canvas), XV_XID);
}

/*****************************************************************************/
init_display_canvas()
/*** Main drawing canvas initialization ***/
{
  int retained,i;

  canvas = xv_create(frame, CANVAS,
  	CANVAS_AUTO_SHRINK,	FALSE,
	/* WIN_X,		256,		
	WIN_Y,			144, */
	WIN_RIGHT_OF,		graph_canvas, 
	WIN_BELOW,		control_panel, 
	XV_WIDTH,		770,
	XV_HEIGHT,		513,
	WIN_DYNAMIC_VISUAL,	TRUE,
        WIN_CMS,		cms, 
	CANVAS_WIDTH,		MAX_COLS, /* 1024 */
	CANVAS_HEIGHT,		MAX_ROWS, /* 1024 */
        CANVAS_X_PAINT_WINDOW,	TRUE,
        CANVAS_REPAINT_PROC,	repaint_proc,
	NULL);
/*
  xv_create(canvas, SCROLLBAR,
	SCROLLBAR_DIRECTION,	SCROLLBAR_HORIZONTAL,
	NULL);

  xv_create(canvas, SCROLLBAR,
	SCROLLBAR_DIRECTION,	SCROLLBAR_VERTICAL,
	NULL);
*/
  xid_canvas = (XID)xv_get(canvas_paint_window(canvas), XV_XID);
}

/*****************************************************************************/
init_control_panel()
/*** Control panel initialization ***/
{
  char *getwd();

  if (!strlen(Filename)) {
    getwd(Filename); /* get current working directory pathname */
    strcat(Filename,"/");
  }

  control_panel = xv_create(frame, PANEL,
	WIN_RIGHT_OF,		graph_canvas,
	WIN_Y,			0,
	XV_HEIGHT,		100,
	XV_WIDTH,		770,
	PANEL_REPAINT_PROC,	panel_repaint_proc,
	NULL);

  /*** Control panel objects initialization ***/
  fname_item = xv_create(control_panel, PANEL_TEXT,
	XV_X,			xv_col(control_panel,0),
	XV_Y,			xv_row(control_panel,0),
	PANEL_VALUE_DISPLAY_LENGTH, 67, /* 67 , original 68 */
	PANEL_LABEL_STRING,	"Image file path: ",
	PANEL_VALUE,		Filename,
	NULL);

  /*** Non-Probabil. and Probabilistic buttons hide following menus ***/
  non_prob_menu = xv_create(NULL, MENU,
	MENU_TITLE_ITEM,	"Non-probabilistic methods",
	MENU_ITEM,
		MENU_STRING,	"SHT for line (rho,theta)",
		MENU_NOTIFY_PROC, select_sht_line,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"CFHT for line y=ax+b",
		MENU_NOTIFY_PROC, select_cfht_line,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"AHT for line y=mx+c",
		MENU_NOTIFY_PROC, select_aht_line,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"CHT for line (rho,theta)",
		MENU_NOTIFY_PROC, select_cht_line,
		NULL,
	NULL); 

  xv_create(control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Non-probabil.",
	XV_X,			xv_col(control_panel,0),
	XV_Y,			xv_row(control_panel,1),
	PANEL_ITEM_MENU,	non_prob_menu,
	NULL);

  prob_menu = xv_create(NULL, MENU,
	MENU_TITLE_ITEM,	"Probabilistic methods",
	MENU_ITEM,
		MENU_STRING,	"RHT for line y=ax+b",
		MENU_NOTIFY_PROC, select_rht_line,
		MENU_PULLRIGHT,
			xv_create(NULL, MENU,
				MENU_TITLE_ITEM,	"RHT methods",
				MENU_ITEM,
					MENU_STRING,	"RHT for line y=ax+b",
					MENU_NOTIFY_PROC, select_rht_line,
					NULL,
				MENU_ITEM,
					MENU_STRING,	"DRHT for line y=ax+b",
					MENU_NOTIFY_PROC, select_drht_line,
					NULL,
				MENU_ITEM,
					MENU_STRING,	"WRHT (CRHT) for line y=ax+b",
					MENU_NOTIFY_PROC, select_wrht_line,
					NULL,
				MENU_ITEM,
					MENU_STRING,	"RWRHT for line y=ax+b",
					MENU_NOTIFY_PROC, select_rwrht_line,
					NULL,
				NULL),
		NULL,
	MENU_ITEM,
		MENU_STRING,	"ProbHT for line (rho, theta)",
		MENU_NOTIFY_PROC, select_probht_line,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"DCHT for line (rho, theta)",
		MENU_NOTIFY_PROC, select_dcht_line,
		NULL,
	NULL); 

  xv_create(control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Probabilistic",
	/*XV_X,			xv_col(control_panel,16),*/
	XV_Y,			xv_row(control_panel,1),
	PANEL_ITEM_MENU,	prob_menu,
	NULL);

  xv_create(control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Start",
	/*XV_X,			xv_col(control_panel,31),*/
	XV_Y,			xv_row(control_panel,1),
	PANEL_NOTIFY_PROC,	start_hough,
	NULL);

  xv_create(control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Pause",
	/*XV_X,			xv_col(control_panel,38),*/
	XV_Y,			xv_row(control_panel,1),
	PANEL_NOTIFY_PROC,	pause_proc,
	NULL);

  xv_create(control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Stop",
	/*XV_X,			xv_col(control_panel,46),*/
	XV_Y,			xv_row(control_panel,1),
	PANEL_NOTIFY_PROC,	StopTraining,
	PANEL_CLIENT_DATA,	frame, 
	NULL);

  xv_create(control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Quit",
	/* XV_X,			xv_col(control_panel,107), */
	XV_Y,			xv_row(control_panel,1),
	PANEL_NOTIFY_PROC,	quit_proc,
	NULL);

  xv_create(control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Help",
	/* XV_X,			xv_col(control_panel,107), */
	XV_Y,			xv_row(control_panel,1),
	PANEL_NOTIFY_PROC,	select_help,
	PANEL_CLIENT_DATA,	frame, 
	NULL);

  xv_create(control_panel, PANEL_SLIDER,
	PANEL_LABEL_STRING,	" Wait states : ",
	/*XV_X,			xv_col(control_panel,61),*/
	XV_Y,			xv_row(control_panel,1),
	PANEL_VALUE,		0,
	PANEL_MAX_VALUE,	60,
	PANEL_SLIDER_WIDTH,	xv_get(control_panel, PANEL_ITEM_X)>400?10:80,
	PANEL_NOTIFY_PROC,	adjust_speed,
	NULL);

  props_menu = xv_create(NULL, MENU, 
	MENU_TITLE_ITEM,	"Properties",
	MENU_ITEM,
		MENU_STRING,	"Background: white/black ",
		MENU_NOTIFY_PROC, change_colormap,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"Grid: enable/disable",
		MENU_NOTIFY_PROC, change_grid_gl,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"Save image",
		MENU_NOTIFY_PROC, save_image,
		NULL,
	NULL);

  xv_create(control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Props",
	XV_X,			xv_col(control_panel,0), 
	XV_Y,			xv_row(control_panel,2),
	PANEL_ITEM_MENU,	props_menu,
	NULL);

  show_menu = xv_create(NULL, MENU, 
	MENU_TITLE_ITEM,	"Show",
	MENU_ITEM,
		MENU_STRING,	"Show in normal way",
		MENU_NOTIFY_PROC, show_norm_proc,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"Show pic as shrunken",
		MENU_NOTIFY_PROC, show_shrink_proc,
		NULL,
	NULL);

  xv_create( control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Show",
	/* XV_X,		xv_col(control_panel,0),*/
	XV_Y,			xv_row(control_panel,2),
	PANEL_ITEM_MENU,	show_menu,
	NULL);

  clear_menu = xv_create(NULL, MENU,
	MENU_TITLE_ITEM,	"Clear",
	MENU_ITEM,
		MENU_STRING,	"Clear selected place",
		MENU_NOTIFY_PROC, clear_sel_proc,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"Clear graph canvas",
		MENU_NOTIFY_PROC, clear_graph_proc,
		NULL,
	MENU_ITEM,
		MENU_STRING,	"Clear all",
		MENU_NOTIFY_PROC, clear_all_proc,
		NULL,
	NULL);

  xv_create( control_panel, PANEL_BUTTON,
	PANEL_LABEL_STRING,	"Clear",
	/*XV_X,			xv_col(control_panel,9),*/
	XV_Y,			xv_row(control_panel,2),
	PANEL_ITEM_MENU,	clear_menu,
	NULL);

  xv_create(control_panel, PANEL_CHOICE,
	/*XV_X,			xv_col(control_panel,19),*/
	XV_Y,			xv_row(control_panel,2),
	PANEL_LABEL_STRING,	" Place : ",
	PANEL_CHOICE_STRINGS,	"1.1", "1.2", "1.3", "2.1", "2.2", "2.3", 0,
	PANEL_NOTIFY_PROC,	select_active_quad_proc, 
	PANEL_FEEDBACK,		PANEL_INVERTED,
	NULL);

  /*** Text space for name of the active method ***/
  message_item = xv_create(control_panel, PANEL_MESSAGE,
	/*XV_X,			xv_col(control_panel,54),*/
	XV_Y,			xv_row(control_panel,2),
	PANEL_LABEL_STRING,	"  No active method",
	XV_SHOW,		TRUE,
	NULL);

  /*window_fit_height(control_panel);*/
}


/*****************************************************************************
	Functions for asking parameter values

 *****************************************************************************/

char value_str[32];

/*****************************************************************************/
Panel_item ask_int_value(panel,col,row,size,label_str,value)
Panel panel;
int col,row,size,value;
char *label_str;
{
  sprintf(value_str,"%d",value);
  return (Panel_item)xv_create(panel,PANEL_TEXT,
	XV_X,			xv_col(panel,col),
	XV_Y,			xv_row(panel,row),
	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT,
	PANEL_VALUE_DISPLAY_LENGTH, size,
	PANEL_LABEL_STRING,	label_str,
	PANEL_VALUE,		value_str,
	NULL);
}

/*****************************************************************************/
Panel_item ask_double_value(panel,col,row,size,label_str,value)
Panel panel;
int col,row,size;
char *label_str;
double value;
{
  sprintf(value_str,"%8.4lf",value);
  return (Panel_item)xv_create(panel, PANEL_TEXT,
	XV_X,			xv_col(panel,col),
	XV_Y,			xv_row(panel,row),
	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT,
	PANEL_VALUE_DISPLAY_LENGTH, size,
	PANEL_LABEL_STRING,	label_str,
	PANEL_VALUE,		value_str,
	NULL);
}

/*****************************************************************************/
Panel_item ask_str_value(panel,col,row,size,label_str,value)
Panel panel;
int col,row,size;
char *label_str, *value;
{
  return (Panel_item)xv_create(panel, PANEL_TEXT,
	XV_X,			xv_col(panel,col),
	XV_Y,			xv_row(panel,row),
	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT,
	PANEL_VALUE_DISPLAY_LENGTH, size,
	PANEL_LABEL_STRING,	label_str,
	PANEL_VALUE,		value,
	NULL);
}

/*****************************************************************************/
Panel_item ask_cycle_value(panel,col,row,label_str,value)
Panel panel;
int col,row,value;
char *label_str;
{
  return (Panel_item)xv_create(panel, PANEL_CYCLE,
	XV_X,			xv_col(panel,col),
	XV_Y,			xv_row(panel,row),
	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT,
	PANEL_LABEL_STRING,	label_str,
	PANEL_CHOICE_STRINGS,	"off","on", NULL,
	PANEL_EVENT_PROC,	cycle_event_proc,
	PANEL_VALUE,		value,
	NULL);

}

/*****************************************************************************/
Panel_item ask_format_value(panel,col,row,label_str,value)
Panel panel;
int col,row,value;
char *label_str;
{
  return (Panel_item)xv_create(panel, PANEL_CHOICE,
	XV_X,			xv_col(panel,col),
	XV_Y,			xv_row(panel,row),
	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT,
	PANEL_LABEL_STRING,	label_str,
	PANEL_CHOICE_STRINGS,	" CVL ", " PGM ",
				" SKE ", " BIN ", NULL,
	PANEL_FEEDBACK,		PANEL_INVERTED,
	PANEL_VALUE,		value,
	NULL);
}

/*****************************************************************************/
Panel_item ask_topic_value(panel,col,row,label_str,value)
Panel panel;
int col,row,value;
char *label_str;
{
  return (Panel_item)xv_create(panel, PANEL_CHOICE,
	XV_X,			xv_col(panel,col),
	XV_Y,			xv_row(panel,row),
	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT,
	PANEL_LABEL_STRING,	label_str,
	PANEL_CHOICE_STRINGS,	" Introduction ", " Command Line Options ",
				" Interface ", " Save ", NULL,
	PANEL_FEEDBACK,		PANEL_INVERTED,
	PANEL_VALUE,		value,
	NULL);
}

/*****************************************************************************/
Panel_item ask_save_value(panel,col,row,label_str,value)
Panel panel;
int col,row,value;
char *label_str;
{
  return (Panel_item)xv_create(panel, PANEL_CHOICE,
	XV_X,			xv_col(panel,col),
	XV_Y,			xv_row(panel,row),
	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT,
	PANEL_LABEL_STRING,	label_str,
	PANEL_CHOICE_STRINGS,	" Result image        ", " Line candidates  ",
				NULL,
	PANEL_FEEDBACK,		PANEL_INVERTED,
	PANEL_VALUE,		value,
	NULL);
}


/*****************************************************************************
	Popup functions for dialog boxes (1 for each method to set parameters)

 *****************************************************************************/

/*****************************************************************************/
create_popup(label_str, method)
char *label_str;
int method;
{
  int row=0;
  char *fname;

  if (popup_frame)
    close_popup_proc();

  popup_frame=xv_create(frame, FRAME_CMD,
	XV_LABEL,		label_str,
        NULL);

  popup_panel=xv_get(popup_frame, FRAME_CMD_PANEL);

  File_item = ask_str_value(popup_panel,0,row++,50,"Edge file: ",Filename);
  Param_item = ask_str_value(popup_panel,0,row++,50,"Parameters save file: ",
	ParametersOut);

  Noise_item = ask_double_value(popup_panel,0,row,9,
                                     "Randomly added noise (%): ",Noise);
/*
  Format_item = ask_format_value(popup_panel,0,row++,"Output image format: ",
                                     OutputFormat);
  PS_item = ask_cycle_value(popup_panel,0,row,"PostScript: ",PostScript);
*/
  Shrink_item = ask_cycle_value(popup_panel,39,row++,"Vertical shrink: ",Shrink);
  NumOfMaxs_item = ask_int_value(popup_panel,0,row++,5,
                                    "Number of maximas (posit int): ",
                                     NumOfMaxs);
  MinSegLen_item = ask_int_value(popup_panel,0,row++,5,
                                     "Minimum line segment length (posit int): ",
                                     MinSegLen);
  MaxSegWidth_item = ask_int_value(popup_panel,0,row++,5,
                                     "Line scanning width (non-negat int): ",
                                     MaxSegWidth);
  MaxSegGap_item = ask_int_value(popup_panel,0,row++,5,
                  "Maximum gap between pixels in line segment (non-negat int): ",
                                     MaxSegGap);
  switch (method) { /* method depended parameter value questions */
    case NO_ACTIVE_METHOD:	break; /* should never happen */
    case SHT_LINE_METHOD:	row++;
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                   "Minimum score accepted for an accumulator max. (posit int): ",
                                     sht_param.Threshold);
      Dim_rho_item = ask_int_value(popup_panel,0,row++,5,
                                   "Accumulator rho quantization (posit int): ",
                                   sht_param.Dim_rho);
      Dim_theta_item = ask_int_value(popup_panel,0,row++,5,
                                 "Accumulator theta quantization (posit int): ",
                                   sht_param.Dim_theta);
				break;
    case RHT_LINE_METHOD:	row++;
      MinDist_item = ask_int_value(popup_panel,0,row++,5,
                                     "Min. dist. for point pair (posit int): ",
                                     rht_param.MinDist);
      MaxDist_item = ask_int_value(popup_panel,0,row++,5,
                        "Max. dist. for point pair (posit int, >=min dist.): ",
                                     rht_param.MaxDist);
      Accuracy_item = ask_double_value(popup_panel,0,row++,9,
                                     "Accumulator accuracy (float, [0, 1]): ",
                                     rht_param.Accuracy);
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                                     "Accumulator threshold (posit int): ",
                                     rht_param.Threshold);
				break;
    case DRHT_LINE_METHOD:	row++;
      MinDist_item = ask_int_value(popup_panel,0,row++,5,
                                     "Min. dist. for point pair (posit int): ",
                                     drht_param.MinDist);
      MaxDist_item = ask_int_value(popup_panel,0,row++,5,
                        "Max. dist. for point pair (posit int, >=min dist.): ",
                                     drht_param.MaxDist);
      Accuracy_item = ask_double_value(popup_panel,0,row++,9,
                                     "Accumulator accuracy (float, [0, 1]): ",
                                     drht_param.Accuracy);
      Accuracy1_item = ask_double_value(popup_panel,0,row++,9,
                                    "2nd Accumulator accuracy (float, [0, 1]): ",
                                     drht_param.Accuracy1);
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                                     "Accumulator threshold (posit int): ",
                                     drht_param.Threshold);
      Threshold1_item = ask_int_value(popup_panel,0,row++,5,
                                     "2nd Accumulator threshold (posit int): ",
                                     drht_param.Threshold1);
      BlockWidth_item = ask_int_value(popup_panel,0,row++,5,
                                    "Block width for 2nd iteration (posit int): ",
                                     drht_param.BlockWidth);
      MaxVar_item = ask_double_value(popup_panel,0,row++,9,
                          "Max. variation of 'a' in degrees (float, [0, 1]): ",
                                   drht_param.MaxVar);
				break;
    case WRHT_LINE_METHOD:	row++;
      Accuracy_item = ask_double_value(popup_panel,0,row++,9,
                                     "Accumulator accuracy (float, [0, 1]): ",
                                     wrht_param.Accuracy);
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                                     "Accumulator threshold (non-negat int): ",
                                     wrht_param.Threshold);
      WinSize_item = ask_int_value(popup_panel,0,row++,5,
                                    "Window size (posit int, (WinSize*2+1)^2): ",
                                     wrht_param.WinSize);
      ThreshForPoints_item = ask_int_value(popup_panel,0,row++,5,
                       "Threshold for min. num. of windowed points (posit int): ",
                                     wrht_param.ThreshForPoints);
      TolForFitting_item = ask_double_value(popup_panel,0,row++,9,
                               "Tolerance for fitting error (float, [0, 1]): ",
                                     wrht_param.TolForFitting);
      DontRmAccu_item = ask_cycle_value(popup_panel,0,row++,
                                    "Do not remove accumulator after line found: ",
                                     wrht_param.DontRmAccu);
      ConnectiveCheck_item = ask_cycle_value(popup_panel,0,row++,
                             "Use only center point connected points in window (CRHT): ",
                                     wrht_param.ConnectiveCheck);
      Sectoring_item = ask_cycle_value(popup_panel,0,row++,
                                   "Do the connective point search as sectored: ",
                                     wrht_param.Sectoring);
				break;
    case RWRHT_LINE_METHOD:	row++;
      MinDist_item = ask_int_value(popup_panel,0,row++,5,
                                     "Min. dist. for point pair (posit int): ",
                                     rwrht_param.MinDist);
      MaxDist_item = ask_int_value(popup_panel,0,row++,5,
                        "Max. dist. for point pair (posit int, >=min dist.): ",
                                     rwrht_param.MaxDist);
      Accuracy_item = ask_double_value(popup_panel,0,row++,9,
                                     "Accumulator accuracy (float, [0, 1]): ",
                                     rwrht_param.Accuracy);
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                                     "Accumulator threshold (posit int): ",
                                     rwrht_param.Threshold);
      MinWinSize_item = ask_int_value(popup_panel,0,row++,5,
                             "Min. window size (posit int, (WinSize*2+1)^2): ",
                                     rwrht_param.MinWinSize);
      MaxWinSize_item = ask_int_value(popup_panel,0,row++,5,
                             "Max. window size (posit int, (WinSize*2+1)^2): ",
                                     rwrht_param.MaxWinSize);
      AccumulationsLimit_item = ask_cycle_value(popup_panel,0,row++,
                                   "Accumulate parameters in a win max 20T times: ",
                                     rwrht_param.AccumulationsLimit);
				break;
    case CFHT_LINE_METHOD:	row++;
      Msize_item = ask_int_value(popup_panel,0,row++,5,
                                    "Window size (posit int, (Msize*2+1)^2): ",
                                     cfht_param.Msize);
      ThreshForPoints_item = ask_int_value(popup_panel,0,row++,5,
                       "Threshold for amount of windowed points (posit int): ",
                                     cfht_param.ThreshForPoints);
      TolForFitting_item = ask_double_value(popup_panel,0,row++,9,
                               "Tolerance for fitting error (float, [0, 1]): ",
                                     cfht_param.TolForFitting);
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                                     "Accumulator threshold (posit int): ",
                                     cfht_param.Threshold);
      Epsilon_item = ask_double_value(popup_panel,0,row++,9,
                    "Tolerance to existing accumulator cells (posit float): ",
                                     cfht_param.Epsilon);
      Gamma_item = ask_double_value(popup_panel,0,row++,9,
                                  "Weight for old score (float, ]0, 1[): ",
                                     cfht_param.Gamma);
      DontRemovePoints_item = ask_cycle_value(popup_panel,0,row++,
                                     "Do not remove edge points: ",
                                     cfht_param.DontRemovePoints);
				break;
    case PROBHT_LINE_METHOD:	row++;
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                  "Minimum accepted score for a accumulator max. (posit int): ",
                                     probht_param.Threshold);
      Dim_rho_item = ask_int_value(popup_panel,0,row++,5,
                                   "Accumulator rho quantization (posit int): ",
                                   probht_param.Dim_rho);
      Dim_theta_item = ask_int_value(popup_panel,0,row++,5,
                                 "Accumulator theta quantization (posit int): ",
                                   probht_param.Dim_theta);
      SampleLevel_item = ask_double_value(popup_panel,0,row++,9,
                                  "Sampling level (%) (float, [0, 100[): ",
                                     probht_param.SampleLevel);
				break;
    case AHT_LINE_METHOD:	row++;
      AccuRangeSel_item = ask_int_value(popup_panel,0,row++,5,
                               "Accumulator range selection (posit int, [1, 2]): ",
                                     aht_param. AccuRangeSel);
      Dim_m_item = ask_int_value(popup_panel,0,row++,5,
                                   "Accumulator slope quantization (posit int): ",
                                   aht_param.Dim_m);
      Dim_c_item = ask_int_value(popup_panel,0,row++,5,
                                "Accumulator intercept quantization (posit int): ",
                                 aht_param.Dim_c);
      BinarizeLevel_item = ask_double_value(popup_panel,0,row++,9,
                                  "Accumulator binarizing level (float, [0, 1[): ",
                                  aht_param.BinarizeLevel);
      Accur_m_item = ask_double_value(popup_panel,0,row++,9,
                                     "Relative slope accuracy (float, [0, 10[): ",
                                     aht_param.Accur_m);
      Accur_c_item = ask_double_value(popup_panel,0,row++,9,
                                     "Relative intercept accuracy (float, [0, 10[): ",
                                     aht_param.Accur_c);
				break;
    case CHT_LINE_METHOD:	row++;
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                   "Minimum score accepted for an accumulator max. (posit int): ",
                                     cht_param.Threshold);
      Segments_item = ask_int_value(popup_panel,0,row++,5,
                   "Number of segments in image (Segments^2, posit int): ",
                                     cht_param.Segments);
      Overlap_item = ask_int_value(popup_panel,0,row++,5,
                   "Number of overlapped pixels in segments (non-negat int): ",
                                     cht_param.Overlap);
      Dim_rho_item = ask_int_value(popup_panel,0,row++,5,
                                   "Accumulator rho quantization (posit int): ",
                                   cht_param.Dim_rho);
      Dim_theta_item = ask_int_value(popup_panel,0,row++,5,
                                 "Accumulator theta quantization (posit int): ",
                                 cht_param.Dim_theta);
      MaskSize_item = ask_int_value(popup_panel,0,row++,5,
                "Local maxima search mask size ((2*MaskSize+1)^2, posit int): ",
                                     cht_param.MaskSize);
				break;
    case DCHT_LINE_METHOD:	row++;
      Threshold_item = ask_int_value(popup_panel,0,row++,5,
                                     "Accumulator threshold (posit int): ",
                                     dcht_param.Threshold);
      Dim_theta_item = ask_int_value(popup_panel,0,row++,5,
                                 "Accumulator theta quantization (posit int): ",
                                 dcht_param.Dim_theta);
				break;
    default:			break; /* probably new method */
  }
/*
  xv_create(popup_panel, PANEL_BUTTON,
	PANEL_NOTIFY_PROC,	close_popup_proc,
	PANEL_LABEL_STRING,	"Close",
	NULL);
*/
  window_fit(popup_panel);
  window_fit(popup_frame);
  xv_set(popup_frame, WIN_SHOW, TRUE, NULL);
}


/*****************************************************************************/
select_help() /* Make help panel visible */
{
	xv_set(popup_help_frame, WIN_SHOW, TRUE, NULL);
}

/*****************************************************************************/
init_help() 
/* Create help panel and read the help file */
{
char 	help_file[50];		
char 	tmp_line[60];
char 	intend[100]; 
char    xhtoolhelpdir[100]; 
char    err_mess[150];
char 	tmp_explanation[6143];
int 	index=0, top=0, len, intend_index=0, max_intend,row=0; 
int     textsw_width, mem_error=0; 
FILE 	*fp_in;
Panel   panel;

 	if (popup_help_frame)		 
  		close_help();
 
  	popup_help_frame=(Frame)xv_create(NULL, FRAME,
		XV_LABEL, "XHoughtool HELP", NULL);
  	panel = (Panel)xv_create(popup_help_frame, PANEL, NULL);
  	(void) xv_create(panel, PANEL_BUTTON,
        	PANEL_LABEL_STRING,     "Close",
        	PANEL_NOTIFY_PROC,      close_help,
        	NULL);

  	Help_item_list=(Panel_item)xv_create(panel, PANEL_LIST, 
		PANEL_LIST_DISPLAY_ROWS, 13,
		PANEL_NOTIFY_PROC, help_topic_selection,
	   	NULL, NULL);

        /* You could definitely program the following code associated with XHTOOLHELPDIR 
           and error messages better. I am sorry for C programming experts. */
        
        if (getenv("XHTOOLHELPDIR")) { 
            strcpy(xhtoolhelpdir, getenv("XHTOOLHELPDIR"));
            sprintf(help_file, "%s/xhoughtool.help", strtok(xhtoolhelpdir,":"));                      
        }
        else {   
            strcpy(help_file, "xhoughtool.help");
        }
	strcpy(intend, "                                                       ");
	max_intend=strlen(intend);
	if((fp_in=fopen(help_file,"r"))==NULL)
	{       
                sprintf(err_mess,"Sorry, can not find the help file %s ! \n Please define the XHTOOLHELPDIR variable correctly or copy the help file to this directory", help_file); 
                ErrMesg(popup_help_frame,err_mess);  
		} else {
		/* reserve space for help information */
		topic=(char**)malloc((unsigned)(sizeof(long)*100));
		definition=(char**)malloc((unsigned)(sizeof(long)*100));
		if( definition==NULL || topic==NULL ) mem_error=1;
		/* parse through the help file */
		while (fgets(tmp_line,4128,fp_in) && !mem_error ) 
		{
			if(!strncmp(tmp_line,"<PUSH>",6)) {
				if(intend_index+1<max_intend) intend_index+=2;
				continue;
			}
			if(!strncmp(tmp_line,"<POP>",5)) {
				if(intend_index) intend_index-=2;
				continue;
			}
			if(!strncmp(tmp_line,"<LB>",4)) {
				
				strcat(&tmp_explanation[index],"\n");
				index++;
				continue;
			}
			if(!strncmp(tmp_line,"<TOPIC>",7)) {
				topic[top]=(char *)malloc(
				    (unsigned)(sizeof(char)*(strlen(tmp_line)+1)));
    				if(topic[top]==NULL) mem_error=1;
				sprintf(topic[top],"%s%s",
					&intend[max_intend-intend_index],&tmp_line[7]); 
				xv_set(Help_item_list, PANEL_LIST_INSERT,row,
					PANEL_LIST_STRING, row++, topic[top], NULL);
				if(top) {
					if(!strlen(tmp_explanation)) {tmp_explanation[0]=' ';
						tmp_explanation[1]='\0'; }
					definition[top-1]=(char *)malloc(
						(unsigned)(sizeof(char)*(strlen(tmp_explanation)+1)));	
					if(definition[top-1]==NULL) mem_error=1;
					else strcpy(definition[top-1],tmp_explanation);			
				}
    				top++; 
				strcpy(tmp_explanation,"");
 				
    				index=0;
 			} else {
    				len=strlen(tmp_line);
				tmp_line[len-1]=' '; 
				strncat(&tmp_explanation[index],tmp_line,len);
    				index+=len;
 			}
		}
		fclose(fp_in);
		if(!mem_error) {
			definition[top-1]=(char *)malloc(
			(	unsigned)(sizeof(char)*strlen(tmp_explanation)+1));
			if(definition[top-1]==NULL) mem_error=1;	
				strcpy(definition[top-1],tmp_explanation);
		 }
		
		if(mem_error) {
	 		warn_prnt ("Can not allocate memory for help! - help disabled \n");
  			xv_destroy(Help_item_list);
			Help_item_list=(Panel_item)xv_create(panel, PANEL_LIST, 
				PANEL_LIST_DISPLAY_ROWS, 1,
				PANEL_NOTIFY_PROC, help_topic_selection,
	   			NULL, NULL);
			window_fit(panel);
			window_fit(popup_help_frame);
		}
		else {
			window_fit(panel);
			textsw2=(Textsw)xv_create(popup_help_frame, TEXTSW, WIN_ROWS, 20,
           			WIN_COLUMNS, 60, TEXTSW_WRAPAROUND_SIZE, 10000,
				/*TEXTSW_BROWSING, TRUE,*/
				TEXTSW_DISABLE_CD, TRUE,
				TEXTSW_DISABLE_LOAD,TRUE,
				NULL);	
			window_fit(popup_help_frame);
	  		xv_set(popup_help_frame, WIN_SHOW, FALSE, NULL);
		}	
	}
}

/*****************************************************************************/
void help_topic_selection(item, string, client_data,op, event)
/* Show help associated with selected topic */ 
Panel_item item;
char *string;
caddr_t client_data;
Panel_list_op op;
Event *event;
{
int i, top=0; 
	textsw_erase(textsw2, 0, TEXTSW_INFINITY);
	while(topic[top]){
		if(!strcmp(topic[top],string)) break;
		top++;
	} 
	textsw_insert(textsw2,definition[top], strlen(definition[top]));
	xv_set(textsw2,TEXTSW_FIRST,0,NULL);
}

/*****************************************************************************/
select_sht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_sht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  SHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for SHT with parameterization (rho, theta) for line",
               SHT_LINE_METHOD);
}

/*****************************************************************************/
select_rht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_rht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  RHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for RHT with parameterization (a, b) for line",
               RHT_LINE_METHOD);
}

/*****************************************************************************/
select_drht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_drht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  DRHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for DRHT with parameterization (a, b) for line",
               DRHT_LINE_METHOD);
}

/*****************************************************************************/
select_wrht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_wrht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  WRHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for WRHT (CRHT) with parameterization (a, b) for line",
               WRHT_LINE_METHOD);
}

/*****************************************************************************/
select_rwrht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_rwrht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  RWRHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for RWRHT with parameterization (a, b) for line",
               RWRHT_LINE_METHOD);
}

/*****************************************************************************/
select_cfht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_cfht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  CFHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for CFHT with parameterization (a, b) for line",
               CFHT_LINE_METHOD);
}

/*****************************************************************************/
select_probht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_probht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  ProbHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for ProbHT with parameterization (rho, theta) for line",
               PROBHT_LINE_METHOD);
}

/*****************************************************************************/
select_aht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_aht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  AHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for AHT with parameterization (m, c) for line",
               AHT_LINE_METHOD);
}

/*****************************************************************************/
select_cht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_cht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  CHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for CHT with parameterization (rho, theta) for line",
               CHT_LINE_METHOD);
}

/*****************************************************************************/
select_dcht_line()
{
  strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
  hough_proc = done_dcht_line_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  DCHT for line: ready and waiting",
         NULL);
  create_popup("Parameters for DCHT with parameterization (rho, theta) for line",
               DCHT_LINE_METHOD);
}

/*****************************************************************************
	Repaint Procedures

 *****************************************************************************/

/*****************************************************************************/
void repaint_proc(canv, paint_window, rects)
/*** Drawing canvas repaint procedure ***/
Canvas        canv;   
Xv_Window     paint_window;
Rectlist     *rects; 
{
  XSetForeground(dpy,gc,pixel_table[grid_graylevel/2]);
  XDrawLine(dpy, xid_canvas, gc,   0, 256, 769, 256);
  XDrawLine(dpy, xid_canvas, gc, 256,   0, 256, 512);
  XDrawLine(dpy, xid_canvas, gc, 513,   0, 513, 512);
  XDrawLine(dpy, xid_graph_canvas, gc,   0, 202, 256, 202);
  XDrawLine(dpy, xid_graph_canvas, gc,   0, 406, 256, 406);
  XSetForeground(dpy,gc,fg_level);
  XFlush(dpy);
}

/*****************************************************************************/
void panel_repaint_proc(panel, paint_window)
/*** Control panel repaint procedure ***/
Panel         panel;
Xv_Window     paint_window;
{
  strcpy(xv_get(fname_item, PANEL_VALUE), Filename);
  XFlush(dpy);
}

/*****************************************************************************
	Procedure functions for control_panel button and slider actions

 *****************************************************************************/

/*****************************************************************************/
Notify_value do_null_processing()
{
  struct tms start_time, wait_time;

  times(&start_time);

  do { /* (almost) null processing for 1/60 seconds */
    notify_dispatch();
    XFlush(dpy);
    times(&wait_time);
  } while ((wait_time.tms_utime-start_time.tms_utime)<1L);

  return NOTIFY_DONE;
}

/*****************************************************************************/
void adjust_speed(item, value)
Panel_item item;
int value;
{
  if (value) { /* do some null processing every now and then, or even more often */
    timer.it_value.tv_usec=999999/value;
    timer.it_interval.tv_usec=999999/value;
    notify_set_itimer_func(frame, do_null_processing,
	ITIMER_REAL,		&timer,
	NULL);
 
  } else /* no wait states */
    notify_set_itimer_func(frame, NOTIFY_FUNC_NULL,
	ITIMER_REAL,		NULL,
	NULL);
}

/*****************************************************************************/
void StopTraining()
/*** Stop button calls this to stop the run of the method ***/
{
  StopFlag=1;
}

/*****************************************************************************/
quit_proc()
/*** Quit button ***/
{
  static int QuitTimes=0;

  if (!QuitTimes) {
    if (xv_destroy_safe(frame) == XV_OK)
      exit(0);
    else
      ErrMesg(frame,"xv_destroy_safe did not return XV_OK\nTry again !!!!!");
  } else {
    xv_destroy(frame);
    exit(0);
  }
  QuitTimes++;
}

/*****************************************************************************/
show_norm_proc()
/*** Show button ***/
{
  char *fname;

  if (!strlen(fname = (char *)xv_get(fname_item,PANEL_VALUE))) return;

  Shrink=0;

  strcpy(Filename,fname);

  ShowImage(active_quad,Filename,Shrink);
}

/*****************************************************************************/
show_shrink_proc()
{
  char *fname;

  if (!strlen(fname = (char *)xv_get(fname_item,PANEL_VALUE))) return;

  Shrink=1;

  strcpy(Filename,fname);

  ShowImage(active_quad,Filename,Shrink);
}

/*****************************************************************************/
clear_sel_proc()
/*** Clear button ***/
{
  ClearQuadrant(active_quad);
}

/*****************************************************************************/
clear_graph_proc()
{
  XClearWindow(dpy,xid_graph_canvas);
  repaint_proc();
}

/*****************************************************************************/
clear_all_proc()
{
  XClearWindow(dpy,xid_canvas);
  XClearWindow(dpy,xid_graph_canvas);
  repaint_proc();
}

/*****************************************************************************/
select_active_quad_proc(item, value, event)
/*** Quad selector ***/
Panel_item item;
int value;
Event *event;
{
  active_quad=quads[value];
}

/*****************************************************************************/
start_hough() {
/*** Start method ***/
  /* change to explicit dispatching while computing */
  notify_stop();
  /* and set up a flag telling why notifier stopped */
  StoppedForTraining = 1;
  if (!popup_frame || xv_get(popup_frame, XV_SHOW)==FALSE)
    strcpy(Filename, (char *)xv_get(fname_item, PANEL_VALUE));
}

/*****************************************************************************/
pause_proc()
/*** Pause button ***/
{
  if (StoppedForTraining)
    notice_prompt(frame, NULL,
	NOTICE_MESSAGE_STRINGS,	"Program paused !", NULL,
	NOTICE_BUTTON_YES,	"Continue",
	NOTICE_NO_BEEPING,	TRUE,
	NULL);
  notify_dispatch();
}

/*****************************************************************************/
void change_colormap()
{
  if (inversed_graylevels) {
    fg_level = WhitePixel(dpy, DefaultScreen(dpy));
    pixel_table=norm_pixel_table;
    cms=norm_cms;
    gc=norm_gc;
    inversed_graylevels=FALSE;
  } else {
    fg_level = BlackPixel(dpy, DefaultScreen(dpy));
    pixel_table=inv_pixel_table;
    cms=inv_cms;
    gc=inv_gc;
    inversed_graylevels=TRUE;
  }
  xv_set(canvas, WIN_CMS, cms, NULL);
  xv_set(graph_canvas, WIN_CMS, cms, NULL);
  XSetForeground(dpy,gc,fg_level);
  repaint_proc();
}

/*****************************************************************************/
void change_grid_gl()
{
  grid_graylevel=grid_graylevel?0:255;
  repaint_proc();
}

/*****************************************************************************/
close_popup_proc()
/*** Close dialog box popup ***/
{
  xv_set(popup_frame, XV_SHOW, FALSE, NULL);
  xv_destroy_safe(popup_frame);
  popup_frame=(Frame)NULL;
}

/*****************************************************************************/
close_help()
/*** Make help dialog box invisible ***/
{
  xv_set(popup_help_frame, XV_SHOW, FALSE, NULL);
}

/*****************************************************************************/
cancel_popup_proc()
/*** Cancelling dialog box popup (not used) ***/
{
  hough_proc = no_active_method_proc;
  xv_set(message_item,  PANEL_LABEL_STRING, "  No active method", NULL);
  close_popup_proc();
}

/*****************************************************************************/
cycle_event_proc(item, event)
Panel_item item;
Event *event;
{
  if ((event_id(event)==ACTION_MENU) || (event_id(event)==MS_LEFT) &&
      event_is_up(event))
    if (xv_get(item, PANEL_VALUE))
      xv_set(item, PANEL_VALUE, FALSE, NULL);
    else
      xv_set(item, PANEL_VALUE, TRUE, NULL);
}

/*****************************************************************************/
void prepare_for_running()
{
  if (popup_frame && xv_get(popup_frame, XV_SHOW)==FALSE)
    close_popup_proc();

  init_pic(rpic,MAX_SIZE,MAX_SIZE);
  init_pic(gpic,MAX_SIZE,MAX_SIZE);
}

/*****************************************************************************/
void show_results()
{
  ClearQuadrant(quads[4]);
  PutPic(quads[4],gpic,dim_x,dim_y);
  ClearQuadrant(quads[5]);
  PutPic(quads[5],rpic,dim_x,dim_y);

/* Close dialog popup after run (someone likes - someone does not)
  if (popup_frame)
    close_popup_proc();
*/
}

/*****************************************************************************
	Procedure functions for reading new parameter values from the dialog box

 *****************************************************************************/
/*****************************************************************************/
read_common_parameters()
{
  strcpy(Filename,(char *)panel_get_value(File_item));
  strcpy(ParametersOut,(char *)panel_get_value(Param_item));
  Noise=atof((char *)panel_get_value(Noise_item));
/*
  OutputFormat=(int)panel_get_value(Format_item);
  PostScript=(int)panel_get_value(PS_item);
*/
  Shrink=(int)panel_get_value(Shrink_item);
  NumOfMaxs=atoi((char *)panel_get_value(NumOfMaxs_item));
  MinSegLen=atoi((char *)panel_get_value(MinSegLen_item));
  MaxSegWidth=atoi((char *)panel_get_value(MaxSegWidth_item));
  MaxSegGap=atoi((char *)panel_get_value(MaxSegGap_item));
}

/*****************************************************************************/
done_sht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    sht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
    sht_param.Dim_rho=atoi((char *)panel_get_value(Dim_rho_item));
    sht_param.Dim_theta=atoi((char *)panel_get_value(Dim_theta_item));

    if (sht_param.Dim_rho<3 || sht_param.Dim_rho>ACC_MAX_SIZE) {
      sht_param.Dim_rho=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
    if (sht_param.Dim_theta<3 || sht_param.Dim_theta>ACC_MAX_SIZE) {
      sht_param.Dim_theta=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item,  PANEL_LABEL_STRING, "  SHT for line: running", NULL);

    copy_pic(pic,pic1,dim_x,dim_y);

    PutPic(quads[0],pic,dim_x,dim_y);
    ClearQuadrant(quads[1]);

    printf("sht_line -f %s -m %d -l %d -w %d -g %d -t %d -n %f -d %d %d ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap, sht_param.Threshold,Noise,sht_param.Dim_theta, 
		sht_param.Dim_rho);

    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");

    sht_line(pic,pic1,gpic,rpic,dim_x,dim_y,sht_param.Dim_theta,sht_param.Dim_rho,
             real_params,NumOfMaxs,sht_param.Threshold,MinSegLen,MaxSegWidth,
             MaxSegGap,FALSE,ParametersOut);

    xv_set(message_item,  PANEL_LABEL_STRING, "  SHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
done_rht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    rht_param.MinDist=atoi((char *)panel_get_value(MinDist_item));
    rht_param.MaxDist=atoi((char *)panel_get_value(MaxDist_item));
    rht_param.Accuracy=atof((char *)panel_get_value(Accuracy_item));
    rht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item,  PANEL_LABEL_STRING, "  RHT for line: running", NULL);

    PutPic(quads[0],pic1,dim_x,dim_y);
    ClearQuadrant(quads[1]);

    printf("rht_line -f %s -m %d -l %d -w %d -g %d -n %f -t %d -d %d %d -a %f ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap,Noise, rht_param.Threshold, rht_param.MinDist, 
		rht_param.MaxDist,rht_param.Accuracy);

    if(Shrink) printf("-s ");
    printf("\n");

    rht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
             MaxSegWidth,MaxSegGap,rht_param.MinDist,rht_param.MaxDist,
            	rht_param.Accuracy,rht_param.Threshold,1,FALSE,FALSE,
		ParametersOut);

    xv_set(message_item,  PANEL_LABEL_STRING, "  RHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
done_drht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    drht_param.MinDist=atoi((char *)panel_get_value(MinDist_item));
    drht_param.MaxDist=atoi((char *)panel_get_value(MaxDist_item));
    drht_param.Accuracy=atof((char *)panel_get_value(Accuracy_item));
    drht_param.Accuracy1=atof((char *)panel_get_value(Accuracy1_item));
    drht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
    drht_param.Threshold1=atoi((char *)panel_get_value(Threshold1_item));
    drht_param.BlockWidth=atoi((char *)panel_get_value(BlockWidth_item));
    drht_param.MaxVar=atof((char *)panel_get_value(MaxVar_item));
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item,  PANEL_LABEL_STRING, "  DRHT for line: running", NULL);

    PutPic(quads[0],pic1,dim_x,dim_y);
    ClearQuadrant(quads[1]);
    ClearQuadrant(quads[2]);
    ClearQuadrant(quads[3]);

    printf("drht_line -f %s -m %d -l %d -w %d -g %d -n %f -t %d %d -d %d %d -a %f %f -v %f ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap,Noise, drht_param.Threshold,drht_param.Threshold1, 
		drht_param.MinDist, drht_param.MaxDist,drht_param.Accuracy,
		drht_param.Accuracy1, drht_param.MaxVar);

    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");

    drht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
              MaxSegWidth,MaxSegGap,drht_param.MinDist,drht_param.MaxDist,
              drht_param.Accuracy,drht_param.Accuracy1,drht_param.Threshold,
              drht_param.Threshold1,drht_param.BlockWidth,drht_param.MaxVar,
              1,FALSE,FALSE,ParametersOut);

    xv_set(message_item, PANEL_LABEL_STRING, "  DRHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
done_wrht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    wrht_param.Accuracy=atof((char *)panel_get_value(Accuracy_item));
    wrht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
    wrht_param.WinSize=atoi((char *)panel_get_value(WinSize_item));
    wrht_param.ThreshForPoints=atoi((char *)panel_get_value(ThreshForPoints_item));
    wrht_param.TolForFitting=atof((char *)panel_get_value(TolForFitting_item));
    wrht_param.DontRmAccu=(int)panel_get_value(DontRmAccu_item);
    wrht_param.ConnectiveCheck=(int)panel_get_value(ConnectiveCheck_item);
    wrht_param.Sectoring=(int)panel_get_value(Sectoring_item);
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item, PANEL_LABEL_STRING, "  WRHT for line: running", NULL);

    PutPic(quads[0],pic1,dim_x,dim_y);
    ClearQuadrant(quads[1]);

    printf("wrht_line -f %s -m %d -l %d -w %d -g %d -n %f -t %d -W %d -T %d -E %f -a %f ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap,Noise, wrht_param.Threshold, wrht_param.WinSize,
		wrht_param.ThreshForPoints,wrht_param.TolForFitting,
		wrht_param.Accuracy);
    if(wrht_param.DontRmAccu) printf("-D ");
    if(wrht_param.ConnectiveCheck) printf("-C ");
    if(wrht_param.Sectoring) printf("-S ");
    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");

    wrht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
             MaxSegWidth,MaxSegGap,wrht_param.Accuracy,wrht_param.Threshold,
             wrht_param.ThreshForPoints,wrht_param.TolForFitting,
             wrht_param.WinSize,wrht_param.DontRmAccu,wrht_param.ConnectiveCheck,
             wrht_param.Sectoring,1,FALSE,FALSE,ParametersOut);

    xv_set(message_item, PANEL_LABEL_STRING, "  WRHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
done_rwrht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    rwrht_param.MinDist=atoi((char *)panel_get_value(MinDist_item));
    rwrht_param.MaxDist=atoi((char *)panel_get_value(MaxDist_item));
    rwrht_param.Accuracy=atof((char *)panel_get_value(Accuracy_item));
    rwrht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
    rwrht_param.MinWinSize=atoi((char *)panel_get_value(MinWinSize_item));
    rwrht_param.MaxWinSize=atoi((char *)panel_get_value(MaxWinSize_item));
    rwrht_param.AccumulationsLimit=(int)panel_get_value(AccumulationsLimit_item);
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item,  PANEL_LABEL_STRING, "  RWRHT for line: running", NULL);

    PutPic(quads[0],pic1,dim_x,dim_y);
    ClearQuadrant(quads[1]);

    printf("rwrht_line -f %s -m %d -l %d -w %d -g %d -n %f -t %d -d %d %d -W %d %d -a %f ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap,Noise, rwrht_param.Threshold, rwrht_param.MinDist, 
		rwrht_param.MaxDist,rwrht_param.MinWinSize,
		rwrht_param.MaxWinSize, rwrht_param.Accuracy);
    if(rwrht_param.AccumulationsLimit) printf("-A ");
    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");

    rwrht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
               MaxSegWidth,MaxSegGap,rwrht_param.MinDist,rwrht_param.MaxDist,
               rwrht_param.Accuracy,rwrht_param.Threshold,
               rwrht_param.MinWinSize,rwrht_param.MaxWinSize,
               rwrht_param.AccumulationsLimit,1,FALSE,FALSE,
		ParametersOut);

    xv_set(message_item, PANEL_LABEL_STRING,
           "  RWRHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
done_cfht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    cfht_param.Msize=atoi((char *)panel_get_value(Msize_item));
    cfht_param.ThreshForPoints=atoi((char *)panel_get_value(ThreshForPoints_item));
    cfht_param.TolForFitting=atof((char *)panel_get_value(TolForFitting_item));
    cfht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
    cfht_param.Epsilon=atof((char *)panel_get_value(Epsilon_item));
    cfht_param.Gamma=atof((char *)panel_get_value(Gamma_item));
    cfht_param.DontRemovePoints=(int)panel_get_value(DontRemovePoints_item);
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item,  PANEL_LABEL_STRING, "  CFHT for line: running", NULL);

    PutPic(quads[0],pic1,dim_x,dim_y);
    ClearQuadrant(quads[1]);
    printf("cfht_line -f %s -m %d -l %d -w %d -g %d -t %d -n %f -W %d -T %d -E %f -e %f -A %f ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap, cfht_param.Threshold,Noise,cfht_param.Msize,
		cfht_param.ThreshForPoints,cfht_param.TolForFitting,
		cfht_param.Epsilon, cfht_param.Gamma);
    if(cfht_param.DontRemovePoints) printf("-D ");
    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");

    cfht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
              MaxSegWidth,MaxSegGap,cfht_param.Msize,
              cfht_param.ThreshForPoints,cfht_param.TolForFitting,
              cfht_param.Threshold,cfht_param.Epsilon,cfht_param.Gamma,
              cfht_param.DontRemovePoints,1,FALSE,FALSE,ParametersOut);

    xv_set(message_item, PANEL_LABEL_STRING, "  CFHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
done_probht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    probht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
    probht_param.Dim_rho=atoi((char *)panel_get_value(Dim_rho_item));
    probht_param.Dim_theta=atoi((char *)panel_get_value(Dim_theta_item));
    probht_param.SampleLevel=atof((char *)panel_get_value(SampleLevel_item));
    if (probht_param.Dim_rho<3 || probht_param.Dim_rho>ACC_MAX_SIZE) {
      probht_param.Dim_rho=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
    if (probht_param.Dim_theta<3 || probht_param.Dim_theta>ACC_MAX_SIZE) {
      probht_param.Dim_theta=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item, PANEL_LABEL_STRING, "  ProbHT for line: running", NULL);

    if (probht_param.SampleLevel>0.0 && probht_param.SampleLevel<100.0)
      take_a_sample_of_pic(pic,pic1,dim_x,dim_y,probht_param.SampleLevel);
    else
      copy_pic(pic1,pic,dim_x,dim_y);

    PutPic(quads[0],pic,dim_x,dim_y);
    ClearQuadrant(quads[1]);

    printf("probht_line -f %s -m %d -l %d -w %d -g %d -n %f -t %d -d %d %d -S %f",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap, Noise,probht_param.Threshold,probht_param.Dim_rho,
		probht_param.Dim_theta,probht_param.SampleLevel);

    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");

    sht_line(pic,pic1,gpic,rpic,dim_x,dim_y,probht_param.Dim_theta,
             probht_param.Dim_rho,real_params,NumOfMaxs,probht_param.Threshold,
             MinSegLen,MaxSegWidth,MaxSegGap,FALSE,ParametersOut);

    xv_set(message_item,  PANEL_LABEL_STRING,
	   "  ProbHT for line: done and waiting",
	   NULL);
  }
}

/*****************************************************************************/
done_aht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    aht_param.AccuRangeSel=atoi((char *)panel_get_value(AccuRangeSel_item));
    aht_param.Dim_m=atoi((char *)panel_get_value(Dim_m_item));
    aht_param.Dim_c=atoi((char *)panel_get_value(Dim_c_item));
    aht_param.BinarizeLevel=atof((char *)panel_get_value(BinarizeLevel_item));
    aht_param.Accur_m=atof((char *)panel_get_value(Accur_m_item));
    aht_param.Accur_c=atof((char *)panel_get_value(Accur_c_item));

    if (aht_param.Dim_m<3 || aht_param.Dim_m>ACC_MAX_SIZE) {
      aht_param.Dim_m=9;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
    if (aht_param.Dim_c<3 || aht_param.Dim_c>ACC_MAX_SIZE) {
      aht_param.Dim_c=9;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item,  PANEL_LABEL_STRING, "  AHT for line: running", NULL);

    copy_pic(pic,pic1,dim_x,dim_y);

    PutPic(quads[0],pic1,dim_x,dim_y);
    ClearQuadrant(quads[1]);
    printf("aht_line -f %s -m %d -l %d -w %d -g %d -n %f -d %d %d -S %d -L %f -a %f %f ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap, Noise,aht_param.Dim_m, 
		aht_param.Dim_c, aht_param.AccuRangeSel,aht_param.BinarizeLevel,
		aht_param.Accur_m,aht_param.Accur_c );

    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");

    aht_line(pic,pic1,gpic,rpic,dim_x,dim_y,aht_param.Dim_m,aht_param.Dim_c,
             real_params,NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,
             aht_param.AccuRangeSel,aht_param.BinarizeLevel,aht_param.Accur_m,
             aht_param.Accur_c,FALSE,ParametersOut);

    xv_set(message_item,  PANEL_LABEL_STRING, "  AHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
done_cht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    cht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
    cht_param.Segments=atoi((char *)panel_get_value(Segments_item));
    cht_param.Overlap=atoi((char *)panel_get_value(Overlap_item));
    cht_param.Dim_rho=atoi((char *)panel_get_value(Dim_rho_item));
    cht_param.Dim_theta=atoi((char *)panel_get_value(Dim_theta_item));
    cht_param.MaskSize=atoi((char *)panel_get_value(MaskSize_item));

    if (cht_param.Dim_rho<3 || cht_param.Dim_rho>ACC_MAX_SIZE) {
      cht_param.Dim_rho=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
    if (cht_param.Dim_theta<3 || cht_param.Dim_theta>ACC_MAX_SIZE) {
      cht_param.Dim_theta=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item,  PANEL_LABEL_STRING, "  CHT for line: running", NULL);

    copy_pic(pic,pic1,dim_x,dim_y);

    PutPic(quads[0],pic,dim_x,dim_y);
    ClearQuadrant(quads[1]);

    printf("cht_line -f %s -m %d -l %d -w %d -g %d -t %d -n %f -d %d %d -S %d -O %d -W %d ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap, cht_param.Threshold,Noise,cht_param.Dim_rho, 
		cht_param.Dim_theta,cht_param.Segments, cht_param.Overlap,
		cht_param.MaskSize);

    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");

    cht_line(pic,pic1,gpic,rpic,dim_x,dim_y,cht_param.Dim_theta,cht_param.Dim_rho,
             real_params,NumOfMaxs,cht_param.Threshold,cht_param.Segments,
             cht_param.Overlap,cht_param.MaskSize,MinSegLen,MaxSegWidth,
             MaxSegGap,FALSE,ParametersOut);

    xv_set(message_item,  PANEL_LABEL_STRING, "  CHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
done_dcht_line_proc()
{
  if (popup_frame) { /* read new values if dialog box is popped up */
    read_common_parameters();
    dcht_param.Threshold=atoi((char *)panel_get_value(Threshold_item));
    dcht_param.Dim_theta=atoi((char *)panel_get_value(Dim_theta_item));

    if (dcht_param.Dim_theta<3 || dcht_param.Dim_theta>ACC_MAX_SIZE) {
      dcht_param.Dim_theta=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator quantization\nUsing default value...");
    }
  }

  if (make_object_from_image(Filename,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    xv_set(message_item,  PANEL_LABEL_STRING, "  DCHT for line: running", NULL);

    PutPic(quads[0],pic1,dim_x,dim_y);
    ClearQuadrant(quads[1]);
    printf("dcht_line -f %s -m %d -l %d -w %d -g %d -t %d -n %f -d %d ",
		Filename, NumOfMaxs, MinSegLen, MaxSegWidth,
		MaxSegGap, dcht_param.Threshold,Noise,dcht_param.Dim_theta);

    if(Shrink) printf("-s ");
    if(strlen(ParametersOut)) printf("-P %s ",ParametersOut);
    printf("\n");


    dcht_line(pic,pic1,gpic,rpic,dim_x,dim_y,dcht_param.Dim_theta,real_params,
              NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,dcht_param.Threshold,
              1,FALSE,ParametersOut);

    xv_set(message_item,  PANEL_LABEL_STRING, "  DCHT for line: done and waiting",
           NULL);
  }
}

/*****************************************************************************/
no_active_method_proc()
{
  ErrMesg(frame,"Please, select first one of the methods !");
}

/*** Error message handling ***/
char err_msg_buff[1024];

/*****************************************************************************/
ErrMesg(frame,s)
Frame frame;
char *s; 
{ 
  int res;

  if (frame == NULL) {
    perror(s);
    return;
  }
  if  (errno > 0 && errno < sys_nerr) 
    sprintf(err_msg_buff,"%s: %s",s,sys_errlist[errno]);
  else
    strcpy(err_msg_buff,s);

  res = notice_prompt(frame, NULL,
	NOTICE_MESSAGE_STRINGS,	" A problem was encountered:", err_msg_buff ,NULL,
	NOTICE_BUTTON_YES,	"Press to continue",
	NULL);
}

/*****************************************************************************/
void err_prnt(s)
char *s;
{
  ErrMesg(frame,s);
}

/*****************************************************************************/
void out_prnt(s)
char *s;
{
  fprintf(stdout,"%s\n",s);
}

/*****************************************************************************/
void warn_prnt(s)
char *s;
{
  ErrMesg(frame,s);
}

/*****************************************************************************/
save_proc()
/* Read selections and save image */
{

  if (popup_frame) { /* read new values if dialog box is popped up */
   if(!dim_x || !dim_y) 
   	err_prnt("No image to save");
   else 
   {
    	strcpy(Filename2,(char *)panel_get_value(File_item));
	if(strlen(Filename2)==0) { 
		err_prnt("Missing filename. Using 'tmp' as filename");  
		strcpy(Filename2, "tmp");
	}	
    	OutputFormat=(int)panel_get_value(Format_item)+ 1;
    	SavePic=(int)panel_get_value(Save_item);
	Ask=(int)panel_get_value(Ask_item);
	if (!Ask) {
    	  if(!SavePic) 
   		pic_to_disk(Filename2,
			rpic,dim_x,dim_y,OutputFormat,proc_name);
    	  else
   		pic_to_disk(Filename2,
			gpic,dim_x,dim_y,OutputFormat,proc_name);
 	  } 
 	  else {	
    	  if(!SavePic) 
   		invert_pic_to_disk(Filename2,
			rpic,dim_x,dim_y,OutputFormat,proc_name);
    	  else
   		invert_pic_to_disk(Filename2,
			gpic,dim_x,dim_y,OutputFormat,proc_name);
	}
    	xv_set(message_item,  PANEL_LABEL_STRING, "  Image saved ",
           NULL);
  	OutputFormat--;
   }
 } 
}

/*****************************************************************************/
save_image()
/* Create popup frame for saving selections */
{
  int row=0;
  char *fname;
 
  if (popup_frame)
    close_popup_proc();
 
  popup_frame=xv_create(frame, FRAME_CMD,
        XV_LABEL,               "Save image",
        NULL);
   popup_panel=xv_get(popup_frame, FRAME_CMD_PANEL);
   
   File_item = ask_str_value(popup_panel,0,row++,50,"Output image name: ",
	Filename2);
   Format_item = ask_format_value(popup_panel,0,row++,"Output image format: ",
                                     OutputFormat);

   Save_item = ask_save_value(popup_panel,0,row++,"Save image: ",
                                     SavePic);

   Ask_item = ask_cycle_value(popup_panel,0,row++,
                                    "Invert colors in saved image: ",
                                     Ask);
  xv_create(popup_panel, PANEL_BUTTON,
	PANEL_NOTIFY_PROC,	save_proc,
	PANEL_LABEL_STRING,	"Save",
	NULL);

  xv_create(popup_panel, PANEL_BUTTON,
	PANEL_NOTIFY_PROC,	close_popup_proc,
	PANEL_LABEL_STRING,	"Cancel",
	NULL);

  window_fit(popup_panel);
  window_fit(popup_frame);
  xv_set(popup_frame, WIN_SHOW, TRUE, NULL);

}

#if 0
/*****************************************************************************/
get_save_params()
/* Read saving parameters */
{
    	strcpy(ParametersOut,(char *)panel_get_value(File_item));
}


/*****************************************************************************/
save_parameters()
/* Create popup frame for saving parameters */
{
  int row=0;
  char *fname;
 
  if (popup_frame)
    close_popup_proc();
 
  popup_frame=xv_create(frame, FRAME_CMD,
        XV_LABEL,               "Save parameters",
        NULL);


   popup_panel=xv_get(popup_frame, FRAME_CMD_PANEL);
   
   File_item = ask_str_value(popup_panel,0,row++,50,"Parameters save file: ",
	ParametersOut);

   xv_create(popup_panel, PANEL_MESSAGE,
	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT,
	PANEL_VALUE_DISPLAY_LENGTH, 4,
	PANEL_LABEL_STRING,	"File must be set before a selected method is run!     ",
	NULL);


  xv_create(popup_panel, PANEL_BUTTON,
	PANEL_NOTIFY_PROC,	get_save_params,
	PANEL_LABEL_STRING,	"Ok",
	NULL);

  window_fit(popup_panel);
  window_fit(popup_frame);
  xv_set(popup_frame, WIN_SHOW, TRUE, NULL);

}

#endif

/*****************************************************************************/
void show_usage()
{
  fprintf(stderr,"Usage: %s ",proc_name);
  /* show_common_usage(stderr); */
  fprintf(stderr,"[-f Edge_pic] [-I] [-G] [-h] [-u]\n");
}

/*****************************************************************************/
void print_options()
{
  fprintf(stdout,"XHoughtool - graphical user inteface for the Houghtool\n\n");
  fprintf(stdout,"Note to set the environment variable XHTOOLHELPDIR \n");
  fprintf(stdout,"to define the path to xhoughtool.help. The default path is \n");
  fprintf(stdout,"a directory where the XHoughtool is started. \n\n");
  fprintf(stdout,"Options: \n");
  fprintf(stdout,"\t-f Edge_pic\n");
  fprintf(stdout,"\t   The file name of a test edge image.\n");
  fprintf(stdout,"\t   The format of the image file can be PGM, CVL, SKE, \n"); 
  fprintf(stdout,"\t   VIS (Visilog), or RAW.\n");
  fprintf(stdout,"\t-I\n");
  fprintf(stdout,"\t   Inverse canvas colors.\n");
  fprintf(stdout,"\t-G \n");
  fprintf(stdout,"\t   Show grid.\n");
  fprintf(stdout,"\t-h \t\t\t\t\n");
  fprintf(stdout,"\t   Print this help.\n");
  fprintf(stdout,"\t-u \n");
  fprintf(stdout,"\t   Usage.\n");
}
