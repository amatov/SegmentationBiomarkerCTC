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
 * File:    graphmacros.h
 * Purpose: header include for graphical macros
 * Date:    Jun 1, 1993
 *****************************************************************************/

#ifdef VISUAL_PACK /* Graphical macros definition, original functions defined
                      in graphics.c */

#define PutPic_MACRO(qd,pic,dx,dy)		PutPic(qd,pic,dx,dy)
#define PutPic_gray_MACRO(qd,pic,dx,dy)		PutPic_gray(qd,pic,dx,dy)
#define PutAccu_gray_MACRO(qd,pic,dx,dy)	PutAccu_gray(qd,pic,dx,dy)
#define PutPixel_MACRO(qd,x,y)			PutPixel(qd,x,y)
#define PutPixel_gray_MACRO(qd,x,y,val)		PutPixel_gray(qd,x,y,val)
#define PutPixel_graph_MACRO(gqd,x,y,val)	PutPixel_graph(gqd,x,y,val)
#define PutVector_MACRO(qd,x1,y1,x2,y2)		PutVector(qd,x1,y1,x2,y2)
#define PutVector_graph_MACRO(gqd,x1,y1,x2,y2)	PutVector_graph(gqd,x1,y1,x2,y2)
#define PutText_MACRO(qd,x,y,str)		PutText(qd,x,y,str)
#define PutLongText_MACRO(qd,x,y,val,cnd)	PutLongText(qd,x,y,val,cnd)
#define PutText_graph_MACRO(gqd,x,y,str)	PutText_graph(gqd,x,y,str)
#define PutText_graph_small_MACRO(gqd,x,y,str)	PutText_graph_small(gqd,x,y,str)
#define ClearQuadrant_MACRO(qd)			ClearQuadrant(qd)
#define ClearQuadrant_graph_MACRO(gqd)		ClearQuadrant_graph(gqd)
#define ShowImage_MACRO(qd,name,shrink)		ShowImage(qd,name,shrink)
#define DrawLines_MACRO(qd,start)		DrawLines(qd,start)
#define DrawCircles_MACRO(qd,cstart)		DrawCircles(qd,cstart)
#define DrawAccumulatorSpace_graph_MACRO(gqd,acc,h,xl,yl,xmi,xma,ymi,yma,accma) DrawAccumulatorSpace_graph(gqd,acc,h,xl,yl,xmi,xma,ymi,yma,accma)
#define DrawLineData_graph_MACRO(gqd,x,y,n,h,xl,yl) DrawLineData_graph(gqd,x,y,n,h,xl,yl)
#define DrawAccumulatorSpace_aht_MACRO(qd,acc,dm,dc,pic,dx,dy,sm,bm,sc,bc,lm,hm,lc,hc) DrawAccumulatorSpace_aht(qd,acc,dm,dc,pic,dx,dy,sm,bm,sc,bc,lm,hm,lc,hc)
#define FlushDisp_MACRO()			FlushDisp()
#define NotifyDispatch_MACRO()			notify_dispatch()

#else /* Null macros for nonvisual use of HT algorithms */

#define PutPic_MACRO(qd,pic,dx,dy)
#define PutPic_gray_MACRO(qd,pic,dx,dy)
#define PutAccu_gray_MACRO(qd,pic,dx,dy)
#define PutPixel_graph_MACRO(gqd,x,y,val)
#define PutPixel_MACRO(qd,x,y)
#define PutPixel_gray_MACRO(qd,x,y,val)
#define PutVector_MACRO(qd,x1,y1,x2,y2)
#define PutVector_graph_MACRO(gqd,x1,y1,x2,y2)
#define PutText_MACRO(qd,x,y,str)
#define PutLongText_MACRO(qd,x,y,val,cnd)
#define PutText_graph_MACRO(gqd,x,y,str)
#define PutText_graph_small_MACRO(gqd,x,y,str)
#define ClearQuadrant_MACRO(qd)
#define ClearQuadrant_graph_MACRO(gqd)
#define ShowImage_MACRO(qd,name,shrink)
#define DrawLines_MACRO(qd,start)
#define DrawCircles_MACRO(qd,cstart)
#define DrawAccumulatorSpace_graph_MACRO(gqd,acc,h,xl,yl,xmi,xma,ymi,yma,accma)
#define DrawLineData_graph_MACRO(gqd,x,y,n,h,xl,yl)
#define DrawAccumulatorSpace_aht_MACRO(qd,acc,dm,dc,pic,dx,dy,sm,bm,sc,bc,lm,hm,lc,hc)
#define FlushDisp_MACRO()
#define NotifyDispatch_MACRO()

#endif
