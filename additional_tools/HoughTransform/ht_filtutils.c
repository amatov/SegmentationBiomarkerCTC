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
 * File:    	ht_filtutils.c
 * Purpose: 	Houghtool filter program utilities
 * Date:    	Jun 1, 1993
 * Last change: Oct 10, 1995
 *****************************************************************************/

#include "ht_hough.h"
#include <strings.h>

extern char proc_name[256];

/*****************************************************************************/
void out_prnt(s)
/*
  Print string to <stdout>.

  Parameters :  s - string to be printed

*/
char *s;
{
  fprintf(stdout,"%s: %s",proc_name,s);
}

/*****************************************************************************/
void err_prnt(s)
/*
  Print error message string to <stderr>.

  Parameters :  s - string to be printed

*/
char *s;
{
  extern int errno;
  char *cptr=rindex(s,'\n'); /* pointing the last NEWLINE if there is one */

  fprintf(stderr,"%s: %s",proc_name,s);
  if (!cptr || strlen(s)>(int)(cptr-s+1))
    fprintf(stderr,"\n");
  if (errno)
    perror(s);
}

/*****************************************************************************/
void warn_prnt(s)
/*
  Print warning message string to <stderr>.

  Parameters :  s - string to be printed

*/
char *s;
{
  char *cptr=rindex(s,'\n'), /* pointing the last NEWLINE if there is one */
       bell='\7';

  fprintf(stderr,"%c%s: WARNING: %s",bell,proc_name,s);
  if (!cptr || strlen(s)>(int)(cptr-s+1))
    fprintf(stderr,"\n");
}

/*****************************************************************************/
void show_common_usage(stream)
/*
  Print usage to stream.

  Parameters :  stream - file stream

*/
FILE *stream;
{
  fprintf(stream,"[-f File] [-G GrayFile] [-R RealParamsFile] [-r ResultFile] ");
  fprintf(stream,"[-F Format] [-s] [-n Noise] [-l MinSegLen] [-w MaxSegWidth] ");
  fprintf(stream,"[-g MaxSegGap] [-m NumOfMaxs] [-M] [-i] [-P ResultParamsFile] [-h] [-u] ");
}

/*****************************************************************************/
void print_common_options(stream)
/*
  Print help on options to stream.

  Parameters :  stream - file stream

*/
FILE *stream;
{
  fprintf(stream,"Common options for all the Houghtool filters:\n");
  fprintf(stream,"\t-f File\n");
  fprintf(stream,"\t   The file name of the test edge image.\n\t   Allowed");
  fprintf(stream," image formats are: CVL, PGM raw, VIS, SKE and BIN.\n");
  fprintf(stream,"\t   If not given, the test image is read from <stdin>.\n");
  fprintf(stream,"\t   In that case, only the PGM raw format is supported.\n");
  fprintf(stream,"\t-G GrayFile\n");
  fprintf(stream,"\t   The file name of gray level image of the ");
  fprintf(stream,"edge image.\n\t   This image is used as");
  fprintf(stream," a background for the result image.\n");
  fprintf(stream,"\t-R RealParamsFile\n");
  fprintf(stream,"\t   The file name for a file containing real parameters ");
  fprintf(stream,"of the\n\t   curves in the input edge image.\n");
  fprintf(stream,"\t-r ResultFile\n");
  fprintf(stream,"\t   The file name for a file to store results.\n\t   ");
  fprintf(stream,"If not given, the result image will be written to <stdout>.\n");
  fprintf(stream,"\t-F FileFormat\t\t\t\t\tdef: 2 (PGM)\n");
  fprintf(stream,"\t   Result file format: CLV (1), PGM raw (2), SKE (3) or");
  fprintf(stream," BIN (4).\n");
  fprintf(stream,"\t-s\t\t\t\t\t\tdef: not in use\n");
  fprintf(stream,"\t   Shrink images vertically by 1/3.\n");
  fprintf(stream,"\t-n Noise\t\t\t\t\tdef: 0.0\n");
  fprintf(stream,"\t   Level of randomly added noise (%%).\n");
  fprintf(stream,"\t-l MinSegLen\t\t\t\t\tdef: 10\n");
  fprintf(stream,"\t   Minimum length of the curve segment.\n");
  fprintf(stream,"\t-w MaxSegWidth\t\t\t\t\tdef: 2\n");
  fprintf(stream,"\t   Line scanning width.\n");
  fprintf(stream,"\t-g MaxSegGap\t\t\t\t\tdef: 5\n");
  fprintf(stream,"\t   Maximum gap between pixels within a curve segment.\n");
  fprintf(stream,"\t-m NumOfMaxs\t\t\t\t\tdef: 999\n");
  fprintf(stream,"\t   The number of curve segments to extract.\n");
  fprintf(stream,"\t-M\t\t\t\t\t\tdef: not in use\n");
  fprintf(stream,"\t   Result file contains lines corresponding the ");
  fprintf(stream,"accumulator maximas.\n\t   If not set, the result ");
  fprintf(stream,"lines are segmented.\n");
  fprintf(stream,"\t-i\t\t\t\t\t\tdef: not in use\n");
  fprintf(stream,"\t   Print detection statistics.\n");
  fprintf(stream,"\t-P ResultParamsFile\t\t\t\tdef: not in use\n");
  fprintf(stream,"\t   The file name for a file to store found parameters. If suffix is\n");
  fprintf(stream,"\t   .rho, parameters are saved in rho theta form.\n");
  fprintf(stream,"\t-h \n");
  fprintf(stream,"\t   Print this help.\n");
  fprintf(stream,"\t-u \n");
  fprintf(stream,"\t   Usage.\n\n");

  fprintf(stream,"Options for this method:\n");
}
