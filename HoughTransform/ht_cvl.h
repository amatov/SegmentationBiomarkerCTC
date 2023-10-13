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
 * File:    ht_cvl.h
 * Purpose: header include for read/write images having cvl format
 * Date:    Jun 1, 1993
 *****************************************************************************/
/*
 * CVL image header include file.
 * Entries consist of 2 byte integer short words.
 * A variable length comments section may follow this header.
 * 0 is an unknown length.
 */
struct cvl_p {
	long cv_tag;    /* 044520 'PI' Identifier for pictures */
                        /* 0152103 'CT'+0100000 continued with parity bit */
	short cv_lcom;  /* Length of comments (default 0) */
	short cv_type;  /* Type of pixel encoding (default 0) */
	short cv_dims;  /* dimensions (default 2) */
	short cv_hpls;  /* hyperplanes (default 1) */
	short cv_plns;  /* planes (default 1) */
	short cv_rows;  /* rows */
	short cv_cols;  /* columns */
	short cv_bnds;  /* bands (default 1) */
	short cv_bp;    /* bits per pixel (default 8) */
/* The following fields are optional (all default 0) */
/* They will be copied without interpretation by most software */
	short cv_ebb;   /* exponent bits per band (0 or 7) */
	short cv_sbb;   /* significant bits per band (no sign bit)*/
	short cv_bb;    /* bits per band (often redundant) */
	short cv_sbp;   /* significant bits per pixel (often 8 or 6) */
	short cv_ebp;   /* exponent bits per pixel (0 or 7) */
};
/*
 * Immediately following the header structure in
 * oldest-to-newest-order, ASCII title and comments
 * consist of lines of less than 80 characters
 * each ended by a newline (012).
 * New comments would be added after the last comment
 * and the comment length field of the header would be incremented.
 */
/*
 * The image data which follows is packed
 * with the least significant bit representing
 * the low number pixel.
 * Bits per pixel which are a power of 2 are encouraged.
 * No alignment is done at the end of planes or hyperplanes,
 * however each row begins on a byte boundary.
 * Images are assumed to raster across the top row
 * from left to right and then downward by rows.
 */
