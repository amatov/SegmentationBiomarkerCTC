.TH xhoughtool l "4 December 1995"
.SH NAME
xhoughtool - a software package for Hough Transform Calculation.
.SH SYNOPSIS
.B xhoughtool
[
.BI \-f \ Edge_pic 
]
[
.BI \-I  
]
[
.BI \-G
]
[
.BI \-h 
]
[
.BI \-u 
]
 

.SH DESCRIPTION
The xhoughtool program opens a graphical user interface, the XHoughtool, 
for visualizing the performance of the Hough Transform algorithm and its variants for line detection. 
The user interface is X Window System based. 
The Houghtool package contains several algorithms that can also be used without the graphical interface 
through UNIX-like commands. 

.SH OPTIONS
.TP 15
.B -f Edge_pic
Specify the file name of a test edge image. The format of the image file can be PGM, CVL, SKE, VIS (Visilog), or RAW. If the name is not given here, it is asked in the package later.
.TP 15
.B  -I
Inverses canvas colors. This changes color to white since the default color is black.
.TP 15
.B -G
Show grid. 
.TP 15
.B -h 
Print the help of options.
.TP 15
.B -u
Print short usage information.
.TP 0
For example: 
.TP 35
.B xhoughtool -f house_edge.pgm -G 
/* load  house_edge.pgm and show grid */

.SH IMPLEMENTED ALGORITHMS
Ten Hough Transform (HT) algorithms are implemented in the Houghtool software package. 
Programs implemented are available both as separate command line commands (i.e. UNIX-like 
filter programs) and as menu buttons in the XHoughtool. The commands are as follows:
 
Non-probabilistic HTs:
.TP 15
.B  aht_line      
Adaptive HT (AHT)
.TP 15
.B  cht_line      
Combinatorial HT (CHT)
.TP 15
.B  cfht_line     
Curve Fitting HT (CFHT)
.TP 15
.B  sht_line      
Standard HT (SHT)
.TP 15
.B Probabilistic HTs:
.TP 15
.B  dcht_line     
Dynamic Combinatorial HT (DCHT)
.TP 15 
.B  probht_line   
Probabilistic HT by Kiryati et al. (ProbHT)
.TP 15
.B  rht_line      
Randomized HT (RHT)
.TP 15
.B  drht_line     
- Dynamic RHT (DRHT)
.TP 15 
.B  rwrht_line    
- Random Window RHT (RWRHT) 
.TP 15
.B  wrht_line     
- Window RHT (WRHT) 
.TP 15
.B 	
  including Connective RHT (CRHT) 
 
.TP 0
Programs provide more information about themselves with 
options -h and -u, e.g., 

rht_line -u  /* for the short description of usage */
.br
rht_line -h  /* for more details about the options */ 

Programs can be compiled by "make".

Makefile usage: make [programname|all|clean|help]
.TP 15
.B programname 
This is one of the names listed above, e.g. make rht_line. 
Only the files associated with the RHT algorithms are now compiled. 
Note that the command xhoughtool is also considered here as programname.
.TP 15
.B all
Compile all programs.
.TP 15
.B clean 
Delete core and all object files. 
.TP 15
.B help 
Print the usage information.

.TP 0
For further information about the algorithms see, e.g., 
Kalviainen, H., Hirvonen, P., Xu, L., and Oja, E.,
Probabilistic and Non-probabilistic Hough Transforms:
Overview and Comparisons, Image and Vision Computing, vol. 13, no. 4, May 1995, pp. 239-252.  

This article is available by  ftp://lut.fi/itlab/rht/reports/Kalviainen_IVC-95.ps.Z or 
http://www.lut.fi/dep/tite/XHoughtool/xhoughtool.html. 
The help window of the Xhoughtool also gives references to the methods. 

.SH BASIC SCREEN LAYOUT
The XHoughtool window is divided into three major areas:
.TP 4
.B \- 
The upper right part of the window contains a control panel that gives access
to hierarchical menu system. 
.TP 4
.B \-
The bottom right part of the window is a graphics subwindow composed of six
regions for showing images.
.TP 4
.B \-
The left-hand side part of the window contains three separate areas for
visualizing the accumulator space.

.SH BUTTONS AND MENUS
.SS Image file path
This is where you define directory where to find input edge pictures. One edge
picture can be defined here, ready for HT methods selected later.
.SS Help
Opens the help window.
.SS Non-probabil. and Probabilistic
With these menus you can select a method and set its parameters.
.SS Common options for all the methods
.TP 10 
.B Edge file
Defines the binary edge image that are used in the HT calculation. Default "Image file path".
.TP 10 
.B Parameters save file
This is where you define, where to save parameters,
which will be found during the next run. NOTE that it will overwrite the previous run
with the same name.
.TP 10 
.B Vertical shrink
Shrink images vertically by 1/3. Default off.
.TP 10 
.B Randomly added noise
To add noise to the source picture. Default no noise.
.TP 10 
.B Number of maxima 
Defines how many lines are searched at maximum.
.TP 10 
.B Minimum line segment length:
Lines shorter than defined here are left out of consideration.
.TP 10 
.B Line scanning width 
Pixels found in the width of the scanning area are considered to belong to a line.
.TP 10 
.B Maximum gap between pixels in line segment 
If the gap between pixels is larger than
defined here, pixels are not considered to belong to the same line.
.SS Non-probabil. Hough Transforms
.SS SHT
.TP 10 
.B Minimum score accepted for a accumulator max. 
When points are accumulated into
the accumulator space, this is the minimum number of hits in a cell to consider
a cell as a maximum which may represent line parameters.
.TP 10 
.B Accumulator rho quantization 
Accumulator array size (max quantization of (rho,theta) is 256 x 256). 
.TP 10 
.B Accumulator theta quantization 
Accumulator array size (max quantization of (rho,theta) is 256 x 256).
.SS CFHT
.TP 10 
.B Window size 
The size of the window where a line is fitted. The size is (2*Msize+1) x (2*Msize+1).
If Msize = 6 (default), there is a 13 x 13 window. 
.TP 10 
.B Threshold for amount of windowed points 
Minimum number of points in the window. If there are less points in the window,
the line fitting is not performed.
.TP 10 
.B Tolerance for fitting error 
Maximum fitting error accepted.
.TP 10 
.B Accumulator threshold 
Determines how many hits an accumulator cell must have to be considered as a maximum
whose parameters may represent a line.
.TP 10 
.B Tolerance to existing accumulator cells 
Tolerance of a new candidate cell to near accumulator cells when updating the accumulator.
A new candidate cell is merged into an old one,
if it satisfies the tolerance; otherwise, a new cell is created. The larger 
the tolerance, the less new cells are created. Tolerance is the Euclidean distance
between positions of cells.
.TP 10 
.B Weight for old score 
Value for weighting accumulator cell position. Old position is
weighted by (1- Gamma) and the new one by Gamma. If Gamma is 0.0, the old position is
weighted by it's score, and the new one by one.
.TP 10 
.B Do not remove edge points  
If set edge point are not removed from image after a found line. Default off. 
.SS AHT
.TP 10 
.B Accumulator range selection 
1 for horizontal lines (0-45 or 135-180 degrees), 2 for vertical lines (45-135 degrees).
.TP 10 
.B Accumulator slope quantization 
The slope quantization in the accumulator array (slope x intercept).
.TP 10 
.B Accumulator intercept quantization 
The intercept quantization in the accumulator array (slope x intercept).
.TP 10 
.B Accumulator binarizing level 
Level to which the accumulator is binarized with respect to the maximum value. If
level is 0.9, a cell whose value is greater or equal to 0.9 * maximum value 
gets value 1; otherwise 0.
.TP 10 
.B Relative slope accuracy 
The relative accuracy of the parameter slope. This is with respect to quantization, i.e.
the iteration of the accumulator is repeated until the difference between the smallest and
largest value of the slope is less than quantization * slope accuracy. 
.TP 10 
.B Relative intercept accuracy 
The relative accuracy of the parameter intercept. This is with respect to quantization, i.e.
the iteration of the accumulator is repeated until the difference between the smallest and
largest value of the intercept is less than quantization * intercept accuracy.
.SS CHT
Minimum score accepted for an accumulator max 
The threshold for the accumulator.
.TP 10 
.B Number of segments in image 
The number of segments in one dimension, i.e., the image is divided into Segments^2.
.TP 10 
.B Number of overlapped pixels in segments 
Defines how much segments are allowed to overlap.     
.TP 10 
.B Accumulator rho quantization 
Accumulator array size (max quantization of (rho,theta) is 256 x 256).
.TP 10 
.B Accumulator theta quantization 
Accumulator array size (max quantization of (rho,theta) is 256 x 256).
.TP 10 
.B Local maxima search mask size 
The size of the local maxima search mask. The real size is (2*MaskSize+1) x (2*MaskSize+1).
.SS Probabilistic Hough Transforms
.SS RHT
.TP 10 
.B Min. dist. for point pair 
Points selected at random must be at least this far from each other.
.TP 10 
.B Max. dist. for point pair 
Points selected at random can be at most this far from each other.
.TP 10 
.B Accumulator accuracy 
The resolution of the accumulator. The accumulator accuracy 0.01 means the accuracy of 
two decimals, e.g., 2.5623 = 2.56.
.TP 10 
.B Accumulator threshold 
Determines how many hits an accumulator cell must have to be considered as a maximum 
whose parameters may represent a line. 
.SS DRHT
.TP 10 
.B Min. dist. for point pair 
Points selected at random must be at least this far from each other.
.TP 10 
.B Max. dist. for point pair 
Points selected at random can be at most this far from each other.
.TP 10 
.B Accumulator accuracy 
The resolution of the accumulator. The accumulator accuracy 0.01 means the accuracy of 
two decimals, e.g., 2.5623 = 2.56.
.TP 10 
.B 2nd Accumulator accuracy 
The resolution of the accumulator. The accumulator accuracy 0.01 means the accuracy of
two decimals, e.g., 2.5623 = 2.56 (for second iteration).
.TP 10 
.B Accumulator threshold 
Determines how many hits an accumulator cell must have to be considered as a maximum
whose parameters may represent a line.  
.TP 10 
.B 2nd Accumulator threshold 
Determines how many hits an accumulator cell must have to be considered as a maximum
whose parameters may represent a line (for second iteration).
.TP 10 
.B Block width for 2nd iteration 
The width of the block (2*BlockWidth+1).
.TP 10 
.B Max. variation of "a" in degrees 
Maximum variation of a slope.
.SS WRHT (CRHT)
.TP 10 
.B Accumulator accuracy 
The resolution of the accumulator. The accumulator accuracy 0.01 means the accuracy of
two decimals, e.g., 2.5623 = 2.56.
.TP 10 
.B Accumulator threshold 
Determines how many hits an accumulator cell must have to be considered as a maximum
whose parameters may represent a line.
.TP 10 
.B Window size 
The size of the window where a line is fitted. The size is (2*size+1) x (2*size+1).
.TP 10 
.B Threshold for min. num. of windowed points 
Minimum number of points in the window. If the number of points is less, no fitting is performed. 
.TP 10 
.B Tolerance for fitting error 
Maximum fitting error accepted.
.TP 10 
.B Do not remove accumulator after line found 
If set, accumulator is not removed after line is found.
.TP 10 
.B Use only center point connected points in window (CRHT) 
If set, the method is CRHT. Only the points which are connected to the center point of the window
are used.
.TP 10 
.B Do the connective point search as sectored 
Do the connective component check for the window as sectored.
.SS RWRHT
.TP 10 
.B Min. dist. for point pair 
Points selected at random must be at least this far from each other.
.TP 10 
.B Max. dist. for point pair 
Points selected at random can be at most this far from each other.
.TP 10 
.B Accumulator accuracy 
The resolution of the accumulator. The accumulator accuracy 0.01 means the accuracy of
two decimals, e.g., 2.5623 = 2.56.
.TP 10 
.B Accumulator threshold 
Determines how many hits an accumulator cell must have to be considered as a maximum
whose parameters may represent a line.
.TP 10 
.B Min. window size 
Lower limit (MinWinSize*2+1) for the size of the window. The size of the window is selected at
random between the lower and upper limit.
.TP 10 
.B Max. window size 
Upper limit (MaxWinSize*2+1) for the size of the window. The size of the window is selected at
random between the lower and upper limit.
.TP 10 
.B Accumulate parameters in a win max 20T times 
Defines if accumulations from the same window are limited (on/off). Is on,
the maximum number of accumulations allowed is 20 * accumulator threshold.
.SS ProbHT
.TP 10 
.B Minimum accepted score for a accumulator max 
Determines how many hits an accumulator cell must have to be considered as a maximum
whose parameters may represent a line.
.TP 10 
.B Accumulator rho quantization 
Accumulator array size (max quantization of (rho,theta) is 256 x 256).
.TP 10 
.B Accumulator theta quantization 
Accumulator array size (max quantization of (rho,theta) is 256 x 256).
.TP 10 
.B Sampling level 
Defines the percentage of edge points which are selected to the HT calculation randomly.
20 % means that the fifth of all the edge points is chosen at random and
all the other points are simply ignored in the HT calculation.
.SS DCHT
.TP 10 
.B Accumulator threshold 
Determines how many hits an accumulator cell must have to be considered as a maximum
whose parameters may represent a line.
.TP 10 
.B Accumulator theta quantization 
The size of the histogram vector into where theta values are accumulated.
Maximum is 256 theta values in the vector.

.SS Start
Starts processing. The method has to be selected before this.
.SS Pause
Pauses processing until user pushes "Continue".
.SS Stop
Stops processing and shows results which are got so far.
.SS Wait states
Wait states lever is made for demonstration purposes. The more wait states you add,
the slower the process goes.
.SS Props
Menus to change the properties of the XHoughtool.
.TP 10
.B Background:white/black 
Select background color (black or white).
.TP 10
.B Grid: enable/disable    
Make the grid visible or invisible.
.TP 10 
.B Save image
.TP 10
.B 	Output image name <filename> 
Defines the name of an image to be saved, e.g. my_image.  
.TP 10 
.B 	Output image format:  (CVL, PGM, SKE, BIN)   
Defines format with what image will be saved and will append suffix accordingly., e.g. my_image.pgm.
.TP 10 
.B 	Save Image: (Result Image, Line Candidates)  
Defines which image to save. Choice is between the final result image and the line candidate image.
.TP 10 
.B 	Invert colors in saved image: (on, off) 
We decided to add this option because you might want 
to use images in your papers and not kill your laser printer by printing a lot of black color!
.SS Show
Options to show pictures.
.TP 10 
.B Show in normal way
Show images without vertical shrinking.
.TP 10 
.B Show pic as shrunken
Shows the image vertically shrunken by 1/3.
.SS Clear
Options to clear selected areas.
.TP 10 
.B Clear selected place
Clears selected place from canvas. Place is selected by selecting Place:
X.Y from interface.
.TP 10 
.B Clear graph canvas
Clears all six places on canvas.
.TP 10 
.B Clear All
Clears graphical canvas and accumulator representation areas (left hand side of canvas).
.SS Place
Indicates the current active area of a graphic subwindow and can be used to
select active place on canvas.
.SS Quit
Exits the XHoughtool and closes windows which are still open. 
This buton was named as "Done" in earlier versions.

.SH IMAGE FORMATS
Images used by HT algorithms must be binary edge pictures, having
a background signed by zeros (0) and edges by values between 1 and
255. Available formats are PGM, CVL, VIS (Visilog), SKE and RAW.
Two last ones have no header in an image file. In SKE format each
pixel is stored in space sized of INT (integer, 4 bytes). In RAW
format (like all the other formats) the pixel is stored in space
of byte.

Images should be 256x256 sized at maximum. To handle larger images,
you should make some changes to Makefile. 
See the file Makefile for further information.

.SH ENVIRONMENT
The filter programs can be compiled in Sun-4 machines. 
The graphical user interface the XHoughtool needs Xlib and
XView Toolkit libraries for being built up.


The help window of the XHoughtool works properly only if the help file
xhoughtool.help exists in the same directory where the XHoughtool
is started, or an environment variable XHTOOLHELPDIR is set
to define the directory where xhoughtool.help exists.

.SH FILES
.TP 20
.B xhoughtool.help 
Used as the help file of the XHoughtool. If XHTOOLHELPDIR not set, 
the file is tried to read from a current directory. 
.TP 20
.B Makefile 
The makefile for compiling Houghtool source programs. 
See IMPLEMENTED ALGORITHMS.

.SH SEE ALSO
X(1), make(1)

.SH BUGS
Since the graphical user interface needs Xlib and 
XView Toolkit libraries for being built up, unfortunately, 
you may find some problems while using Solaris, instead of SunOS. 
Binary codes should be compatible but problems, 
while compiling the source files, are possible.

However, notice that if you compile the codes of the filter programs or the XHoughtool in a SunOS machine but execute them in a Solaris machine, 
the proper behavior of library functions is not guaranteed. 
There is a bug in XHoughtool 1.0 due to the different behavior of 
a same library function between SunOS and Solaris, if codes are executed in Solaris environment. The seed of the random generator is not set properly.
This problem has been fixed in newer versions.

You can not use the source code directly, for example, in HP9000 machines 
without making necessary changes in Makefile. 

.SH COPYRIGHT
This program package is copyrighted in sense that it may be used 
for scientific purposes. THE PACKAGE as a whole, or parts thereof, 
CAN NOT BE INCLUDED OR USED IN ANY COMMERCIAL APPLICATION WITHOUT 
WRITTEN PERMISSION GRANTED BY ITS PRODUCERS. No programs contained 
in this package may be copied for commercial distribution. 
However, when used for scientific purposes, we wish that this 
package will be refereed in corresponding papers and reports.

Please refer to the following article when applying the Houghtool:
Kalviainen, H., Hirvonen, P., Oja, E., Houghtool-a Software
Package for Hough Transform Calculation, Proc. of the 9th
Scandinavian Conference on Image Analysis, Uppsala, Sweden,
June 1995, pp. 841-848.

The paper is available by 
ftp://lut.fi/itlab/rht/reports/Kalviainen_SCIA-95_htool.ps.Z, 
or http://www.lut.fi/dep/tite/XHoughtool/xhoughtool.html.

.SH TRADEMARKS
X Window System is a trademark of MIT. 
UNIX is a registered trademark of AT&T Bell Laboratories. 
Xview is a trademark of Sun Microsystems, Inc.

.SH AUTHORS
Petri Hirvonen, Jouni Ikonen, Pekka Kultanen, and Heikki Kalviainen, 
Lappeenranta University of Technology, Department of Information Technology, 
Laboratory of Information Processing, P.O. Box 20, FIN-53851 Lappeenranta, Finland. 
Telephone: +358-53-62111, Telefax:   +358-53-6212899.

Comments about       
.br
Houghtool software to:  Jouni.Ikonen@lut.fi
                        http://www.lut.fi/~jikonen
.br
RHT algorithms to:      Heikki.Kalviainen@.lut.fi
                        http://www/lut.fi/~kalviai>
.br
Other HT algorithms to: Original authors of the methods 
                        (See IMPLEMENTED ALGORITHMS)
