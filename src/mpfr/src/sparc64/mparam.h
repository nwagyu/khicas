/* Various Thresholds of MPFR, not exported.  -*- mode: C -*-

Copyright 2005-2020 Free Software Foundation, Inc.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

/* Generated by MPFR's tuneup.c, 2018-02-22, gcc 7.3.0 */
/* gcc202.fsffrance.org (UltraSparc T5 (Niagara5)) with gmp snapshot
   gmp-6.1.99-20180221 (gmp-6.1.2 is broken on this machine),
   which defines -m64 -mptr64 -mcpu=niagara4 -Wa,-Av9d */

#define MPFR_MULHIGH_TAB  \
 -1,0,0,0,0,0,-1,0,0,0,-1,0,0,0,-1,0, \
 0,0,0,14,14,15,16,15,18,18,16,18,18,19,20,22, \
 22,23,22,23,26,27,26,27,26,27,27,27,27,26,26,27, \
 27,29,30,27,32,33,34,35,36,37,38,38,40,40,40,44, \
 44,44,46,46,44,44,46,46,44,44,46,46,44,53,52,54, \
 52,53,54,54,52,54,54,54,52,53,54,54,54,54,58,54, \
 58,58,58,64,72,72,70,72,72,72,72,72,76,76,72,76, \
 76,76,76,76,76,76,76,72,88,88,76,88,88,88,88,88, \
 88,92,92,92,88,92,88,88,88,92,92,92,88,92,92,96, \
 104,108,104,108,92,108,104,104,104,108,108,108,108,108,108,108, \
 108,108,104,108,108,108,108,108,108,106,106,108,108,108,108,108, \
 108,108,108,116,108,108,108,108,116,108,116,108,108,116,108,108, \
 116,108,108,108,116,112,116,116,116,136,144,144,136,148,136,144, \
 136,143,144,152,148,148,144,144,152,152,144,144,148,148,152,144, \
 151,152,144,144,152,144,152,144,150,152,152,152,172,152,148,176, \
 176,176,176,176,176,176,176,176,176,176,176,184,184,184,184,184, \
 184,184,184,176,184,176,176,176,184,184,184,184,184,184,184,184, \
 184,183,184,184,184,184,184,184,184,184,184,213,184,184,216,213, \
 213,184,213,219,213,207,213,213,213,212,213,213,213,216,216,216, \
 225,219,219,225,213,225,216,225,225,225,228,225,225,213,213,237, \
 237,225,225,237,237,261,261,237,237,261,261,259,261,261,225,261, \
 261,261,261,225,261,255,260,261,261,261,261,261,261,264,264,264, \
 261,273,264,273,273,273,273,273,273,273,272,273,273,273,273,273, \
 273,273,273,273,272,273,273,273,273,273,273,273,273,261,273,273, \
 273,273,273,273,273,264,264,264,271,273,273,267,273,273,273,273, \
 273,272,273,273,273,272,273,273,273,273,273,273,273,272,273,273, \
 273,276,273,273,272,273,273,273,273,273,273,273,273,273,273,273, \
 273,273,273,273,273,273,318,273,318,273,318,318,318,273,273,273, \
 318,318,318,318,318,318,318,318,318,318,318,318,318,317,318,318, \
 318,318,317,318,318,318,317,318,318,318,318,318,318,318,318,318, \
 342,342,318,342,318,342,342,342,342,318,318,342,342,342,342,342, \
 342,342,342,342,342,342,342,342,342,342,342,342,342,342,424,424, \
 424,424,424,424,424,424,424,424,424,400,408,408,424,408,400,424, \
 424,424,400,423,424,424,424,424,424,424,424,424,424,424,424,423, \
 424,424,423,424,424,424,424,424,424,424,424,424,424,424,424,424, \
 424,424,424,424,424,424,424,423,423,424,424,424,423,424,424,424, \
 423,424,424,424,424,424,424,424,424,424,424,424,424,424,424,424, \
 424,424,424,424,424,424,424,424,424,424,424,456,424,424,456,456, \
 456,456,456,456,456,456,456,456,456,456,456,456,456,456,456,456, \
 456,456,456,456,456,456,456,456,456,456,456,456,456,456,456,456, \
 456,455,456,456,455,456,456,456,456,456,456,455,456,455,456,456, \
 456,488,487,488,488,488,488,488,488,488,488,487,488,488,488,488, \
 456,456,488,488,488,488,488,552,552,456,552,552,552,560,552,568, \
 552,560,560,552,551,560,560,568,559,560,560,568,559,568,560,583, \
 567,567,551,584,600,568,568,599,600,600,600,567,568,584,600,584, \
 600,568,583,584,560,568,600,568,583,600,600,584,568,584,583,584, \
 584,584,584,584,584,584,584,584,600,600,600,600,583,584,599,584, \
 600,600,600,600,600,600,584,600,597,598,599,600,600,598,600,599, \
 600,598,599,600,600,600,600,584,600,599,600,599,600,599,615,600, \
 616,600,600,616,600,599,600,599,600,600,600,599,600,616,648,648, \
 632,648,648,648,648,648,680,656,680,680,680,680,680,680,680,680, \
 696,696,680,680,696,696,696,696,696,679,696,695,696,680,696,696, \
 696,678,696,696,696,696,696,680,680,680,688,688,680,680,680,712, \
 696,696,695,696,696,696,696,696,696,728,720,728,696,728,696,696, \
 727,728,727,728,712,728,728,696,696,696,728,720,728,728,727,728, \
 728,728,727,728,728,728,727,728,728,728,728,727,728,728,728,728, \
 728,728,727,727,728,728,728,728,728,728,728,727,728,728,728,728, \
 728,728,728,728,728,728,728,728,728,728,728,728,728,727,728,727, \
 728,728,728,728,726,727,728,728,728,727,728,728,726,727,728,728, \
 728,728,728,728,728,727,728,726,727,728,728,728,728,727,728,728, \
 727,712,728,728,720,726,727,728,728,728,728,728,728,727,727,824, \
 824,728,727,728,728,727,728,728,728,727,728,824,728,824,824,824, \
 824,856,808,856,824,824,824,824,824,856,855,856,855,856,856,856, \
 856,856,856,856,855,856,856,856,855,856,855,856,856,856,856,856 \

#define MPFR_SQRHIGH_TAB  \
 -1,0,0,0,0,0,-1,0,0,0,-1,7,8,8,10,10, \
 10,10,11,11,12,12,14,14,14,14,16,16,16,16,18,18, \
 18,18,20,20,20,20,22,22,22,24,24,26,26,24,26,26, \
 26,26,28,28,28,30,30,29,30,31,32,31,32,32,34,34, \
 34,34,36,36,36,36,38,38,38,38,40,39,40,40,42,42, \
 44,42,44,44,46,44,46,46,46,46,48,48,48,50,54,49, \
 52,50,52,52,54,54,58,56,58,56,58,58,56,58,58,60, \
 58,58,59,60,60,60,62,72,62,62,66,64,64,66,70,72, \
 72,72,72,72,72,72,70,72,72,72,72,76,72,76,76,76, \
 76,76,76,76,80,76,78,80,80,80,80,80,80,84,84,84, \
 84,84,88,84,84,84,88,88,88,88,88,88,88,88,92,92, \
 92,92,92,92,92,92,96,100,96,96,96,96,96,96,100,100, \
 100,100,104,104,104,104,104,104,104,104,104,112,120,104,112,119, \
 108,119,124,108,108,108,120,112,112,112,132,112,124,112,124,119, \
 120,119,118,119,117,117,119,117,124,120,130,120,120,120,126,125, \
 124,128,125,124,129,131,126,126,128,136,130,128,131,136,130,130, \
 131,131,137,131,136,134,137,136,137,136,138,142,136,138,137,137, \
 141,138,140,142,144,142,143,142,148,142,144,144,147,154,154,146, \
 147,149,149,154,160,148,156,160,154,153,152,153,152,154,155,155, \
 160,154,155,155,156,168,159,161,166,166,160,160,160,162,162,161, \
 162,166,166,172,173,166,167,166,172,166,167,167,178,174,178,172, \
 174,171,172,172,173,174,178,174,179,180,178,179,184,177,179,178, \
 178,180,182,180,180,185,190,184,184,184,190,213,184,213,190,213, \
 186,213,213,213,213,213,190,213,213,213,213,213,213,213,213,213, \
 213,213,213,213,216,213,213,213,225,213,213,213,225,225,225,213, \
 213,225,213,213,213,225,225,225,237,237,237,237,225,237,225,225, \
 225,237,237,237,237,237,228,225,237,237,237,237,237,237,237,237, \
 237,237,237,237,237,249,240,225,225,225,225,225,228,225,237,237, \
 234,237,237,237,237,237,237,237,237,237,240,233,237,245,237,237, \
 237,237,237,249,242,249,249,243,244,249,249,249,249,261,249,261, \
 249,248,256,249,249,259,249,261,249,249,249,249,249,252,252,255, \
 261,257,261,257,257,261,261,261,261,261,261,261,261,261,261,261, \
 261,261,261,261,261,261,273,273,269,273,268,273,273,273,273,273, \
 273,273,273,273,273,273,273,273,273,273,273,273,273,273,273,273, \
 278,285,285,285,285,285,285,285,285,285,285,285,285,285,288,288, \
 285,285,285,285,285,285,297,297,297,297,296,297,297,297,297,297, \
 297,297,297,292,297,293,295,296,297,297,297,297,297,297,303,297, \
 321,306,300,300,304,302,303,305,308,309,307,309,309,309,307,309, \
 309,309,309,321,309,321,319,319,321,321,319,321,321,321,321,321, \
 321,321,321,316,321,321,321,320,321,321,321,320,321,321,321,321, \
 329,333,333,328,333,333,333,332,333,333,333,333,333,333,333,333, \
 333,333,342,333,333,333,333,333,342,342,342,354,342,342,342,342, \
 342,341,342,342,342,342,342,342,342,342,348,354,354,348,348,426, \
 354,354,426,354,426,353,354,354,353,354,354,354,354,354,354,354, \
 354,354,450,365,365,450,449,366,365,366,366,366,450,425,426,450, \
 426,426,426,425,366,426,426,426,426,424,450,426,426,426,426,426, \
 425,426,450,424,425,426,426,450,449,450,450,449,450,450,450,449, \
 450,449,449,450,449,450,449,450,449,450,450,450,448,450,426,449, \
 425,450,426,450,425,426,426,426,426,426,426,426,425,426,426,426, \
 426,426,426,426,425,426,450,448,449,450,450,449,450,450,450,450, \
 450,450,450,450,450,450,450,450,450,450,450,450,450,450,450,450, \
 450,449,450,450,450,450,449,450,448,449,450,450,449,448,450,450, \
 448,450,450,450,450,450,449,450,450,448,448,450,450,448,449,450, \
 450,449,450,449,450,450,450,449,450,450,450,449,450,450,450,450, \
 450,450,450,450,450,450,450,450,450,450,450,450,449,450,450,450, \
 450,450,450,450,450,450,450,450,449,450,450,450,450,450,450,450, \
 488,450,488,488,488,487,488,488,568,568,560,560,599,600,567,568, \
 568,600,568,568,568,568,568,568,568,568,568,568,567,568,567,568, \
 568,568,567,568,568,567,568,568,568,568,568,568,567,568,599,567, \
 600,568,599,600,568,600,600,600,568,600,600,600,599,599,600,600, \
 600,599,568,600,600,600,600,600,599,600,600,600,568,568,600,598, \
 600,600,568,600,600,600,600,599,600,599,600,598,600,599,600,600, \
 600,599,600,599,600,600,600,600,600,600,600,600,600,599,600,600, \
 600,598,600,600,632,599,632,631,632,632,631,632,631,632,632,632 \

#define MPFR_DIVHIGH_TAB  \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /*0-15*/ \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /*16-31*/ \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /*32-47*/ \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /*48-63*/ \
 0,0,0,0,37,0,37,37,44,41,43,44,44,45,45,44, /*64-79*/ \
 44,45,45,45,45,51,51,52,51,53,52,49,49,52,50,50, /*80-95*/ \
 51,50,51,52,53,52,53,56,58,54,57,59,59,66,64,66, /*96-111*/ \
 64,75,74,68,72,72,72,72,72,72,74,72,72,68,72,67, /*112-127*/ \
 75,72,80,88,72,88,88,88,86,88,88,72,88,88,88,88, /*128-143*/ \
 88,84,88,88,88,88,88,88,84,90,88,88,88,84,88,87, /*144-159*/ \
 88,91,92,91,88,89,91,88,90,88,92,101,90,100,98,104, /*160-175*/ \
 104,104,102,104,104,104,104,106,104,104,104,104,104,104,100,100, /*176-191*/ \
 106,108,107,105,102,106,108,104,104,108,104,104,116,104,120,107, /*192-207*/ \
 112,116,112,114,114,116,116,116,116,116,128,136,116,144,116,128, /*208-223*/ \
 124,132,128,136,132,136,144,148,136,144,144,144,152,144,144,144, /*224-239*/ \
 148,148,144,148,144,144,144,144,144,148,160,132,135,136,136,136, /*240-255*/ \
 144,136,144,176,144,144,148,148,152,152,152,144,148,152,176,148, /*256-271*/ \
 176,144,144,140,150,176,152,176,176,176,148,144,148,144,176,152, /*272-287*/ \
 176,169,176,176,176,176,176,180,180,176,175,176,176,176,176,176, /*288-303*/ \
 180,176,175,176,182,176,176,176,175,176,176,175,176,176,180,176, /*304-319*/ \
 180,183,184,176,176,176,176,176,176,175,176,176,176,176,176,176, /*320-335*/ \
 176,182,200,180,176,176,186,180,180,176,180,200,184,180,184,188, /*336-351*/ \
 180,180,200,208,204,180,192,212,196,208,216,184,208,208,208,213, /*352-367*/ \
 208,208,208,207,208,215,208,208,212,208,212,215,208,212,210,208, /*368-383*/ \
 208,208,210,207,208,208,216,207,208,208,207,208,208,208,208,216, /*384-399*/ \
 215,215,215,215,214,216,216,208,208,216,216,232,232,214,224,216, /*400-415*/ \
 230,216,232,224,231,216,232,224,232,232,232,232,232,230,232,232, /*416-431*/ \
 232,240,232,231,232,232,240,228,231,232,232,232,232,232,231,256, /*432-447*/ \
 232,264,232,232,232,272,256,232,256,264,232,272,272,288,272,288, /*448-463*/ \
 288,280,272,272,288,271,256,296,288,288,288,288,296,288,296,288, /*464-479*/ \
 296,296,288,304,296,296,288,303,288,304,296,288,288,288,296,296, /*480-495*/ \
 296,304,296,264,300,288,288,303,288,296,304,304,296,288,288,304, /*496-511*/ \
 288,300,288,296,296,296,336,352,352,296,288,288,304,288,288,296, /*512-527*/ \
 303,272,296,287,352,288,296,352,352,352,304,288,352,296,352,352, /*528-543*/ \
 352,304,288,360,288,288,352,296,352,352,296,352,288,288,296,296, /*544-559*/ \
 288,288,288,304,296,304,304,352,368,304,304,352,352,296,312,320, /*560-575*/ \
 304,304,336,351,360,352,304,352,304,352,351,352,352,352,351,352, /*576-591*/ \
 352,352,352,352,352,350,352,352,352,352,368,352,360,360,352,352, /*592-607*/ \
 352,352,351,352,352,352,352,352,360,352,352,352,352,367,352,352, /*608-623*/ \
 352,352,368,352,356,368,360,352,352,360,352,368,368,352,368,368, /*624-639*/ \
 360,352,368,360,352,352,352,350,368,352,352,368,352,352,352,364, /*640-655*/ \
 359,352,351,352,360,352,352,352,352,352,352,352,352,352,359,351, /*656-671*/ \
 352,352,368,352,352,352,352,352,352,350,352,352,352,349,351,352, /*672-687*/ \
 352,352,360,360,352,360,360,368,352,364,368,352,360,368,362,368, /*688-703*/ \
 360,360,366,367,360,360,360,359,367,368,368,367,368,368,361,367, /*704-719*/ \
 368,368,367,368,368,368,365,368,367,368,367,367,368,368,416,432, /*720-735*/ \
 373,432,432,399,424,424,432,432,414,432,414,424,432,414,424,416, /*736-751*/ \
 400,424,424,424,424,424,423,426,424,432,422,423,424,416,426,432, /*752-767*/ \
 424,426,430,432,424,424,426,414,424,432,414,432,416,416,416,417, /*768-783*/ \
 414,414,416,416,416,414,416,432,426,432,432,415,432,414,432,424, /*784-799*/ \
 432,424,418,414,424,424,425,414,425,416,412,414,416,424,422,430, /*800-815*/ \
 432,424,424,430,432,432,424,430,416,423,416,424,432,425,432,426, /*816-831*/ \
 432,426,432,432,432,424,428,432,431,431,432,432,432,448,432,426, /*832-847*/ \
 448,432,432,432,456,456,462,462,448,464,448,462,432,464,461,456, /*848-863*/ \
 456,464,462,462,456,461,462,462,462,462,463,464,462,462,460,464, /*864-879*/ \
 528,456,464,462,460,461,462,460,464,464,464,528,544,544,464,528, /*880-895*/ \
 528,462,528,528,512,544,544,528,528,544,544,544,544,528,542,544, /*896-911*/ \
 544,544,544,544,544,528,512,544,512,542,545,528,544,542,528,544, /*912-927*/ \
 544,544,543,528,544,541,542,542,528,544,544,543,545,544,544,544, /*928-943*/ \
 544,545,542,544,544,546,542,544,544,544,544,546,542,544,544,544, /*944-959*/ \
 542,544,544,543,544,543,542,544,576,544,544,544,576,576,544,543, /*960-975*/ \
 544,543,544,544,544,544,576,544,545,576,576,541,592,544,544,576, /*976-991*/ \
 576,575,576,544,543,544,576,527,576,528,608,576,544,608,544,544, /*992-1007*/ \
 542,544,544,544,544,539,608,608,543,542,544,608,543,544,528,544 /*1008-1023*/ \

#define MPFR_MUL_THRESHOLD 13 /* limbs */
#define MPFR_SQR_THRESHOLD 13 /* limbs */
#define MPFR_DIV_THRESHOLD 4 /* limbs */
#define MPFR_EXP_2_THRESHOLD 1712 /* bits */
#define MPFR_EXP_THRESHOLD 3213 /* bits */
#define MPFR_SINCOS_THRESHOLD 21539 /* bits */
#define MPFR_AI_THRESHOLD1 -6778 /* threshold for negative input of mpfr_ai */
#define MPFR_AI_THRESHOLD2 549
#define MPFR_AI_THRESHOLD3 8223
/* Tuneup completed successfully, took 686 seconds */