/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file bulge.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Azzam Haidar
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/

/***************************************************************************//**
 *  bulge chasing global definition for all L/U HE/HB/GE matrices.
 **/
#ifndef _MORSE_BULGE_H_
#define _MORSE_BULGE_H_

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  internal common routines to all bulgechasing function
 **/
inline static void findVTpos(int N, int NB, int Vblksiz, int sweep, int st, int *Vpos, int *TAUpos, int *Tpos, int *myblkid);
inline static void findVTsiz(int N, int NB, int Vblksiz, int *blkcnt, int *LDV);
inline static int morse_ceildiv(int a, int b);

////////////////////////////////////////////////////////////////////////////////////////////////////
inline static int morse_ceildiv(int a, int b)
{
  double r = (double)a/(double)b;
  r = (r-(int)r)==0? (int)r:(int)r+1;
  return (int) r;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
inline static void findVTpos(int N, int NB, int Vblksiz, int sweep, int st, int *Vpos, int *TAUpos, int *Tpos, int *myblkid)
{
  int prevcolblknb, prevblkcnt, prevcolblkid;
  int curcolblknb, nbprevcolblk, mastersweep;
  int blkid, locj, LDV;

  prevcolblknb = 0;
  prevblkcnt   = 0;
  curcolblknb  = 0;

  nbprevcolblk = sweep/Vblksiz;
  for (prevcolblkid = 0; prevcolblkid < nbprevcolblk; prevcolblkid++)
  {
       mastersweep  = prevcolblkid * Vblksiz;
       prevcolblknb = morse_ceildiv((N-(mastersweep+2)),NB);
       prevblkcnt   = prevblkcnt + prevcolblknb;
  }
  curcolblknb = morse_ceildiv((st-sweep),NB);
  blkid       = prevblkcnt + curcolblknb -1;
  locj        = sweep%Vblksiz;
  LDV         = NB + Vblksiz -1;

  *myblkid= blkid;
  *Vpos   = blkid*Vblksiz*LDV  + locj*LDV + locj;
  *TAUpos = blkid*Vblksiz + locj;
  *Tpos   = blkid*Vblksiz*Vblksiz + locj*Vblksiz + locj;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
inline static void findVTsiz(int N, int NB, int Vblksiz, int *blkcnt, int *LDV)
{
  int colblk, nbcolblk;
  int curcolblknb, mastersweep;

  *blkcnt   = 0;
  nbcolblk = morse_ceildiv((N-1),Vblksiz);
  for (colblk = 0; colblk<nbcolblk; colblk++)
  {
        mastersweep = colblk * Vblksiz;
        curcolblknb = morse_ceildiv((N-(mastersweep+2)),NB);
        *blkcnt      = *blkcnt + curcolblknb;
  }
  *blkcnt = *blkcnt +1;
  *LDV= NB+Vblksiz-1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

#endif
