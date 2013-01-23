#!/bin/bash

rm -Rf resid.vis t.*

uvmodel vis=hd61005/hd61005_all.vis model=images/model.mp options=subtract out=resid.vis
invert vis=resid.vis map=t.mp beam=t.bm cell=0.2 imsize=256 robust=0.5 options=systemp,mfs
clean map=t.mp beam=t.bm out=t.cl cutoff=2e-4 niters=1000
restor model=t.cl map=t.mp beam=t.bm out=t.cm
cgdisp in=t.cm,t.cm labtyp=arcsec device=/xs slev=a,3.411e-4 region=arcsec,box'(-5,-5,5,5)' levs1=2,3,4,5,6,7 type=con,pix options=mirr,beambl,wedge,full
