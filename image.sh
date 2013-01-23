#!/bin/bash

rm -Rf images/t.*

invert vis=images/model.vis map=images/t.mp beam=images/t.bm cell=0.2 imsize=256 robust=0.5 options=systemp,mfs
clean map=images/t.mp beam=images/t.bm out=images/t.cl cutoff=2e-4 niters=1000
restor model=images/t.cl map=images/t.mp beam=images/t.bm out=images/t.cm
cgdisp in=images/t.cm,images/t.cm labtyp=arcsec device=/xs slev=a,3.411e-4 levs1=3,5,7 type=con,pix region=arcsec,box'(-5,-5,5,5)' options=mirr,beambl,wedge,full
