#!/bin/sh

pdftk bachelor-thesis.pdf cat 1-6 output tmp1.pdf
pdftk bachelor-thesis.pdf cat 7-end output tmp2.pdf
convert -density 600 +antialias tmp1.pdf tmp3.pdf
pdftk tmp3.pdf tmp2.pdf output thesis-zban.pdf
