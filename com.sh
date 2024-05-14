gfortran wood_anderson.f90 pgplotl.for BPFILTER.f  -o wood_anderson
./wood_anderson
gmt psconvert  wood_anderson.ps -A -Tg 