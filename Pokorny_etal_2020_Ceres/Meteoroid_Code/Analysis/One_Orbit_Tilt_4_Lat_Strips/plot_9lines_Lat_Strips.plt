ifile=1
CANX=1050
CANY=1500
res=2

#scale = 3*3 # This is done on 6 degree grid so we are missing 3 x 3 
scale = 1
scale = scale/(6378000.)**2.0/(1.0*pi) # Cross-section of Earth is seen by each meteoroid -> converts to kg per m^2 per day
scale = scale/86400.0 # converts to kg per second
CRAT_DIA=63.0
#scale = scale * 0.1 # from kg/m^2 to g/cm^2
set encoding iso_8859_1
load "~/viridis.pal"
  mwhigh(x) = 1.0 + -0.919061*(abs(1.0-x))**0.674635
  mwmed(x) = mwhigh(x) -0.0807765 + 0.161136*x -0.0796311*x**2
  mwlow(x) = mwmed(x) -3.53551e-05 + 0.000645749*x
  mwst(x) = x < 0.2 ? mwlow(x) : (x < 0.9 ? mwmed(x) : mwhigh(x) )
  mwt(b) = b>0 ? asin(mwst(sin(abs(b)))) : -asin(mwst(sin(abs(b))))
  mwx(b,l) = 2.0*sqrt(2.0) * l * cos(mwt(b))*(180/2/pi/sqrt(2))
  mwy(b,l) = sqrt(2.0) * sin(mwt(b))*(90/sqrt(2))


  bconv(b) = b*pi/180
  lconv(l) = (l<180) ? l*pi/180 : (l-360)*pi/180



set term pngcairo enhanced font "Times,40" size CANX*res,CANY*res lw 3

set xr[0:500]
set xtics 0,30,360 rotate by 30 right
set obj 1 rectangle from 360,0 to 500,graph 1 fs solid fc rgb "#eeeeee" lw 0 back
set key center right
set out "Line_Plot_9areas_Lat_Strip.png"
set multiplot layout 3,1
set grid
set xl "True Anomaly Angle (deg)"
set format x "%.0f\260"

set key Left
set yl "Mass flux (kg m^{-2} s^{-1})"

set obj 601 circle at graph 0.035, graph 0.915 size graph 0.03 fc rgb "gray" fs solid 0.55 border lw -1
set label 11 "A" at graph 0.035, graph 0.92  font "Times,56" center front

set key samplen 1
p "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_1.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_2.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_3.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_4.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_5.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_6.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_7.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_8.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_9.txt | sort -k1 -g" u 1:($4*scale):($6*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_1.txt | sort -k1 -g" u 1:($2*scale) w l ls 1 lw 2  t "-90{\260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-70{/E \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_2.txt | sort -k1 -g" u 1:($2*scale) w l ls 2 lw 2  t "-70{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-50{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_3.txt | sort -k1 -g" u 1:($2*scale) w l ls 3 lw 2  t "-50{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-30{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_4.txt | sort -k1 -g" u 1:($2*scale) w l ls 4 lw 2  t "-30{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-10{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_5.txt | sort -k1 -g" u 1:($2*scale) w l ls 5 lw 2  t "-10{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}10{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_6.txt | sort -k1 -g" u 1:($2*scale) w l ls 6 lw 2  t "10{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}30{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_7.txt | sort -k1 -g" u 1:($2*scale) w l ls 7 lw 2  t "30{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}50{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_8.txt | sort -k1 -g" u 1:($2*scale) w l ls 8 lw 2  t "50{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}70{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_9.txt | sort -k1 -g" u 1:($2*scale) w l ls 9 lw 2  t "70{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}90{/Times \260}"

set label 11 "B"
set yl "Ejecta mass production (kg m^{-2} s^{-1})"
p "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_1.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_2.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_3.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_4.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_5.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_6.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_7.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_8.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_9.txt | sort -k1 -g" u 1:($14*scale):($16*scale) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_1.txt | sort -k1 -g" u 1:($12*scale) w l ls 1 lw 2  t "-90{\260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-70{/E \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_2.txt | sort -k1 -g" u 1:($12*scale) w l ls 2 lw 2  t "-70{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-50{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_3.txt | sort -k1 -g" u 1:($12*scale) w l ls 3 lw 2  t "-50{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-30{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_4.txt | sort -k1 -g" u 1:($12*scale) w l ls 4 lw 2  t "-30{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-10{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_5.txt | sort -k1 -g" u 1:($12*scale) w l ls 5 lw 2  t "-10{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}10{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_6.txt | sort -k1 -g" u 1:($12*scale) w l ls 6 lw 2  t "10{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}30{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_7.txt | sort -k1 -g" u 1:($12*scale) w l ls 7 lw 2  t "30{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}50{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_8.txt | sort -k1 -g" u 1:($12*scale) w l ls 8 lw 2  t "50{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}70{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_9.txt | sort -k1 -g" u 1:($12*scale) w l ls 9 lw 2  t "70{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}90{/Times \260}"

set label 11 "C"
set lmargin 13.5
set yl offset -0.4
set format y "%.1f {\327} 10^6"
set yl "E-folding timescale (yr)"
p "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_1.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_2.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_3.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_4.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_5.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_6.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_7.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_8.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_9.txt | sort -k1 -g" u 1:(1.0/($19*scale*CRAT_DIA*86400.0*365.25*1e6)):(1.0/($21*scale*CRAT_DIA*86400.0*365.25*1e6)) w filledcurves fs transparent solid 0.1 fc rgb "black" notitle,\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_1.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 1 lw 2  t "-90{\260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-70{/E \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_2.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 2 lw 2  t "-70{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-50{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_3.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 3 lw 2  t "-50{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-30{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_4.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 4 lw 2  t "-30{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}-10{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_5.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 5 lw 2  t "-10{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}10{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_6.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 6 lw 2  t "10{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}30{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_7.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 7 lw 2  t "30{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}50{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_8.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 8 lw 2  t "50{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}70{/Times \260}",\
  "<awk 'BEGIN {print \"#\"} {print $0}' ./TAA_169 > tmp; paste tmp Area_Results_9.txt | sort -k1 -g" u 1:(1.0/($17*scale*CRAT_DIA*86400.0*365.25*1e6)) w l ls 9 lw 2  t "70{/Times \260}{/Symbol=32 \243}{/Symbol=32 b}_{}{/Symbol=32 \243}90{/Times \260}"
