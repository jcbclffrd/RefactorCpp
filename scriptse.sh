latex plot.tex
dvips plot.dvi
#gv plot.ps
ps2pdf14 plot.ps oData/plot.pdf
rm plot.aux plot.dvi plot.log plot.ps
