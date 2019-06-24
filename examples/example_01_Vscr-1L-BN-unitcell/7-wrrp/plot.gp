set term png enhanced size 1024,768 font "Helvetica,30"
set output "vscr.plavg.png"

set title "Planar average along z of Vscr(r,r) in 1L BN"
set xlabel "z [au]"
set ylabel "Vscr [eV]"

set arrow from 6.141609463,graph 0 to 6.141609463,graph 1 nohead
set grid
plot "wrrp.real.out" u (($1-1)*0.708647245725000):($2*13.6) w lp 
