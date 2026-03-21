set datafile separator ","

# -------- grafica solucion --------
set terminal png size 1200,800
set output "solucion.png"

set title "Solucion numerica vs solucion exacta"
set xlabel "x"
set ylabel "y"
set grid

plot "datos.csv" using 1:2 with lines lw 2 title "Numerica", \
     "datos.csv" using 1:3 with lines lw 2 title "Exacta"

# -------- grafica error --------
set output "error.png"

set title "Error "
set xlabel "x"
set ylabel "Error"
set grid

plot "datos.csv" using 1:4 with lines lw 2 title "Error"

set output
