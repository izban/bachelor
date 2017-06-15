f(x) = a*x + b
fit f(x) 'AverageBinaryLCAF' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f', a, b)
set xlabel "n" 0.0,0.0
set ylabel "lcaf" 0.0,0.0
plot "AverageBinaryLCAF" w l, f(x) t title_f(a,b) with lines
pause mouse
