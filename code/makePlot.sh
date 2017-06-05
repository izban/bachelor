f(x) = a*x + b
fit f(x) 'avlcas.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f', a, b)
plot "avlcas.dat" u 1:2 w l, f(x) t title_f(a,b) with points
pause mouse
