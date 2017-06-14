f(x) = a*x*x*log(x)
g(x) = t*x*x*x
fit f(x) 'lcasFast.dat' u 1:2 via a,b,c
fit g(x) 'lcasAttabi.dat' u 1:2 via t
title_f(a,b,c) = sprintf('f(x) = %fx^{1.86}', a, b, c)
plot "lcasFast.dat" u 1:2 w l, f(x) t title_f(a,b,c), "lcasAttabi.dat" u 1:2 w l, g(x)
pause mouse
