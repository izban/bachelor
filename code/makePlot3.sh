f(x) = a*x*x + b*x + c
fit f(x) '3SumRandom.dat' u 1:2 via a, b, c
title_f(a,b,c) = sprintf('f(x) = %fx^2 + %fx + %f', a, b, c)
plot "3SumRandom.dat" u 1:2 w l, f(x) t title_f(a,b,c)
pause mouse
