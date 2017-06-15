f(x) = a*(x**1.86)
fit f(x) '3SumUnaryTest' u 1:2 via a
title_f(a) = sprintf('f(x) = %ex^{1.86}', a)
set xlabel "n"
set ylabel "time, seconds"
plot "3SumUnaryTest" w l, f(x) t title_f(a)
pause mouse
