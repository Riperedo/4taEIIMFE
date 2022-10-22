set xrange [0 to 30]
set yrange [0 to 2.5]

set xlabel "k"
set ylabel "F(k,{/Symbol t})"

plot 1 lt "--" lc rgb "#fde725" t "", "sdk.dat" w l lc rgb "#fde725" t "{/Symbol t} = 0", "ISF_tau_0p001.dat" u 1:2 w l lt "--" lc rgb "#5ec962" t "", "" u 1:3 w l lc rgb "#5ec962" t "{/Symbol t} = 10^{-3}", "ISF_tau_0p01.dat" u 1:2 w l lc rgb "#21918c" lt "--" t "", "" u 1:3 w l lc rgb "#21918c" t "{/Symbol t} = 10^{-2}", "ISF_tau_0p1.dat" u 1:2 w l lt "--" lc rgb "#3b528b" t "", "" u 1:3 w l lc rgb "#3b528b" t "{/Symbol t} = 10^{-1}", "ISF_tau_1p0.dat" u 1:2 w l lt "--" lc rgb "#440154" t "", "" u 1:3 w l lc rgb "#440154" t "{/Symbol t} = 10^{0}", "ISF_tau_10p0.dat" u 1:2 w l lt "--" lc rgb "#000000" t "", "" u 1:3 w l lc rgb "#000000" t "{/Symbol t} = 10^{1}"
