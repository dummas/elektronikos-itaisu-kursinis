ANI = 3*10**20;
ANZ0 = 2.5*10**(19)
TNZ = 1000
DNZ = 1.9*10**(-13)
WZ = 2.2
TDZ = 24
TNR = 1000
DNR = 1.4*10**(-13)
WR = 3.7
TDR = 77
TZ = 942
TR = 1056
ANE0 = 1.4*10**(21)
TNE = 1000
DNE = 3.8*10**(-13)
WE = 2.9
TDE = 42
TE = 318


QE = 0.16E-18
AK = 0.138E-22
S = QE/AK
TDZ = TDZ*60
TDR = TDR*60
TDE = TDE*60
DZ = DNR/exp(S*WZ*(1/TZ-1/TNZ))
DR = DNR/exp(S*WZ*(1/TR-1/TNR))
DE = DNR/exp(S*WZ*(1/TE-1/TNE))

Q = 1.13*ANZ0*sqrt(DZ*TDZ)

ANB(x) = Q/sqrt(3.14*(DR*TDR+DE*TDE)) * exp(-(x/(2*sqrt(DR*TDR+DE*TDE)))**2)
ANE(x) = ANE0*erf(x/(2*sqrt(DE*TDE)))
AN(x) = ANI-ANB(x)+ANE(x)
set xr [0.0:10]
set term table
set output "dvitakte_difuzija.plt"

plot(log10(abs(AN(x))))
