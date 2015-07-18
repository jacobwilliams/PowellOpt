#  please change the variables FC and FFLAGS to match your fortran compiler

FC = gfortran
FFLAGS = -g

main: addcon.o adjtol.o conres.o delcon.o eqcons.o fgcalc.o getd.o getfes.o \
      getmin.o initzu.o ktvec.o lsrch.o minflc.o minfun.o newcon.o satact.o \
      sdegen.o sdirn.o stepbd.o zbfgs.o main.o
	$(FC) $(FFLAGS) -o main addcon.o adjtol.o conres.o delcon.o eqcons.o \
          fgcalc.o getd.o getfes.o getmin.o initzu.o ktvec.o lsrch.o minflc.o \
          minfun.o newcon.o satact.o sdegen.o sdirn.o stepbd.o zbfgs.o main.o
test: main
	main > output ; diff output output.txt ; rm -f output
addcon.o: addcon.f
	$(FC) $(FFLAGS) -c addcon.f
adjtol.o: adjtol.f
	$(FC) $(FFLAGS) -c adjtol.f
conres.o: conres.f
	$(FC) $(FFLAGS) -c conres.f
delcon.o: delcon.f
	$(FC) $(FFLAGS) -c delcon.f
eqcons.o: eqcons.f
	$(FC) $(FFLAGS) -c eqcons.f
fgcalc.o: fgcalc.f
	$(FC) $(FFLAGS) -c fgcalc.f
getd.o: getd.f
	$(FC) $(FFLAGS) -c getd.f
getfes.o: getfes.f
	$(FC) $(FFLAGS) -c getfes.f
getmin.o: getmin.f
	$(FC) $(FFLAGS) -c getmin.f
initzu.o: initzu.f
	$(FC) $(FFLAGS) -c initzu.f
ktvec.o: ktvec.f
	$(FC) $(FFLAGS) -c ktvec.f
lsrch.o: lsrch.f
	$(FC) $(FFLAGS) -c lsrch.f
main.o: main.f
	$(FC) $(FFLAGS) -c main.f
minflc.o: minflc.f
	$(FC) $(FFLAGS) -c minflc.f
minfun.o: minfun.f
	$(FC) $(FFLAGS) -c minfun.f
newcon.o: newcon.f
	$(FC) $(FFLAGS) -c newcon.f
satact.o: satact.f
	$(FC) $(FFLAGS) -c satact.f
sdegen.o: sdegen.f
	$(FC) $(FFLAGS) -c sdegen.f
sdirn.o: sdirn.f
	$(FC) $(FFLAGS) -c sdirn.f
stepbd.o: stepbd.f
	$(FC) $(FFLAGS) -c stepbd.f
zbfgs.o: zbfgs.f
	$(FC) $(FFLAGS) -c zbfgs.f
clean:
	rm -f *.o main
