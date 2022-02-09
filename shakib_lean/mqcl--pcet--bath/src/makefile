# This is an commentary line in a makefile
# Start of the makefile
pcet.x: mdnonadiab.o adiab_nve.o adiab_nveI.o force.o gaussn.o quantum.o runnvt_NHC.o runnve.o functions.o  
	ifort -o pcet.x mdnonadiab.o adiab_nve.o adiab_nveI.o force.o gaussn.o quantum.o runnvt_NHC.o runnve.o functions.o -mkl=sequential
	mv pcet.x ../
mdnonadiab.o: mdnonadiab.f90
	ifort -c mdnonadiab.f90 -mkl=sequential
adiab_nve.o: adiab_nve.f90
	ifort -c adiab_nve.f90 -mkl=sequential
adiab_nveI.o: adiab_nveI.f90
	ifort -c adiab_nveI.f90 -mkl=sequential
force.o: force.f90
	ifort -c force.f90 -mkl=sequential
gaussn.o: gaussn.f90
	ifort -c gaussn.f90 -mkl=sequential
quantum.o: quantum.f90
	ifort -c quantum.f90 -mkl=sequential
runnvt_NHC.o: runnvt_NHC.f90
	ifort -c runnvt_NHC.f90 -mkl=sequential
#runnvt.o: runnvt.f90
#	ifort -c runnvt.f90 -mkl=sequential
runnve.o: runnve.f90
	ifort -c runnve.f90 -mkl=sequential
functions.o: functions.f90
	ifort -c functions.f90 -mkl=sequential
clean:
	rm mdnonadiab.o adiab_nve.o adiab_nveI.o force.o gaussn.o quantum.o runnvt_NHC.o runnve.o functions.o 
# End of the makefile

