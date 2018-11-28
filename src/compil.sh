nagfor -c Types_mod.f90
nagfor -c RHS_mod.f90
nagfor -c CFL_mod.f90
nagfor -c IO_mod.f90 -I/usr/local/netcdf-4.6.1/include
nagfor -c Solver_mod.f90
nagfor -c -I/usr/local/netcdf-4.6.1/include fd1d_heat_explicit.f90
nagfor -L/usr/local/netcdf-4.6.1/lib -lnetcdff -lnetcdf fd1d_heat_explicit.o Types_mod.o RHS_mod.o CFL_mod.o IO_mod.o Solver_mod.o -o fd1d_heat_explicit.exe

./fd1d_heat_explicit.exe


diff h_test01.txt h_test01.txt_bak
