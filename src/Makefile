FFLAGS = 
# add any other compilation flag
LDFLAGS = -L. -lfmw # add any other link flag

fd1d_heat_explicit.exe: fd1d_heat_explicit.o
	nagfor $^ $(LDFLAGS) -o $@

fd1d_heat_explicit.o: fd1d_heat_explicit.f90
	nagfor $(FFLAGS) $(LDFLAGS)  -c $<

.PHONY: clean
clean:
	rm -rf *.o *.mod main.exe
