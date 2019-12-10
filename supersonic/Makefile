# This file manages the compilation of the Euler code. 
# All lines start with "#" are comment lines not considered 
# by the operation system

OUT = output.o output_hg.o output_mat.o

OBS = common_block.o euler.o apply_bconds.o  set_others.o check_conv.o set_timestep.o check_grid.o smooth.o crude_guess.o flow_guess.o read_data.o  sum_fluxes.o  generate_grid.o  set_fluxes.o new_guess.o $(OUT) 

# Choice of comipler
FC = gfortran

# Option of Compilation

# debug option
FFLAGS = -g

# select the optimising option by delete # in the following line
#FFLAGS = -O2

# comiple the executable and name it "Euler" using "-o"
euler : $(OBS)
	$(FC)  $(FFLAGS) -o Euler  $(OBS)  

# remove all .o files and Euler for a new start
clean:
	rm -f *.o  Euler

# compile .f90 to .o
%.o:%.f90
	$(FC) $(FFLAGS) -c $^ -o $@
