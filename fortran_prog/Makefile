FC = gfortran

FFLAGS = -O2
DEBUGFLAGS = -g -fcheck=all -Wall -fbacktrace

SRC = spline_cubic.F90 DUO_inst.F90 main_intsDUO.F90
OBJ = $(SRC:.F90=.o)

TARGET = test

# --------------------------------------------------
# NORMAL
# --------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $(TARGET)

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

# --------------------------------------------------
# DEBUG
# --------------------------------------------------
debug: FFLAGS = $(DEBUGFLAGS)
debug: clean $(TARGET)

# --------------------------------------------------
# CLEAN
# --------------------------------------------------
clean:
	rm -f $(OBJ) *.mod $(TARGET)