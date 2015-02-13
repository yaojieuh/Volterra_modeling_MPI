TARGET  = FWPSModeling
CC      = mpicc
CFLAGS  = -Wall -g -openmp
OPT     = -O3
LIBS    = -lm

DEFINES = -DVERB -DTEST_VEL

INCLUDE = -I./include

BIN_PTH = ./bin
SRC_PTH = ./source
SOURCE  = $(SRC_PTH)/Test_VolterraModeling.c $(SRC_PTH)/params.c\
          $(SRC_PTH)/comput_fourierb.c $(SRC_PTH)/comput_volterra.c\
          $(SRC_PTH)/comput_reflection.c $(SRC_PTH)/matrice.c

OBJ     = $(SOURCE:.c=.o) 

## Default rule executed
all: $(TARGET)
	@true

$(TARGET): $(OBJ)
	@echo
	@echo "=> Linking the target $@"
	@$(CC) $(CFLAGS) $(OPT) -o $(BIN_PTH)/$@ $^ $(LIBS)
	@echo '=> Link done'
	@-rm -f $(OBJ)

%.o: %.c
	@mkdir -p $(dir $@)
	@echo
	@echo "=> Compiling $<"
	@$(CC) $(CFLAGS) $(OPT) $(DEFINES) -c $< -o $@ $(INCLUDE)

# the rule to clean binaries
clean:
	@-rm -f bin/$(TARGET) $(OBJ)
	@echo "=> Clean done"
