#Template of makefile with multiple files
# Vincent Ouazan-Reboul, 20/11/2023
CC = g++
CFLAGS = -pedantic-errors -Weffc++ -Wextra -Wconversion -Wsign-conversion -Werror -std=c++2a
# Target file
BUILD_DIR = ./build
TRGT = lattice_mc
SRC_DIRS = ./app ./src/utility ./src/engine ./src/models ./src/models/particles
# SRCS := $(shell find $(SRC_DIRS) -name '*.cc' -or -name '*.c' -or -name '*.s')

SRCS = main.cc lattice_particles_state.cc lattice_particles_parameters.cc lattice_particles_geometry.cc lattice_particles_hexagonal.cc lattice_particles_interactions.cc lattice_particles_update.cc ../../utility/vector_utils.cc
#OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
OBJS = $(patsubst %.cpp,$(BUILD_DIR)/%.o, $(SRCS))
INCL_PARENT = ./include
INCL = -I$(INCL_PARENT)/models
INCL += -I$(INCL_PARENT)/models/particles
INCL += -I$(INCL_PARENT)/thirdparty
INCL += -I$(INCL_PARENT)/engine
INCL += -I$(INCL_PARENT)/utility

# $(shell echo $(INCL))

# Ensure the build directory exists
$(shell mkdir -p $(BUILD_DIR))

$(OBJS): $(SRCS) $(INCL_PARENT)
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(TRGT): $(OBJS)
	$(CC) $(CFLAGS) $(INCL) -ggdb -o $(TRGT) $(OBJS)
	# Comment out if you don't want the compiled file here

.PHONY: rollout
rollout: $(OBJS)
	$(CC) $(CFLAGS) $(INCL) -o $(TRGT) $(OBJS)
	# Comment out if you don't want the compiled file here
	cp $(TRGT)


.PHONY: clean
clean:
	rm -rf build
