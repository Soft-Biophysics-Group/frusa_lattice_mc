#Template of makefile with multiple files
# Vincent Ouazan-Reboul, 20/11/2023
CC = g++
CFLAGS = -pedantic-errors -Weffc++ -Wextra -Wconversion -Wsign-conversion -Werror -std=c++2a
# Target file
BUILD_DIR = build
TRGT = lattice_mc
SRCS = main.cc lattice_particles_state.cc lattice_particles_parameters.cc
OBJS = $(patsubst %.cpp,$(BUILD_DIR)/%.o, $(SRCS))
INCL_PARENT = ../../../include
INCL = -I$(INCL_PARENT)/models/lattice_particles
INCL += -I$(INCL_PARENT)/thirdparty
INCL += -I$(INCL_PARENT)/engine
INCL += -I$(INCL_PARENT)/utility

# Ensure the build directory exists
$(shell mkdir -p $(BUILD_DIR))

$(TRGT): $(OBJS)
	$(CC) $(CFLAGS) $(INCL) -ggdb -o $(TRGT) $(OBJS)
	# Comment out if you don't want the compiled file here

.PHONY: rollout
rollout: $(OBJS)
	$(CC) $(CFLAGS) $(INCL) -o $(TRGT) $(OBJS)
	# Comment out if you don't want the compiled file here
	cp $(TRGT)

$(BUILD_DIR)/%.o: %.cpp $(INCL_PARENT)
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

.PHONY: clean
clean:
	rm -rf build
