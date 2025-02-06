build:
	mkdir -p build
	cmake . -DMODEL_TYPE="lattice_particles" -B build
	cmake --build build

python/.venv:
	(\
	python -m venv "python/.venv";\
	source python/.venv/bin/activate;\
	python -m pip install requirements.txt\
	python -m pip install -e python\
	)

all: build python/.venv

.PHONY clean:
	rm -rf build
	rm -rf ./python/.venv/
