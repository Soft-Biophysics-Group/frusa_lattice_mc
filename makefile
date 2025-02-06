build:
	mkdir -p build
	cmake . -DMODEL_TYPE="lattice_particles" -B build
	cmake --build build

python/.venv:
	(\
	cd python;\
	python -m venv ".venv";\
	source .venv/bin/activate;\
	python -m pip install -r requirements.txt;\
	python -m pip install -e .;\
	)

all: build python/.venv

.PHONY clean:
	rm -rf build
	rm -rf ./python/.venv/
