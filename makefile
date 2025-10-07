all: build python/.venv compile_command.json

build: src/* include/*
	mkdir -p build
	cmake . -DMODEL_TYPE="lattice_particles" -B build
	cmake --build build

compile_command.json: build
	cmake . -DMODEL_TYPE="lattice_particles" -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -B build

python/.venv: python/requirements.txt
	(\
	cd python;\
	python -m venv ".venv";\
	source .venv/bin/activate;\
	python -m pip install -r requirements.txt;\
	python -m pip install -e .;\
	)

.PHONY clean:
	rm -rf build
	rm -rf ./python/.venv/
