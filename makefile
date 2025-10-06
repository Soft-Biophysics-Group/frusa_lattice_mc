all: build python/.venv compile_command.json data/ 3dFigures/

build: src/* include/*
	mkdir -p build
	cmake . -DMODEL_TYPE="lattice_particles" -B build
	cmake --build build
	cp build/app/frusa_mc .

compile_command.json: build
	cmake . -DMODEL_TYPE="lattice_particles" -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -B build

.venv: python/requirements.txt
	(\
	python -m venv ".venv";\
	source .venv/bin/activate;\
	python -m pip install -r python/requirements.txt;\
	python -m pip install -e python;\
	)

data/:
	mkdir -p data

3dFigures/:
	mkdir -p 3dFigures

.PHONY clean:
	rm -rf build
	rm -rf ./.venv/
