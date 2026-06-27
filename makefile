all: build .venv compile_command.json data/ 3dFigures/

build: src/* include/*
	mkdir -p build
	cmake . -DMODEL_TYPE="lattice_particles" -B build
	cmake --build build
	cp build/app/frusa_mc .

compile_command.json: build
	cmake . -DMODEL_TYPE="lattice_particles" -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -B build

# Bootstrap uv via pip if absent, then let uv own ./.venv at 3.11.
.venv: pyproject.toml uv.lock
	@command -v uv >/dev/null 2>&1 || python3 -m pip install --user uv
	uv sync

data/:
	mkdir -p data

3dFigures/:
	mkdir -p 3dFigures

.PHONY clean:
	rm -rf build
	rm -rf ./.venv/
