PYTHON_VERSION = 3.8
VENV_NAME = venv

all: install

install:
	python3 -m venv $(VENV_NAME)
	. $(VENV_NAME)/bin/activate && pip install --upgrade pip
	. $(VENV_NAME)/bin/activate && pip install numpy scipy

clean:
	rm -rf $(VENV_NAME)
