PYTHON_VERSION = 3.8
VENV_NAME = venv

TESTSCRIPT=create_tests.py

all: install

install:
	python3 -m venv $(VENV_NAME)
	. $(VENV_NAME)/bin/activate && pip install --upgrade pip
	. $(VENV_NAME)/bin/activate && pip install pandas numpy scipy matplotlib
	. $(VENV_NAME)/bin/activate && 	python3 $(TESTSCRIPT)

clean:
	rm -rf $(VENV_NAME) test*
