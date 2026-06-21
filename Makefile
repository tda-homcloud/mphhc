.PHONY: build build-cpp build-python test test-cpp test-python format format-cpp format-python

build: build-cpp build-python

test: test-cpp test-python

build-cpp:
	cmake --build build

test-cpp: build-cpp
	ctest --output-on-failure --test-dir build

build-python:
	pip install -e .

test-python: build-python
	pytest python-tests

format: format-cpp format-python

CPPFILES := $(shell find src -type f -name '*.cpp') $(shell find tests -type f -name '*.hpp')

format-cpp:
	clang-format -i --style=Google src/*.cpp include/mphhc/*.hpp tests/*.cpp python/*.cpp

format-python:
	black python-tests/*.py benchmark/*.py
