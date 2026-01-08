from setuptools import setup, Extension
import os

# Define the extension module
mphhc_module = Extension(
    'mphhc',
    sources=['python/mphhc_module.cpp', 'src/core.cpp'],
    include_dirs=['include'],
    language='c++',
    extra_compile_args=['-std=c++17'],
)

setup(
    name='mphhc',
    version='0.1.0',
    description='Python interface for mphhc',
    ext_modules=[mphhc_module],
)
