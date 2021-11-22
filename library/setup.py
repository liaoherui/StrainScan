from setuptools import setup, Extension

seqpy_module = Extension('seqpy', sources = ["seqpy.c"])

setup(
    name = "seqpy",
    description = "seqpy",
    ext_modules = [seqpy_module])
