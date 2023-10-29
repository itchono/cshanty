from pathlib import PurePath
from setuptools import setup

if __name__ == "__main__":
    # Make a relative path for the cffi module
    ext_dir = PurePath("cshanty", "_builder.py:ffi")
    setup(cffi_modules=[str(ext_dir)])
