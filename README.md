# Cshanty

Solar sailing, slightly speedier.

Cshanty is a re-implementation of my undergraduate thesis code in C (with a Python wrapper) for greater performance.

# Requirements
* Python 3.10+
* A C compiler installed on your system (Tested with GCC on Linux and MSVC on Windows)

# Installation
1. Clone the repository
2. (Optional) Create a virtual environment
3. Run `pip install -e .` from the root directory to automatically install the package and its dependencies (compiles C code in the process)
4. Run the script `cshanty/python_demo.py` to run a test propagation