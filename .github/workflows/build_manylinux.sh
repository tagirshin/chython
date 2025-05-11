#!/bin/bash
set -e -x

PY_VERSION_SHORT=$1 # e.g., "3.10"
PROJECT_DIR=$2    # e.g., /workdir

# Map to full path in manylinux image, e.g., /opt/python/cp310-cp310/bin/python
PY_VERSION_NODOTS=$(echo "${PY_VERSION_SHORT}" | tr -d '.')
PYTHON_EXE="/opt/python/cp${PY_VERSION_NODOTS}-cp${PY_VERSION_NODOTS}/bin/python"
PIP_EXE="/opt/python/cp${PY_VERSION_NODOTS}-cp${PY_VERSION_NODOTS}/bin/pip"

echo "Attempting to use Python at ${PYTHON_EXE} for version ${PY_VERSION_SHORT}"
echo "Available Pythons in /opt/python/:"
ls -l /opt/python/

# Check if Python exists
if [ ! -f "$PYTHON_EXE" ]; then
    echo "Python executable $PYTHON_EXE not found in Docker image."
    exit 1
fi

# Install poetry using this specific pip
echo "Installing Poetry using ${PIP_EXE}..."
"$PIP_EXE" install poetry

# Navigate to project
cd "$PROJECT_DIR"

# Add the specific python's bin to PATH so poetry uses it
export PATH="/opt/python/cp${PY_VERSION_NODOTS}-cp${PY_VERSION_NODOTS}/bin:$PATH"
echo "Updated PATH: $PATH"
echo "Using python: $(which python)"
echo "Python version: $(python --version)"
echo "Using poetry: $(which poetry)"
echo "Poetry version: $(poetry --version)"

# Build wheel
echo "Building wheel with Poetry..."
poetry build -f wheel

echo "Build complete. Wheels should be in $PROJECT_DIR/dist"
# If needed, add chown command here to fix permissions of files in dist/
# e.g., chown -R $(stat -c '%u:%g' .) dist/
# For now, assuming default permissions are handled by subsequent GHA steps. 