#!/usr/bin/env bash
# Build libinchi from source for platforms with no committed binary in INCHI/.
#
# Every platform except linux aarch64 ships a prebuilt binary, so this is a no-op
# there. Runs inside the cibuildwheel manylinux container via CIBW_BEFORE_ALL_LINUX.
set -euo pipefail

arch="$(uname -m)"
if [ "$arch" != "aarch64" ]; then
    echo "libinchi: $arch uses the committed binary, nothing to build"
    exit 0
fi

target="INCHI/libinchi_aarch64.so"
if [ -f "$target" ]; then
    echo "libinchi: $target already present"
    exit 0
fi

# pinned to the release upstream chython builds against
version="v1.07.5"

command -v cmake >/dev/null 2>&1 || pip install cmake
command -v git >/dev/null 2>&1 || { yum install -y git || apt-get update && apt-get install -y git; }

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

git clone --depth 1 --branch "$version" https://github.com/IUPAC-InChI/InChI "$tmp/InChI"

src="$tmp/InChI/INCHI-1-SRC/INCHI_API/libinchi/src"
[ -d "$src" ] || { echo "libinchi: source tree missing at $src" >&2; exit 1; }

cmake -S "$src" -B "$tmp/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$tmp/build" --config Release --target libinchi -j "$(nproc)"

built="$(find "$tmp/build" -name 'libinchi*.so*' -type f | head -1)"
[ -n "$built" ] || { echo "libinchi: build produced no shared library" >&2; exit 1; }

mkdir -p INCHI
cp "$built" "$target"
echo "libinchi: built $target ($(uname -m))"
