#!/usr/bin/env bash
# Regenerate the API documentation (HTML, with search) into docs/.
# Run from the repository root:  ./docs/build.sh

pdoc -o docs \
    coltess \
    coltess.core \
    coltess.catalog \
    coltess.download \
    coltess.photometry \
    coltess.analysis \
    coltess.parallel \
    coltess.utils
