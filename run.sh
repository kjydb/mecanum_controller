#!/usr/bin/env bash

set -e
SCRIPT_DIR=$(dirname ${BASH_SOURCE})

DEBUG=false
for arg in "$@"; do
    if [ "$arg" == "--debug" ]; then
        DEBUG=true
        break
    fi
done

if [ "$DEBUG" == "true" ]; then
    gdb $SCRIPT_DIR/build/"$1"
else
    $SCRIPT_DIR/build/"$1" "$@"
fi
