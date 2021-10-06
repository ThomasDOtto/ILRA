#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

gawk -f $SCRIPT_DIR/VSlistTo1HitPerLine.awk "$@"

