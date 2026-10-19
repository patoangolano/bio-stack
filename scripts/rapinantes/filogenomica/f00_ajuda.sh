#!/usr/bin/env bash
datasets download gene gene-id --help 2>&1 | sed -n 1,80p
echo ===========
datasets summary gene taxon --help 2>&1 | sed -n 1,60p
