#!/bin/sh
mkdir -p ../data
curl https://archive.eso.org/ASTROM/TYC-2/data/catalog.dat > ../data/tycho2_catalog.dat
curl https://archive.eso.org/ASTROM/TYC-2/data/suppl_1.dat > ../data/tycho2_suppl_1.dat
