#!/bin/bash
cd ..
rm -rf PeakSegDP-release
cp -r PeakSegDP PeakSegDP-release
#grep -v Remotes PeakSegDP/DESCRIPTION > PeakSegDP-release/DESCRIPTION
PKG_TGZ=$(R CMD build PeakSegDP-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
