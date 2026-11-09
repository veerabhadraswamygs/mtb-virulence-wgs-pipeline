#!/bin/bash
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
cd $P/src
# CD-HIT
( cd cdhit-4.8.1 && make -j2 openmp=no >$P/logs/build_cdhit.log 2>&1 && for b in cd-hit cd-hit-est cd-hit-2d; do cp $b $P/bin/; done && echo CDHIT_OK >>$P/logs/build_cdhit.log )
# MAFFT
( cd mafft-main/core && sed -i "s|^PREFIX = .*|PREFIX = $P|" Makefile && make -j2 >$P/logs/build_mafft.log 2>&1 && make install >>$P/logs/build_mafft.log 2>&1 && echo MAFFT_OK >>$P/logs/build_mafft.log )
# HyPhy
( cd hyphy-2.5.101 && mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$P -DNOAVX=ON .. >$P/logs/build_hyphy.log 2>&1 && make -j2 HYPHYMP >>$P/logs/build_hyphy.log 2>&1 && make install >>$P/logs/build_hyphy.log 2>&1 && echo HYPHY_OK >>$P/logs/build_hyphy.log )
