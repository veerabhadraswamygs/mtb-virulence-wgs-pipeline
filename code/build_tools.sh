set -u
SP=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad
P=$SP/p2; U="/mnt/user-data/uploads/Downloads/limitation fix files"; S=$P/src; B=$P/bin
cd $S
echo "== unpack =="
tar xjf "$U/htslib-1.24.tar.bz2"
tar xjf "$U/samtools-1.24.tar.bz2"
tar xjf "$U/New folder/got during crosschecks/bcftools-1.24.tar.bz2"
tar xzf "$U/New folder/got during crosschecks/bwa-0.7.18.tar.gz"
tar xjf "$U/minimap2-2.31_x64-linux.tar.bz2"
tar xzf "$U/diamond-linux64.tar.gz"
tar xzf "$U/New folder/got during crosschecks/iqtree-2.4.0-Linux-intel.tar.gz"
gunzip -c "$U/freebayes-1.3.10-linux-amd64-static.gz" > $B/freebayes
cp "$U/prodigal.linux" $B/prodigal
chmod +x $B/freebayes $B/prodigal
cp $S/minimap2-2.31_x64-linux/minimap2 $B/ 2>/dev/null
cp $S/diamond $B/ 2>/dev/null
cp $S/iqtree-2.4.0-Linux-intel/bin/iqtree2 $B/ 2>/dev/null
chmod +x $B/* 2>/dev/null
echo "== htslib =="
cd $S/htslib-1.24 && ./configure --prefix=$P --disable-libcurl --disable-s3 --disable-gcs >/dev/null 2>&1 && make -j2 >/dev/null 2>&1 && make install >/dev/null 2>&1 && echo htslib_ok || echo htslib_FAIL
echo "== samtools =="
cd $S/samtools-1.24 && ./configure --prefix=$P --with-htslib=$S/htslib-1.24 --without-curses >/dev/null 2>&1 && make -j2 >/dev/null 2>&1 && make install >/dev/null 2>&1 && echo samtools_ok || echo samtools_FAIL
echo "== bcftools =="
cd $S/bcftools-1.24 && ./configure --prefix=$P --with-htslib=$S/htslib-1.24 >/dev/null 2>&1 && make -j2 >/dev/null 2>&1 && make install >/dev/null 2>&1 && echo bcftools_ok || echo bcftools_FAIL
echo "== bwa =="
cd $S/bwa-0.7.18 && make -j2 >/dev/null 2>&1 && cp bwa $B/ && echo bwa_ok || echo bwa_FAIL
echo "BUILD_DONE"
