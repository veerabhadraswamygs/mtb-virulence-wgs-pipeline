set -u
SP=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad
P=$SP/p2; U="/mnt/user-data/uploads/Downloads/limitation fix files"; S=$P/src; B=$P/bin
cd $S && tar xjf "$U/New folder/got during crosschecks/htslib-1.24.tar.bz2" && echo unpack_ok
cd $S/htslib-1.24 && ./configure --prefix=$P --disable-libcurl --disable-s3 --disable-gcs >/dev/null 2>&1 && make -j2 >/dev/null 2>&1 && make install >/dev/null 2>&1 && echo htslib_ok || echo htslib_FAIL
cd $S/samtools-1.24 && ./configure --prefix=$P --with-htslib=$S/htslib-1.24 --without-curses >/dev/null 2>&1 && make -j2 >/dev/null 2>&1 && make install >/dev/null 2>&1 && echo samtools_ok || echo samtools_FAIL
cd $S/bcftools-1.24 && ./configure --prefix=$P --with-htslib=$S/htslib-1.24 >/dev/null 2>&1 && make -j2 >/dev/null 2>&1 && make install >/dev/null 2>&1 && echo bcftools_ok || echo bcftools_FAIL
echo BUILD2_DONE
