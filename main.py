# https://www.kaggle.com/wolfram77/puzzlef-pagerank-ordered-openmp-barrier-vs-barrierfree
import os
from IPython.display import FileLink
src="pagerank-ordered-openmp-barrier-vs-barrierfree"
inp="/kaggle/input/graphs"
out="{}.txt".format(src)
!printf "" > "$out"
display(FileLink(out))
!echo ""

# Download program
!rm -rf $src
!git clone https://github.com/puzzlef/$src

# Run
!g++ -std=c++17 -O3 -fopenmp $src/main.cxx
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/web-Stanford.mtx      2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/web-BerkStan.mtx      2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/web-Google.mtx        2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/web-NotreDame.mtx     2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/soc-Slashdot0811.mtx  2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/soc-Slashdot0902.mtx  2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/soc-Epinions1.mtx     2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/coAuthorsDBLP.mtx     2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/coAuthorsCiteseer.mtx 2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/soc-LiveJournal1.mtx  2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/coPapersCiteseer.mtx  2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/coPapersDBLP.mtx      2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/indochina-2004.mtx    2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/italy_osm.mtx         2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/great-britain_osm.mtx 2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/germany_osm.mtx       2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out $inp/asia_osm.mtx          2>&1 | tee -a "$out"
