python -m cProfile -o hotspot.pstats main.py -p 0
../tool/profile/gprof2dot.py -f pstats hotspot.pstats|dot -Tsvg -o hotspot.svg
