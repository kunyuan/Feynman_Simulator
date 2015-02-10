#clang -bundle -undefined dynamic_lookup -L/usr/local/lib -L/usr/local/opt/sqlite/lib -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk build/temp.macosx-10.10-x86_64-2.7/lu_fast.o -o /Users/kchen/project/Feynman_Simulator/dyson/solver/lu_fast.so -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
python setup.py build_ext --inplace

