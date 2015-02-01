if [ -n "$1" ]
then
  if [ $1 = "-d" ] || [ $1 = "--diag" ]; then
      rm -rf ./diagram/*
  elif [ $1 = "-a" ] || [ $1 = "--all" ]; then
      rm *.hkl
      rm *.txt
  fi
fi
rm Message.txt
rm statis_total.hkl
rm *_statis.hkl
rm *_para.txt
rm Coordinates.txt
rm *.log
rm *.gv
rm -rf infile
rm -rf outfile
rm *.jpg
rm -rf ./diagram/
mkdir ./diagram
