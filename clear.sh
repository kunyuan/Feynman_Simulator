if [ -n "$1" ]
then
  if [ $1 = "-d" ] || [ $1 = "--diag" ]; then
      rm -rf ./diagram/*
  elif [ $1 = "-a" ] || [ $1 = "--all" ]; then
      rm *.pkl
      rm *.txt
  fi
fi
rm *_statis.pkl
rm *_para.txt
rm *.log
rm *.gv
rm -rf infile
rm -rf outfile
rm *.jpg
rm -rf ./diagram/
mkdir ./diagram
