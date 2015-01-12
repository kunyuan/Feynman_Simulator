if [ -n "$1" ]
then
  if [ $1 = "-d" ] || [ $1 = "d" ]; then
      rm -rf ./diagram/*
  fi
else
  rm *_statis.pkl
  rm *_para.txt
  rm *.log
  rm *.gv
  rm -rf infile
  rm -rf outfile
  rm *.jpg
  rm -rf ./diagram/
  mkdir ./diagram
fi
