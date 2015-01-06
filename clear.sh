if [ -n "$1" ]
then
  if [ $1 = "-d" ] || [ $1 = "d" ]; then
      rm -rf ./data/diagram/*
  fi
else
  rm -rf ./data/*
  mkdir ./data/diagram
fi
