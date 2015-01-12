aa=$(ps axu | grep -E "simulator|main.py" | grep -v "grep" | awk '{printf "%s,", $2}')
echo $aa
top -p ${aa:0:-1}
