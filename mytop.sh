aa=$(ps axu | grep -E "run_loop|gamma3" | grep -v "grep" | awk '{printf "%s,", $2}')
echo $aa
top -p ${aa:0:-1}
