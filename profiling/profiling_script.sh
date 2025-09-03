#!/bin/bash

echo " "
echo "Profiling $1 workflow."
echo " "

# Extract path from configuration file
result_path=$(grep 'result_path:' ../config/config.yml | sed 's/.*: //; s/"//g')
test_path=$(grep 'test_path:' ../config/config.yml  | sed 's/.*: //; s/"//g')
psrecord_path=$(grep 'psrecord_path:' ../config/config.yml  | sed 's/.*: //; s/"//g')

echo "Profiling results saved in: $result_path/$1/profiling/"
echo "Test found in: $test_path$1/"
echo "psrecord imported from $psrecord_path"
echo " "

# export path to psrecord
export PATH=$PATH:$psrecord_path

# check if results path exists
cd $result_path
if [ -d "$1" ]; then
  echo "$1 directory already exists."
else
  mkdir $1
fi

# run psrecord profiling
echo "Starting the profiling session..."
cd $test_path
# python $1/$1.py
psrecord "python $1/$1.py" --log $result_path$1/profiling/activity.txt --include-io --interval 0.1 --include-children 

echo "Formatting outputs..."
# output results in required format
cd profiling/
python ./monitor_activity.py $1 

echo "Done."
exit 0
