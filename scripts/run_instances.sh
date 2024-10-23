#!/bin/bash

OUTPUT_DIR=$HOME/Q-HOPE/logs9_zipped/
mkdir -p $OUTPUT_DIR

#FILES=$HOME/Q-HOPE/instances_auth/nodes_*
FILES=$HOME/Q-HOPE/instances/turkey2/nodes_*
for f in $FILES
do
  echo "Processing $f file..."
  f1=${f/nodes/edges}

  OUTPUT_FILE=$OUTPUT_DIR/$(basename $f)-system-warmstart.log
  if [ -f "$OUTPUT_FILE" ]; then
    echo "$OUTPUT_FILE does exist."
  else
    echo "Running for $OUTPUT_FILE"
    timeout 200000 ./leblanc_solver $f $f1 FR > $OUTPUT_FILE
    #timeout 17000 ./leblanc_solver $f $f1 FR > $OUTPUT_FILE
  fi
done
