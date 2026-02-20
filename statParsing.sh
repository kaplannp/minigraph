FILE="gwf-ed-stats.log"
cat $FILE | grep n: | grep -oP "\d*$" > ns.txt
cat $FILE | grep -P "SIZE A" | grep -oP "\d*$" > as.txt
cat $FILE | grep -P "INTV_BYTES" | grep -oP "\d*$" > intv.txt
cat $FILE | grep -P "SCORE" | grep -oP "\d*$" > scores.txt
cat $FILE | grep -P "N_INACTIVE" | grep -oP "\d*$" > nInactive.txt
cat $FILE | grep -P "HA_ENTRIES" | grep -oP "\d*$" > ha_sizes.txt
cat $FILE | grep -P "HA_NODES" | grep -oP "\d*$" > ha_nodes.txt
cat $FILE | grep -P "BYTES_FULL" | grep -oP "\d*$" > graphSizeFull.txt
cat $FILE | grep -P "BYTES_ESSENTIAL" | grep -oP "\d*$" > graphSizeEssential.txt
cat $FILE | grep -P "^BYTES_COMPRESSED" | grep -oP "\d*$" > graphSize2BitSeq.txt
cat $FILE | grep -P "TOTAL_SEGMENTS" | grep -oP "\d*$" > graphNodes.txt
cat $FILE | grep -P "TOTAL_EDGES" | grep -oP "\d*$" > graphEdges.txt
#cat $FILE | grep "LENGTH_BP:" | grep -oP "\d*$" > graphNodeLen.txt

#BFS stuff
cat $FILE | grep -P "BFS_TIME_MS" | grep -oP "\d*$" > bfsTime.txt
cat $FILE | grep -P "BFS_SEGMENTS" | grep -oP "\d*$" > bfsNodes.txt
cat $FILE | grep -P "BFS_EDGES" | grep -oP "\d*$" > bfsEdges.txt
cat $FILE | grep -P "BFS_BYTES_COMPRESSED" | grep -oP "\d*$" > bfsSize2BitSeq.txt


