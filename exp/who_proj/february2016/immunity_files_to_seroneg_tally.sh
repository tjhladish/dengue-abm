# directories/file locations need to be adjusted; this is just copied from history

grep end auto_output/burn80-sero.err-* | cut -d' ' -f3,12 | sort -k2 -n | uniq > foi_lookup
for filename in ./immunity.*; do echo $filename; grep -H ' 0 0 0 0 $' $filename | tr [.:] ' ' | cut -f3,5 -d' ' | sort -k2 -n | uniq -c >> ../seroneg_tally; done
awk '{print $2, $3, $1}' seroneg_tally > seroneg_tmp
join <(sort -k 1b,1 seroneg_tmp) <(sort -k 1b,1 foi_lookup) > seroneg_w_foi
