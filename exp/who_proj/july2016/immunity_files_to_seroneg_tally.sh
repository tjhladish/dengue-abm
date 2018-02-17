# directories/file locations need to be adjusted; this is just copied from history

#grep end ./auto_output/rep*.err | cut -d' ' -f3,12 | sort -k2 -n | uniq > foi_lookup_120
#for filename in /ufrc/longini/tjhladish/imm_who-baseline-seroprev-july2016/immunity.*$1; do echo $filename; grep -H ' 0 0 0 0 $' $filename | tr [.:] ' ' | cut -f2,4 -d' ' | sort -k2 -n | uniq -c >> /ufrc/longini/tjhladish/seroneg_tally-july2016_120-$1; done
awk '{print $2, $3, $1}' /ufrc/longini/tjhladish/seroneg_tally-july2016_120 > /ufrc/longini/tjhladish/seroneg_tmp-july2016_120
join <(sort -k 1b,1 /ufrc/longini/tjhladish/seroneg_tmp-july2016_120) <(sort -k 1b,1 foi_lookup_120) > seroneg_w_foi_120
