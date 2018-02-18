grep end auto_output/abc_merida.err-* | cut -d' ' -f4,22,27- | sort -k 1b,1 | uniq > all_seroprevs_set2
echo 'select serial from jobs where posterior > -1 and smcSet=2;' | sqlite3 merida-all_historical_serotypes.sqlite | sort -k 1b,1 > post_serials_set2
join post_serials_set2 all_seroprevs_set2 > tmp2
