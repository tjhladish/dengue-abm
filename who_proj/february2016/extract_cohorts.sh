rm cohort.out
touch cohort.out
for f in auto_output/newjunk.err*
#while read f;
do
  seed=`grep SCENARIO auto_output/$f | cut -f2 -d' '` 
  awk -v seed=$seed '{if ($1=="COHORT") {$1=seed; print $0}}' FS="," OFS="," auto_output/$f >> cohort.out
  echo $seed
done #<alt_hits
