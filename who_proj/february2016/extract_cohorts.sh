rm cohort.out
touch cohort.out
for f in auto_output/junk.err*
do
  seed=`grep SCENARIO $f | cut -f2 -d' '` 
  awk -v seed=$seed '{if ($1=="COHORT") {$1=seed; print $0}}' FS="," OFS="," $f >> cohort.out
  echo $seed
done
