rm cohort.out
touch cohort.out
for f in auto_output/reg*.err
#while read f;
do
  #seed=`grep SCENARIO $f | cut -f2 -d' '` 
  #echo $f, $seed
  #awk -v seed=$seed '{if ($1=="COHORT") {$1=seed; print $0}}' FS="," OFS="," $f
  grep COHORT $f | cut -f2-7 -d',' >> cohort.out
  #awk '{if ($1=="COHORT") {$1=""; print $0}}' FS="," OFS="," $f >> cohort.out
  #echo $seed
done #<alt_hits
