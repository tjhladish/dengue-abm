attach '/home/tjhladish/work/dengue/exp/tirs_study/transfer/tirs_trial-trans_loc_sensitivity-v3.0b.sqlite' as toMerge;
BEGIN; 
insert into job select * from toMerge.job; 
insert into par select * from toMerge.par; 
insert into met select * from toMerge.met; 
COMMIT; 
detach toMerge;

attach '/home/tjhladish/work/dengue/exp/tirs_study/transfer/tirs_trial-trans_loc_sensitivity-v3.0c.sqlite' as toMerge;
BEGIN; 
insert into job select * from toMerge.job; 
insert into par select * from toMerge.par; 
insert into met select * from toMerge.met; 
COMMIT; 
detach toMerge;
