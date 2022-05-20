delete from par where gen_vc = 0 and vc_coverage != 1.0;
delete from par where gen_vc = 1 and not (vc_coverage = 0.5 or vc_coverage = 0.7 or vc_coverage = 0.9);
delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);
update par set seed = (select posterior + 1 from job where par.serial = job.serial);
vacuum;
analyze;
