--delete from par where vac = 0 and (vac_rate != 0);
--delete from job where serial not in (select serial from par);
--delete from met where serial not in (select serial from par);
delete from par where wide_vc_campaign = 0 and vc_coverage != 1;
delete from par where wide_vc_campaign = 1 and (vc_coverage != 0.5 and vc_coverage != 0.7 and vc_coverage != 0.9)
