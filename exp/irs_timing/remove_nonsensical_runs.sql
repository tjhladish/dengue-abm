delete from par where vector_control=1 and campaign_duration=2 and timing != 0;
delete from par where vector_control=0 and (campaign_duration!=0 or timing != 0 or vc_coverage != 0.25);
delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);
vacuum;
