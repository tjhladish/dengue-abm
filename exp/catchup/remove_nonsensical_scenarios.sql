delete from par where (vac_mech=0 and vac_target=2) or (vac_mech=1 and vac_target=9) or (vac_mech=0 and catchup_to=100) or (vac_mech=1 and catchup_to=50);
delete from par where vac=0 and not (vac_mech=0 and catchup=0);
delete from par where vector_control=0 and not vc_coverage=0.75;
delete from job where serial not in (select serial from par);
delete from met where serial not in (select serial from par);
vacuum;
