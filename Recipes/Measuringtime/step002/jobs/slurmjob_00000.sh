#! /bin/bash
#SBATCH --job-name=Measuringtime_002
#SBATCH --output=/rwthfs/rz/cluster/home/bo791269/Software/AntiprotonAnalysis/Recipes/Measuringtime/step002/output/output_00000.txt
#SBATCH --error=/rwthfs/rz/cluster/home/bo791269/Software/AntiprotonAnalysis/Recipes/Measuringtime/step002/output/output_00000.txt
#SBATCH --time=23:00:00
#SBATCH --partition=c18m
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --mem-per-cpu=3882M
#SBATCH --account=jara0052

source /rwthfs/rz/cluster/home/bo791269/Software/AntiprotonAnalysis/Recipes/Measuringtime/step002/sandbox/environment.sh
export PSI_EXPORTS=ACQTDATADIR,ACSOFTDIR,ACSOFTLOOKUPS,ACSOFT_ADDITIONAL_LOGONS,AMSDataDir,AMSROOTDATADIR,AMSROOTPREFIX,AMSWD,DYLD_LIBRARY_PATH,HPCHIGHENERGYDATADIR,HPCHIGHENERGYRESULTDIR,HPCHIGHENERGYWORKDIR,HPCINTERMEDIATEDIR,HPCLOWENERGYDIR,JUAMSHIGHENERGYDATADIR,JUAMSINTERMEDIATEENERGYDIR,LD_LIBRARY_PATH,MY_ANALYSIS,PATH,PRODUCTION_SETTINGS,PYTHONPATH,ROOTSYS,ROOT_INCLUDE_PATH
export PMI_BARRIER_ROUNDS=5

srun Measuringtime --startTime 1305800000 --endTime 1620024849 --cutoffmode GEOMETRIC --degree 35 --safetyfactor 1.2
RET=${?}
echo
echo job 00000 exited with return code ${RET}

exit ${RET}
