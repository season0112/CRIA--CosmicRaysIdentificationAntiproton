# Source this script to set up the ExampleAnalysis build that this script is part of.

drop_from_path()
{
   # Assert that we got enough arguments
   if test $# -ne 2 ; then
      echo "drop_from_path: needs 2 arguments"
      return 1
   fi

   p=$1
   drop=$2

   newpath=`echo $p | sed -e "s;:${drop}:;:;g" \
                          -e "s;:${drop};;g"   \
                          -e "s;${drop}:;;g"   \
                          -e "s;${drop};;g"`
}

if [ -n "${MY_ANALYSIS}" ] ; then
   old_analysis=${MY_ANALYSIS}
fi

if [ "x${BASH_ARGV[0]}" = "x" ]; then
    if [ ! -f Scripts/thisanalysis.sh ]; then
        echo ERROR: must "'cd where/analysis/is'" before calling "'source Scripts/thisanalysis.sh'" for this version of bash!
        MY_ANALYSIS=; export MY_ANALYSIS
        return 1
    fi
    MY_ANALYSIS="$PWD"; export MY_ANALYSIS
else
    # get param to "."
    thisanalysis=$(dirname ${BASH_ARGV[0]})
    MY_ANALYSIS=$(cd ${thisanalysis}/..;pwd); export MY_ANALYSIS
fi

if [ -n "${old_analysis}" ] ; then
   if [ -n "${PATH}" ]; then
      drop_from_path "$PATH" ${old_analysis}/bin
      PATH=$newpath
      drop_from_path "$PATH" ${old_analysis}/Scripts
      PATH=$newpath
   fi
   if [ -n "${LD_LIBRARY_PATH}" ]; then
      drop_from_path $LD_LIBRARY_PATH ${old_analysis}/lib
      LD_LIBRARY_PATH=$newpath
   fi
   if [ -n "${DYLD_LIBRARY_PATH}" ]; then
      drop_from_path $DYLD_LIBRARY_PATH ${old_analysis}/lib
      DYLD_LIBRARY_PATH=$newpath
   fi
   if [ -n "${SHLIB_PATH}" ]; then
      drop_from_path $SHLIB_PATH ${old_analysis}/lib
      SHLIB_PATH=$newpath
   fi
   if [ -n "${LIBPATH}" ]; then
      drop_from_path $LIBPATH ${old_analysis}/lib
      LIBPATH=$newpath
   fi
   if [ -n "${AC_COOKBOOK_MACRO_PATH}" ]; then
      drop_from_path $AC_COOKBOOK_MACRO_PATH ${old_analysis}
      AC_COOKBOOK_MACRO_PATH=$newpath
   fi
   if [ -n "${ACSOFT_JOB_SETTINGS}" ]; then
      drop_from_path $ACSOFT_JOB_SETTINGS ${old_analysis}/Configuration/analysis-job-settings.cfg
      ACSOFT_JOB_SETTINGS=$newpath
   fi
   if [ -n "${CLUSTERS_JOB_SETTINGS}" ]; then
      drop_from_path $CLUSTERS_JOB_SETTINGS ${old_analysis}/Configuration/analysis-job-settings.cfg
      CLUSTERS_JOB_SETTINGS=$newpath
   fi
   if [ -n "${ACSOFT_ADDITIONAL_LOGONS}" ]; then
      drop_from_path $ACSOFT_ADDITIONAL_LOGONS ${old_analysis}/rootlogon.C
      ACSOFT_ADDITIONAL_LOGONS=$newpath
   fi
fi

if [ -z "${PATH}" ]; then
   PATH=$MY_ANALYSIS/bin:$MY_ANALYSIS/Scripts; export PATH
else
   PATH=$MY_ANALYSIS/bin:$MY_ANALYSIS/Scripts:$PATH; export PATH
fi

if [ -z "${LD_LIBRARY_PATH}" ]; then
   LD_LIBRARY_PATH=$MY_ANALYSIS/lib; export LD_LIBRARY_PATH       # Linux, ELF HP-UX
else
   LD_LIBRARY_PATH=$MY_ANALYSIS/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
fi

if [ -z "${DYLD_LIBRARY_PATH}" ]; then
   DYLD_LIBRARY_PATH=$MY_ANALYSIS/lib; export DYLD_LIBRARY_PATH   # Mac OS X
else
   DYLD_LIBRARY_PATH=$MY_ANALYSIS/lib:$DYLD_LIBRARY_PATH; export DYLD_LIBRARY_PATH
fi

if [ -z "${SHLIB_PATH}" ]; then
   SHLIB_PATH=$MY_ANALYSIS/lib; export SHLIB_PATH                 # legacy HP-UX
else
   SHLIB_PATH=$MY_ANALYSIS/lib:$SHLIB_PATH; export SHLIB_PATH
fi

if [ -z "${LIBPATH}" ]; then
   LIBPATH=$MY_ANALYSIS/lib; export LIBPATH                       # AIX
else
   LIBPATH=$MY_ANALYSIS/lib:$LIBPATH; export LIBPATH
fi

if [ -z "${AC_COOKBOOK_MACRO_PATH}" ]; then
   AC_COOKBOOK_MACRO_PATH=$MY_ANALYSIS; export AC_COOKBOOK_MACRO_PATH
else
   AC_COOKBOOK_MACRO_PATH=$MY_ANALYSIS:$AC_COOKBOOK_MACRO_PATH; export AC_COOKBOOK_MACRO_PATH
fi

if [ -z "${ACSOFT_JOB_SETTINGS}" ]; then
   ACSOFT_JOB_SETTINGS=$MY_ANALYSIS/Configuration/analysis-job-settings.cfg; export ACSOFT_JOB_SETTINGS
else
   ACSOFT_JOB_SETTINGS=$ACSOFT_JOB_SETTINGS:$MY_ANALYSIS/Configuration/analysis-job-settings.cfg; export ACSOFT_JOB_SETTINGS
fi

if [ -z "${CLUSTERS_JOB_SETTINGS}" ]; then
   CLUSTERS_JOB_SETTINGS=$MY_ANALYSIS/Configuration/analysis-job-settings.cfg; export CLUSTERS_JOB_SETTINGS
else
   CLUSTERS_JOB_SETTINGS=$CLUSTERS_JOB_SETTINGS:$MY_ANALYSIS/Configuration/analysis-job-settings.cfg; export CLUSTERS_JOB_SETTINGS
fi

if [ -z "${ACSOFT_ADDITIONAL_LOGONS}" ]; then
   ACSOFT_ADDITIONAL_LOGONS=$MY_ANALYSIS/rootlogon.C; export ACSOFT_ADDITIONAL_LOGONS
else
   ACSOFT_ADDITIONAL_LOGONS=$ACSOFT_ADDITIONAL_LOGONS:$MY_ANALYSIS/rootlogon.C; export ACSOFT_ADDITIONAL_LOGONS
fi

unset old_analysis
unset thisanalysis
