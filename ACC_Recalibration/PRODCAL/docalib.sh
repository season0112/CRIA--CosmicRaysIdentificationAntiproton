#!/bin/sh

root -b -q makeChargeCalibrations.C+\(1\) &
root -b -q makeChargeCalibrations.C+\(2\) &
root -b -q makeChargeCalibrations.C+\(3\) &
root -b -q makeChargeCalibrations.C+\(4\) &
root -b -q makeChargeCalibrations.C+\(5\) &
root -b -q makeChargeCalibrations.C+\(6\) &
root -b -q makeChargeCalibrations.C+\(7\) &
root -b -q makeChargeCalibrations.C+\(8\) &
wait

root -b -q makeZetaCalibrations.C+\(1\) \
&& root -b -q makeZetaCalibrations.C+\(2\) \
&& root -b -q makeZetaCalibrations.C+\(3\) \
&& root -b -q makeZetaCalibrations.C+\(4\) \
&& root -b -q makeZetaCalibrations.C+\(5\) \
&& root -b -q makeZetaCalibrations.C+\(6\) \
&& root -b -q makeZetaCalibrations.C+\(7\) \
&& root -b -q makeZetaCalibrations.C+\(8\) \
&& root -b -q makeTimeCalibrations.C+\(1\) \
&& root -b -q makeTimeCalibrations.C+\(2\) \
&& root -b -q makeTimeCalibrations.C+\(3\) \
&& root -b -q makeTimeCalibrations.C+\(4\) \
&& root -b -q makeTimeCalibrations.C+\(5\) \
&& root -b -q makeTimeCalibrations.C+\(6\) \
&& root -b -q makeTimeCalibrations.C+\(7\) \
&& root -b -q makeTimeCalibrations.C+\(8\)

hadd AntiPG_calib_2.root OUTCAL/*.root
exit



