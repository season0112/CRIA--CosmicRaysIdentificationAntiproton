#!/usr/bin/env python
#### above4.43: intermediate

import sys
import numpy as np
import ROOT
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad
from PyQt5 import QtCore, QtGui, QtWidgets
import argparse
import os
import binning

parser = argparse.ArgumentParser()
parser.add_argument('-C','--cluster', help='which cluster you choose')
arguments = parser.parse_args()

if arguments.cluster == "JUAMS":
    intermediateworkpath = os.getenv('JUAMSINTERMEDIATEENERGYDIR')
elif arguments.cluster == "HPC":
    intermediateworkpath = os.getenv('HPCINTERMEDIATEDIR')
    lowworkpath = os.getenv('HPCLOWENERGYDIR')

WIDTH  = 1000
HEIGHT = 1000
BORDER = 350
GAP    = 200
SIZE   = 3.5

SCALE = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
COLOR = (QtCore.Qt.magenta, QtCore.Qt.blue, QtCore.Qt.red, QtCore.Qt.black, QtCore.Qt.darkYellow, QtCore.Qt.darkGreen)
#PLOTRANGE = ["7.76_9.26" , "6.47_7.76" ,"5.37_6.47", "4.43_5.37", "2.97_3.64", "1.51_1.92"]
#PLOTRANGE = ["4.43_5.37", "2.97_3.64", "1.92_2.4"]
PLOTRANGE = [ "4.43_5.37", "2.97_3.64", "1.92_2.4", "1.51_1.92" ]

# method to convert a UTC timestamp into an angle (2011 - 2021)
def convertTime(utc):
    UTC_2011 = 1293840000.0
    UTC_2022 = 1640995200.0
    return (utc - UTC_2011) / (UTC_2022 - UTC_2011) * 2.0 * np.pi

def drawPlot(painter, arrays, averatio):
    painter.setRenderHint(QtGui.QPainter.Antialiasing)
    # background
    painter.setPen(QtCore.Qt.NoPen)
    painter.setBrush(QtCore.Qt.white)
    painter.drawRect(0, 0, WIDTH + BORDER, HEIGHT)
    # circle
    pen = QtGui.QPen(QtCore.Qt.black)
    pen.setWidth(5)
    painter.setPen(pen)
    RADIUS = (min(WIDTH, HEIGHT) - GAP) // 2
    painter.drawEllipse(WIDTH // 2 - RADIUS, HEIGHT // 2 - RADIUS, 2 * RADIUS, 2 * RADIUS)
    # find maximum ratio and derive radius scale
    maxRatio = max(arrays[0][1])  #arrays[0][1] is R_max ratio. first index is R, second index 0 is timestamp, 1 is ratio.
    RSCALE = RADIUS * 0.99 / maxRatio
    # draw average ratio
    for i in range (len(averatio)):
        pen2 = QtGui.QPen(COLOR[i])
        pen2.setWidth(5)
        painter.setPen(pen2)
        painter.setBrush(QtCore.Qt.white);
        #RADIUS2 = (min(WIDTH, HEIGHT) - GAP) // 2 // 1.5
        RADIUS2 = averatio[i] * RSCALE
        #painter.drawArc(WIDTH // 2.0 - RADIUS2, HEIGHT // 2.0 - RADIUS2, 2 * RADIUS2, 2 * RADIUS2, 700, 4850)
    # year lines
    pen.setWidth(3)
    pen.setStyle(QtCore.Qt.DashLine)
    painter.setPen(pen)
    for i in range(11):
        angle = i / 11.0 * 2.0 * np.pi
        painter.drawLine(WIDTH // 2, HEIGHT // 2, int(WIDTH / 2.0 + RADIUS * np.cos(angle)), int(HEIGHT / 2.0 + RADIUS * np.sin(angle)))  #Draws a line from (x1, y1) to (x2, y2).
    # year labels
    font = QtGui.QFont()
    font.setPointSize(20)
    font.setBold(True)
    painter.setFont(font)
    for i in range(11):
        year = 2011 + i
        angle = (i + 0.5) / 11.0 * 2.0 * np.pi
        x = WIDTH  / 2.0 + RADIUS * 1.1 * np.cos(angle)
        y = HEIGHT / 2.0 + RADIUS * 1.1 * np.sin(angle)
        rect = QtCore.QRect(int(x - 50), int(y - 15), 100, 30)
        painter.drawText(rect, QtCore.Qt.AlignCenter, str(year))
    # draw radius scale
    for val in (0, 0.5, 1.0):  # polar X axis
        x = WIDTH / 2.0 + val * RSCALE
        rect = QtCore.QRect(int(x - 30), HEIGHT // 2, 60, 30) #QRect(left, top, width, height)
        painter.drawText(rect, QtCore.Qt.AlignCenter, str(val))
    painter.drawText(QtCore.QRect(WIDTH // 2 + RADIUS + 20, HEIGHT // 2, BORDER, 30), int(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignLeft), 'Ratio(*10^4)')
    # draw ratio
    for i in range (len(arrays)):
        painter.setBrush(COLOR[i])
        for x,y,z in zip(arrays[i][0], arrays[i][1], arrays[i][2]):
            angle = convertTime(x)
            r     = y * RSCALE * SCALE[i]
            LINE  = z * RSCALE 
            painter.setPen(QtCore.Qt.NoPen)
            painter.drawEllipse(int(WIDTH // 2.0 + r * np.cos(angle)) - SIZE, int(HEIGHT // 2.0 + r * np.sin(angle)) - SIZE, 2 * SIZE, 2 * SIZE) #(int x, int y, int width, int height)，
            pen3 = QtGui.QPen(COLOR[i])
            painter.setPen(pen3)
            # painter.drawLine( int(WIDTH // 2.0 + r * np.cos(angle)) - LINE * np.cos(angle), int(HEIGHT // 2.0 + r * np.sin(angle)) - LINE * np.sin(angle), int(WIDTH // 2.0 + r * np.cos(angle)) + LINE * np.cos(angle), int(HEIGHT // 2.0 + r * np.sin(angle)) + LINE * np.sin(angle)  ); #(int x1, int y1, int x2, int y2)，
    # draw legend
    y = 20
    for i, arr in enumerate(arrays):
        painter.setBrush(COLOR[i])
        painter.setPen(COLOR[i])
        painter.drawEllipse(WIDTH, y + 10, 10, 10)
        text = '%s GV' % PLOTRANGE[i]
        if SCALE[i] > 1.0:
            text += ' * %.0f' % SCALE[i]
        painter.drawText(QtCore.QRect(WIDTH + 20, y, BORDER - 20, 30), int(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignLeft), text)
        y += 40
    
def main():
    #Load time dependent ratio
    arrays=list([])
    '''
    PLOTRANGE_intermediate = ["7.76_9.26" , "6.47_7.76" ,"5.37_6.47", "4.43_5.37", ] # "15.3_18", "13_15.3", "11_13", "9.26_11", "7.76_9.26", "6.47_7.76", "5.37_6.47", "4.43_5.37", "3.64_4.43", "2.97_3.64"
    for index in PLOTRANGE_intermediate:
        with open (intermediateworkpath + "/total/" + "results/fit_results_" + "0124" + "_" + "free" + "_" + index + "_3BartalRotation.txt") as f3b:
            lines3b=f3b.readlines()
        with open (intermediateworkpath + "/total/" + "results/fit_results_error_" + "0124" + "_" + "free" + "_" + index + "_3BartalRotation.txt") as e3b:
            error3b=e3b.readlines()
        B3result = np.array(list(map(lambda s: float(s.strip()), lines3b)))*10000
        E3result = np.array(list(map(lambda s: float(s.strip()), error3b)))*10000
        dates_unixsecond = binning.Bartals3Unixtime[0:B3result.shape[0]]
        arrays.append(np.array([list(dates_unixsecond), list(B3result), list(E3result)]))
    '''
    
    PLOTRANGE_low = ["4.43_5.37", "2.97_3.64", "1.92_2.4", "1.51_1.92" ]  # "4.43_5.37", "3.64_4.43", "2.97_3.64", "2.4_2.97", "1.92_2.4", "1.51_1.92",
    for index in PLOTRANGE_low:
        with open (lowworkpath + "/totalall/" + "results/fit_results_" + index + "_3BartalRotation.txt") as f3b:
            lines3b=f3b.readlines()
        with open (lowworkpath + "/totalall/" + "results/fit_results_error_" + index + "_3BartalRotation.txt") as e3b:
            error3b=e3b.readlines()
        B3result = np.array(list(map(lambda s: float(s.strip()), lines3b)))*10000
        E3result = np.array(list(map(lambda s: float(s.strip()), error3b)))*10000
        dates_unixsecond = binning.Bartals3Unixtime[0:B3result.shape[0]]
        arrays.append(np.array([list(dates_unixsecond), list(B3result), list(E3result)]))
    

    '''
    PLOTRANGE_low = ["4.43_5.37", "2.97_3.64", "1.92_2.4", "1.51_1.92" ]  # "4.43_5.37", "3.64_4.43", "2.97_3.64", "2.4_2.97", "1.92_2.4", "1.51_1.92",
    for index in PLOTRANGE_low:
        with open (lowworkpath + "/totalall/" + "results/fit_results_" + index + "_6months.txt") as f6m:
            lines6m=f6m.readlines()
        with open (lowworkpath + "/totalall/" + "results/fit_results_error_" + index + "_6months.txt") as e6m:
            error6m=e6m.readlines()
        M6result = np.array(list(map(lambda s: float(s.strip()), lines6m)))*10000
        EM6result = np.array(list(map(lambda s: float(s.strip()), error6m)))*10000
        dates_unixsecond = binning.Monthes6Unixtime[0:19] # 
        arrays.append(np.array([list(dates_unixsecond), list(M6result[0:19]), list(EM6result[0:19])]))
    '''

    # Load time averaged ratio
    averatio=list([])
    '''
    # Option1: Average Fit
    averaged_intermedaite = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/" + "intermediate_0124_free_pass7.8binmerge2.root")
    averagedratio_intermediate = averaged_intermedaite.Get("g_ratio")
    for i in [5,4,3,2]:
        averatio.append(averagedratio_intermediate.GetY()[i]*10000)
    averaged_low = TFile(lowworkpath + "/totalall/Time_Averaged_ratio_tof/Ratio_pass7.8.root")
    averagedratio_low = averaged_low.Get("ratio_tof_with_effective")
    for i in [3,0]:  
        averatio.append(averagedratio_low.GetY()[i]) 
    '''
    # Option2: Average value
    for i in range (4):
        averatio.append(np.array(np.mean(arrays[i][1])))

    # Create QApplication without X server
    app = QtWidgets.QApplication(['-platform', 'offscreen'])
    # create plot
    image = QtGui.QImage(WIDTH + BORDER, HEIGHT, QtGui.QImage.Format_RGB32)
    painter = QtGui.QPainter(image)
    drawPlot(painter, arrays, averatio)
    painter.end()
#    image.save(intermediateworkpath + "/total/results/plot/Polarplot/" + 'AntiprotonRatio_Polar.png', 'png')
    image.save(lowworkpath + "/totalall/results/plot/Polarplot/" + 'AntiprotonRatio_Polar.png', 'png')

if __name__ == '__main__':
    main()

