#! /usr/bin/env python3

import numpy as np
import xml.etree.ElementTree as ET

class FluxData:
    """
    Simple class for cosmic-ray flux data, consisting of data points
    (energy, flux), and the corresponding statistical, systematic, and
    total uncertainties.
    """
    def __init__(self):
        self.E = None
        self.rigidity_min = None
        self.rigidity_max = None
        self.kinetic_energy_min = None
        self.kinetic_energy_max = None
        self.flux = None
        self.fluxratio = None
        self.err_stat = None
        self.err_syst = None
        self.err_total = None
        self.err_total_low = None
        self.err_total_high = None
        self.E_unit = 'void'
        self.flux_unit = 'void'
        
    def __str__(self):
        return '\n'.join([f'E= {self.E[i]:.4g} flux= {self.flux[i]:.4g}' for i,_ in enumerate(self.E)])

class ASIReader_E:
    '''
    Parser for XML files downloaded from ASI cosmic-ray database.
    '''
    def __init__(self):
        pass
    
    def read(self, inputfile):
        
        # initialize XML parser
        root = ET.parse(inputfile).getroot()

        # lists for data points and error bars
        v_E = []
        v_kinetic_energy_min = []
        v_kinetic_energy_max = []
        v_flux = []
        v_fluxratio = []
        v_err_stat = []
        v_err_syst = []

        result = FluxData()

        # read header information
        header = root.find('HEADER')
        mission = header.findtext('MISSION')
        units = header.find('UNITS')
        result.E_unit = units.findtext('kinetic_energy')
        result.flux_unit = units.findtext('flux')

        # extract individual data points
        for d in root.findall('DATA'):

            # X-axis
            E = float(d.findtext('kinetic_energy', default=0.0))
            E_min = float(d.findtext('kinetic_energy_min', default=0.0))
            E_max = float(d.findtext('kinetic_energy_max', default=0.0))

            # Y-axis
            #flux = float(d.findtext('flux'))
            fluxratio = float(d.findtext('fluxratio'))
            fluxratio_stat_down = float(d.findtext('fluxratio_statistical_error_low'))
            fluxratio_stat_up = float(d.findtext('fluxratio_statistical_error_high'))
            fluxratio_syst_down = float(d.findtext('fluxratio_systematical_error_low', default=0.0))
            fluxratio_syst_up = float(d.findtext('fluxratio_systematical_error_high', default=0.0))

            v_E.append(E)
            v_kinetic_energy_min.append(E_min)
            v_kinetic_energy_max.append(E_max)
            #v_flux.append(flux)
            v_fluxratio.append(fluxratio)
            v_err_stat.append(0.5 * (np.abs(fluxratio_stat_down) + np.abs(fluxratio_stat_up)))
            v_err_syst.append(0.5 * (np.abs(fluxratio_syst_down) + np.abs(fluxratio_syst_up)))

        # convert to numpy arrays and return result as FluxData class
        result.E = np.array(v_E)
        result.kinetic_energy_min = np.array(v_kinetic_energy_min)
        result.kinetic_energy_max = np.array(v_kinetic_energy_max)
        #result.flux = np.array(v_flux)
        result.fluxratio = np.array(v_fluxratio)
        result.err_stat = np.array(v_err_stat)
        result.err_syst = np.array(v_err_syst)
        result.err_total = np.sqrt(result.err_stat**2 + result.err_syst**2)
        return result


class ASIReader_R:
    '''
    Parser for XML files downloaded from ASI cosmic-ray database.
    '''
    def __init__(self):
        pass

    def readflux(self, inputfile):

        # initialize XML parser
        root = ET.parse(inputfile).getroot()

        # lists for data points and error bars
        v_E = []
        v_rigidity_min = []
        v_rigidity_max = []
        v_flux = []
        v_err_stat = []
        v_err_syst = []

        result = FluxData()

        # read header information
        header = root.find('HEADER')
        mission = header.findtext('MISSION')
        units = header.find('UNITS')
        result.E_unit = units.findtext('kinetic_energy')
        result.flux_unit = units.findtext('flux')

        # extract individual data points
        for d in root.findall('DATA'):

            # X-axis
            R_min = float(d.findtext('rigidity_min', default=0.0))
            R_max = float(d.findtext('rigidity_max', default=0.0))

            # Y-axis
            flux = float(d.findtext('flux'))
            flux_stat_down = float(d.findtext('flux_statistical_error_low'))
            flux_stat_up = float(d.findtext('flux_statistical_error_high'))
            flux_syst_down = float(d.findtext('flux_systematical_error_low', default=0.0))
            flux_syst_up = float(d.findtext('flux_systematical_error_high', default=0.0))

            v_rigidity_min.append(R_min)
            v_rigidity_max.append(R_max)
            v_flux.append(flux)
            v_err_stat.append(0.5 * (np.abs(flux_stat_down) + np.abs(flux_stat_up)))
            v_err_syst.append(0.5 * (np.abs(flux_syst_down) + np.abs(flux_syst_up)))

        # convert to numpy arrays and return result as FluxData class
        result.E = np.array(v_E)
        result.rigidity_min = np.array(v_rigidity_min)
        result.rigidity_max = np.array(v_rigidity_max)
        result.flux = np.array(v_flux)
        result.err_stat = np.array(v_err_stat)
        result.err_syst = np.array(v_err_syst)
        result.err_total = np.sqrt(result.err_stat**2 + result.err_syst**2)
        return result

    def readfluxratio(self, inputfile):

        # initialize XML parser
        root = ET.parse(inputfile).getroot()

        # lists for data points and error bars
        v_E = []
        v_rigidity_min = []
        v_rigidity_max = []
        v_fluxratio = []
        v_err_stat = []
        v_err_syst = []

        result = FluxData()

        # read header information
        header = root.find('HEADER')
        mission = header.findtext('MISSION')
        units = header.find('UNITS')
        result.E_unit = units.findtext('kinetic_energy')
        result.flux_unit = units.findtext('flux')

        # extract individual data points
        for d in root.findall('DATA'):

            # X-axis
            R_min = float(d.findtext('rigidity_min', default=0.0))
            R_max = float(d.findtext('rigidity_max', default=0.0))

            # Y-axis
            fluxratio = float(d.findtext('fluxratio'))
            fluxratio_stat_down = float(d.findtext('fluxratio_statistical_error_low'))
            fluxratio_stat_up = float(d.findtext('fluxratio_statistical_error_high'))
            fluxratio_syst_down = float(d.findtext('fluxratio_systematical_error_low', default=0.0))
            fluxratio_syst_up = float(d.findtext('fluxratio_systematical_error_high', default=0.0))

            v_rigidity_min.append(R_min)
            v_rigidity_max.append(R_max)
            v_fluxratio.append(fluxratio)
            v_err_stat.append(0.5 * (np.abs(fluxratio_stat_down) + np.abs(fluxratio_stat_up)))
            v_err_syst.append(0.5 * (np.abs(fluxratio_syst_down) + np.abs(fluxratio_syst_up)))

        # convert to numpy arrays and return result as FluxData class
        result.E = np.array(v_E)
        result.rigidity_min = np.array(v_rigidity_min)
        result.rigidity_max = np.array(v_rigidity_max)
        result.fluxratio = np.array(v_fluxratio)
        result.err_stat = np.array(v_err_stat)
        result.err_syst = np.array(v_err_syst)
        result.err_total = np.sqrt(result.err_stat**2 + result.err_syst**2)
        return result


class ASIReader_R_PAMELAData:
    '''
    Parser for XML files downloaded from ASI cosmic-ray database.
    '''
    def __init__(self):
        pass

    def readfluxratio(self, inputfile):

        # initialize XML parser
        root = ET.parse(inputfile).getroot()

        # lists for data points and error bars
        v_E = []
        v_rigidity_min = []
        v_rigidity_max = []
        v_fluxratio = []
        v_toterr_low = []
        v_toterr_high = []

        result = FluxData()

        # read header information
        header = root.find('HEADER')
        mission = header.findtext('MISSION')
        units = header.find('UNITS')
        result.E_unit = units.findtext('kinetic_energy')
        result.flux_unit = units.findtext('flux')

        # extract individual data points
        for d in root.findall('DATA'):

            # X-axis
            R_min = float(d.findtext('rigidity_min', default=0.0))
            R_max = float(d.findtext('rigidity_max', default=0.0))

            # Y-axis
            fluxratio = float(d.findtext('fluxratio'))
            fluxratio_total_error_low = float(d.findtext('fluxratio_total_error_low'))
            fluxratio_total_error_high = float(d.findtext('fluxratio_total_error_high'))

            v_rigidity_min.append(R_min)
            v_rigidity_max.append(R_max)
            v_fluxratio.append(fluxratio)
            v_toterr_low.append(fluxratio_total_error_low)
            v_toterr_high.append(fluxratio_total_error_low)

        # convert to numpy arrays and return result as FluxData class
        result.E = np.array(v_E)
        result.rigidity_min = np.array(v_rigidity_min)
        result.rigidity_max = np.array(v_rigidity_max)
        result.fluxratio = np.array(v_fluxratio)
        result.err_total_low = np.array(v_toterr_low)
        result.err_total_high = np.array(v_toterr_high)
        return result








def LaffertyWyatt(E1,E2,gamma=2.7):
    denom = E2**(1-gamma)-E1**(1.0-gamma)
    num = (E2-E1)*(1.0-gamma)
    return (num / denom)**(1.0/gamma)


