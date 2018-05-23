# Module for guiding 1D Wavelength Calibration
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from importlib import reload


from pypit import msgs
from pypit import ardebug as debugger
from pypit import masterframe
from pypit import ararc
from pypit import ararclines
from pypit.core import arsort

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

frametype = 'wv_calib'

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
default_settings = dict(calibrate={'nfitpix': 5}
                        )
#settings_spect[dnum]['saturation']*settings_spect[dnum]['nonlinear'])  -- For satmask (echelle)

#  See save_master() for the data model for output


class WaveCalib(masterframe.MasterFrame):
    """Class to guide slit/order tracing

    Parameters
    ----------
    mstrace : ndarray
      Trace image
    pixlocn : ndarray
      Pixel location array
    binbpx : ndarray, optional
      Bad pixel mask
      If not provided, a dummy array with no masking is generated
    settings : dict, optional
      Settings for trace slits
    det : int, optional
      Detector number
    ednum : int, optional
      Edge number used for indexing

    Attributes
    ----------
    frametype : str
      Hard-coded to 'wv_calib'

    steps : list
      List of the processing steps performed
    """
    def __init__(self, msarc, spectrograph=None, settings=None, det=None, setup=None, fitstbl=None, sci_ID=None):

        # Required parameters (but can be None)
        self.msarc = msarc

        # Optional parameters
        self.det = det
        self.fitstbl = fitstbl
        self.setup = setup
        self.sci_ID = sci_ID
        if settings is None:
            self.settings = default_settings.copy()
        else:
            self.settings = settings.copy()
            if 'calibrate' not in self.settings.keys():
                self.settings.update(default_settings)
        self.spectrograph = spectrograph

        # Attributes
        self.frametype = frametype
        self.steps = []

        # Main outputs
        self.wv_calib = {}

        # Key Internals

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)

    def _build_wv_calib(self, method, skip_QA=False):
        reload(ararclines)
        # Loop
        self.wv_calib = {}
        ok_mask = np.where(self.maskslit == 0)[0]
        for slit in ok_mask:
            ###############
            # Extract arc and identify lines
            #if settings.argflag['arc']['calibrate']['method'] == 'simple':
            #elif settings.argflag['arc']['calibrate']['method'] == 'arclines':
            if method == 'simple':
                iwv_calib = ararc.simple_calib(self.det, self.msarc, self.arcparam,
                                               censpec=self.arccen[:, slit], slit=slit)
            elif method == 'arclines':
                iwv_calib = ararc.calib_with_arclines(slit, self.arcparam, self.arccen[:, slit])
            self.wv_calib[str(slit)] = iwv_calib.copy()
            # QA
            if not skip_QA:
                ararc.arc_fit_qa(self.setup, iwv_calib, slit)
        # Return
        return self.wv_calib

    def _extract_arcs(self, lordloc, rordloc, pixlocn):
        self.arccen, self.maskslit, _ = ararc.get_censpec(lordloc, rordloc, pixlocn,
                                                          self.msarc, self.det, self.settings,
                                                          gen_satmask=False)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.arccen, self.maskslit

    def _load_arcparam(self, calibrate_lamps=None):
        """

        Parameters
        ----------
        calibrate_lamps : str, optional
           List of lamps used

        Returns
        -------

        """
        # Setup arc parameters (e.g. linelist)
        arc_idx = arsort.ftype_indices(self.fitstbl, 'arc', self.sci_ID)
        self.arcparam = ararc.setup_param(self.spectrograph, self.msarc.shape,
                                          self.fitstbl, arc_idx[0],
                                          calibrate_lamps=calibrate_lamps)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.arcparam

    def master(self, lordloc, rordloc, pixlocn, method='arclines', nonlinear=None):
        """ Main driver for wavelength calibration

        Parameters
        ----------

        Returns
        -------
        """
        # Attempt to load the Master Frame
        self.wv_calib, _, _ = self.load_master_frame(self, "wv_calib")
        if self.wv_calib is None:

            ###############
            # Extract an arc down each slit
            #   The settings here are settings.spect (saturation and nonlinear)
            _, _ = self._extract_arcs(lordloc, rordloc, pixlocn)

            # Load the arcparam
            _ = self._load_arcparam()

            # Fill up the calibrations and generate QA
            self.wv_calib = self._build_wv_calib(method)
            self.wv_calib['steps'] = self.steps
            # Save to Masters
            self.save_master(self.wv_calib)
        # Finish
        return self.wv_calib

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt



