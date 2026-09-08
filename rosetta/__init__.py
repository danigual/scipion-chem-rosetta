# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import pwem
import os
import fnmatch
import pyworkflow.utils as pwutils
from pwem import Config as emConfig

from pyworkflow.utils import Environ
from .constants import *



_version_ = "3.12"
_logo = "rosetta_icon.png"
_references = ['LeaverFay2011']


ROSETTA_DIC = {'name': 'rosetta', 'version': '3.12', 'home': 'ROSETTA_HOME'}

# FRODOCK is a separate external tool (not Rosetta), used by the PROTAC-Model pipeline
# for the initial global protein-protein docking step. No 'version' key: we don't pin/
# validate a specific FRODOCK version, only that FRODOCK_HOME points somewhere real.
FRODOCK_DIC = {'name': 'frodock', 'home': 'FRODOCK_HOME'}

# The remaining four are used later, by filterPosesStep (PROTAC-Model's filter_frodock()):
# reduce/obabel/obenergy/prepare_receptor/prepare_ligand (ADFRsuite), vina (Vina), the
# voronota-voromqa binary (Voromqa), and the FCC clustering scripts (plain Python scripts,
# not a compiled binary, but still resolved the same way via an env-var home directory).
ADFRSUITE_DIC = {'name': 'adfrsuite', 'home': 'ADFRSUITE_HOME'}
VINA_DIC = {'name': 'vina', 'home': 'VINA_HOME'}
VOROMQA_DIC = {'name': 'voromqa', 'home': 'VOROMQA_HOME'}
FCC_DIC = {'name': 'fcc', 'home': 'FCC_HOME'}

class Plugin(pwem.Plugin):
    _homeVar = ROSETTA_DIC['home']
    _pathVars = [ROSETTA_DIC['home']]
    _supportedVersions = [ROSETTA_DIC['version']]


    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file. Set the Rosetta path on the computer
        """
        cls._defineVar(ROSETTA_DIC['home'], cls.getRosettaDir())

        # FRODOCK_HOME: unlike ROSETTA_HOME, there's no auto-detection (no fixed folder
        # name to search for under EM_ROOT). defaultValue=None means: if the user hasn't
        # set FRODOCK_HOME in scipion.conf or their shell environment, this stays None.
        cls._defineVar(FRODOCK_DIC['home'], None)

        # Same story for the four tools used by filterPosesStep: no auto-detection, the
        # user has to point each one at their local install (or a shared one on the CNB
        # machine) via scipion.conf or the shell environment.
        cls._defineVar(ADFRSUITE_DIC['home'], None)
        cls._defineVar(VINA_DIC['home'], None)
        cls._defineVar(VOROMQA_DIC['home'], None)
        cls._defineVar(FCC_DIC['home'], None)


    @classmethod
    def defineBinaries(cls, env):
        pass


    @classmethod  #  Test that
    def getEnviron(cls):
        """ Setup the environment variables needed to launch programs from Rosetta """
        environ = pwutils.Environ(os.environ)
        environ.update({
            'ROSETTA': cls.getHome(''),
            'PATH': cls.getHome(''),
        }, position=Environ.BEGIN)
        return environ


    @classmethod
    def getProgram(cls, progName, path=ROSETTA_BINARIES_PATH):
        """ Return the program binary that will be used. """
        return os.path.join(cls.getHome(),
                            path,
                            progName)

    @classmethod
    def _requireToolHome(cls, toolDic):
        """ Return toolDic['home']'s configured value, raising if the user never set it.
        Shared by the getXProgram() helpers below so the "not configured" error is
        consistent across tools instead of duplicated once per tool. """
        home = cls.getVar(toolDic['home'])
        if home is None:
            raise FileNotFoundError(
                '%s is not set. Point it to your %s installation '
                '(e.g. in scipion.conf or as a shell environment variable).'
                % (toolDic['home'], toolDic['name']))
        return home

    @classmethod
    def getFrodockProgram(cls, progName):
        """ Return the FRODOCK binary that will be used, trying the intel build first and
        falling back to the gcc build (FRODOCK ships both, e.g. frodockgrid/frodockgrid_gcc,
        and only one is guaranteed to work on a given machine). """
        home = cls._requireToolHome(FRODOCK_DIC)

        # The two candidate paths, same convention as FRODOCK's own <name>/<name>_gcc pair.
        intel = os.path.join(home, 'bin', progName)
        gcc = os.path.join(home, 'bin', '%s_gcc' % progName)
        if os.path.exists(intel):
            return intel
        elif os.path.exists(gcc):
            return gcc
        else:
            # Neither build is present: fail loudly with both checked paths, instead of the
            # original PROTAC-Model script's print() + sys.exit() (which would kill the whole
            # Scipion process, not just this step).
            raise FileNotFoundError(
                '%s not found under FRODOCK_HOME/bin (%s). Checked %s and %s.'
                % (progName, home, intel, gcc))

    @classmethod
    def getADFRSuiteProgram(cls, progName):
        """ Return an ADFRsuite binary (reduce, obabel, obenergy, prepare_receptor,
        prepare_ligand...). Unlike FRODOCK, ADFRsuite ships a single build: no intel/gcc
        fallback needed. """
        home = cls._requireToolHome(ADFRSUITE_DIC)
        path = os.path.join(home, 'bin', progName)
        if not os.path.exists(path):
            raise FileNotFoundError(
                '%s not found under ADFRSUITE_HOME/bin (%s).' % (progName, home))
        return path

    @classmethod
    def getVinaProgram(cls):
        """ Return the Vina binary. No progName parameter: VINA_HOME only ever provides
        this one program, unlike ADFRsuite/FRODOCK which bundle several. """
        home = cls._requireToolHome(VINA_DIC)
        path = os.path.join(home, 'bin', 'vina')
        if not os.path.exists(path):
            raise FileNotFoundError('vina not found under VINA_HOME/bin (%s).' % home)
        return path

    @classmethod
    def getVoromqaProgram(cls):
        """ Return the Voromqa binary (voronota-voromqa). Same single-program case as
        getVinaProgram. """
        home = cls._requireToolHome(VOROMQA_DIC)
        path = os.path.join(home, 'bin', 'voronota-voromqa')
        if not os.path.exists(path):
            raise FileNotFoundError('voronota-voromqa not found under VOROMQA_HOME/bin (%s).' % home)
        return path

    @classmethod
    def getFCCScript(cls, scriptName):
        """ Return the path to an FCC clustering script (make_contacts.py,
        calc_fcc_matrix.py, cluster_fcc.py, ppretty_clusters.py...). These are plain
        Python 2 scripts, not compiled binaries: unlike the other getXProgram() helpers,
        the caller must still prepend its own interpreter (e.g. 'python2')
        when building the command, this only resolves the script path. """
        home = cls._requireToolHome(FCC_DIC)
        path = os.path.join(home, 'scripts', scriptName)
        if not os.path.exists(path):
            raise FileNotFoundError('%s not found under FCC_HOME/scripts (%s).' % (scriptName, home))
        return path

    @classmethod
    def runProgram(cls, program, args=None, extraEnvDict=None, cwd=None):
        """ Internal shortcut function to launch an external program (Rosetta or, e.g.,
        FRODOCK). Not tool-specific: only builds the environment and launches the process. """
        env = cls.getEnviron()
        if extraEnvDict is not None:
            env.update(extraEnvDict)
        pwutils.runJob(None, program, args, env=env, cwd=cwd)

    @classmethod
    def runRosettaProgram(cls, program, args=None, extraEnvDict=None, cwd=None):
        """ Kept for backwards compatibility with existing call sites; delegates to the
        tool-agnostic runProgram. """
        cls.runProgram(program, args, extraEnvDict, cwd)

    @classmethod
    def validateInstallation(cls):
        """ Check if the installation of this protocol is correct.
             Returning an empty list means that the installation is correct
             and there are not errors. If some errors are found, a list with
             the error messages will be returned.
             """
        missingPaths = []
        if not os.path.exists(os.path.expanduser(cls.getVar(ROSETTA_DIC['home']))):
            missingPaths.append("Path of Rosetta does not exist (%s) : %s " % (ROSETTA_DIC['home'],
                                                                               cls.getVar(ROSETTA_DIC['home'])))
        return missingPaths



    # ---------------------------------- Utils functions  -----------------------
    @staticmethod
    def find(path, pattern):  # This function is analogous to linux find (case of folders)
        paths = []
        for root, dirs, files in os.walk(path):
            for name in dirs:
                if fnmatch.fnmatch(name, pattern):
                    paths.append(os.path.join(root, name))
        return paths


    @classmethod
    def getRosettaDir(cls, fn=""):
        fileList = cls.find(emConfig.EM_ROOT, "rosetta_bin_linux*")
        if len(fileList) == 0:
            return None
        else:
            if fn == "":
                return fileList[0]
            else:
                return os.path.join(fileList[0], fn)
