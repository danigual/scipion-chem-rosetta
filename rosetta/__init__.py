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
from scipion.install.funcs import InstallHelper

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
# voronota-voromqa binary (Voromqa), and the FCC clustering scripts. Unlike FRODOCK/Rosetta
# (license-gated downloads, see below), none of these four require a license click-through,
# so they get real defineBinaries() support via InstallHelper - 'version' is needed now
# (InstallHelper/getEnvName use it to name the conda env / package folder).
ADFRSUITE_DIC = {'name': 'adfrsuite', 'version': '1.0', 'home': 'ADFRSUITE_HOME'}
VINA_DIC = {'name': 'vina', 'version': '1.2.7', 'home': 'VINA_HOME'}
VOROMQA_DIC = {'name': 'voromqa', 'version': '1.29.4816', 'home': 'VOROMQA_HOME'}
FCC_DIC = {'name': 'fcc', 'version': 'latest', 'home': 'FCC_HOME'}

# PROTAC-Model itself (gaoqiweng/PROTAC-Model): its own driving code (main.py, utils/*) is
# called directly by the protocol instead of being reimplemented here - see
# protocol_protac_model.py and rosetta/scripts/run_*.py. PROTAC_MODEL_HOME just needs a
# checkout of that repo (main.py + utils/), nothing to build. PROTAC_MODEL_PYTHON_HOME is
# separate: that repo's code is genuine Python 2 (print statements, dict.has_key(), etc.),
# so it needs its own Python 2.7 conda env (with rdkit for py2 installed) to run under -
# never Scipion's own scipion3 env, same isolation reasoning as pwchem's RDKIT_DIC/
# OPENBABEL_DIC. Both are public/open source (no license gate), so both get InstallHelper
# support too.
PROTAC_MODEL_DIC = {'name': 'protac-model', 'version': 'latest', 'home': 'PROTAC_MODEL_HOME'}
PROTAC_MODEL_PYTHON_DIC = {'name': 'protac-model-python', 'version': '2.7',
                           'home': 'PROTAC_MODEL_PYTHON_HOME'}

class Plugin(pwem.Plugin):
    _homeVar = ROSETTA_DIC['home']
    _pathVars = [ROSETTA_DIC['home']]
    _supportedVersions = [ROSETTA_DIC['version']]


    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file. Set the Rosetta path on the computer
        """
        cls._defineVar(ROSETTA_DIC['home'], cls.getRosettaDir())

        # FRODOCK_HOME: like ROSETTA_HOME, this one stays manual (license click-through
        # required on chaconlab.org before download) - no auto-detection, defaultValue
        # =None means: if the user hasn't set FRODOCK_HOME in scipion.conf or their shell
        # environment, this stays None until they point it at their own download.
        cls._defineVar(FRODOCK_DIC['home'], None)

        # The five tools below have no license gate, so their home is wherever
        # defineBinaries()/InstallHelper installs them (software/em/<name>-<version>/...)
        # instead of a manually-set path - _defineEmVar (not _defineVar) is what wires a
        # variable to that Scipion-managed package location.
        cls._defineEmVar(ADFRSUITE_DIC['home'], cls.getEnvName(ADFRSUITE_DIC))
        cls._defineEmVar(VINA_DIC['home'], cls.getEnvName(VINA_DIC))
        cls._defineEmVar(VOROMQA_DIC['home'], cls.getEnvName(VOROMQA_DIC))
        cls._defineEmVar(FCC_DIC['home'], cls.getEnvName(FCC_DIC))
        cls._defineEmVar(PROTAC_MODEL_DIC['home'], cls.getEnvName(PROTAC_MODEL_DIC))
        cls._defineEmVar(PROTAC_MODEL_PYTHON_DIC['home'], cls.getEnvName(PROTAC_MODEL_PYTHON_DIC))


    @classmethod
    def defineBinaries(cls, env):
        # Rosetta and FRODOCK stay out of here on purpose: both require accepting a
        # license (personal academic login for Rosetta, a click-through agreement for
        # FRODOCK) that can't be scripted, so they're never auto-installed. The
        # remaining five have no such gate.
        cls.addADFRSuitePackage(env)
        cls.addVinaPackage(env)
        cls.addVoromqaPackage(env)
        cls.addFCCPackage(env)
        cls.addProtacModelPackage(env)
        cls.addProtacModelPythonPackage(env)

    # ---------------------------- Conda env naming helpers -----------------------
    # Mirrors pwchemPlugin.getEnvName()/getEnvActivationCommand() (scipion-chem/pwchem):
    # kept as our own small copies instead of switching rosetta.Plugin's base class from
    # pwem.Plugin to pwchemPlugin, since getCondaActivationCmd() (what these build on) is
    # already available from pwem.Plugin's own parent, pyworkflow.plugin.Plugin.
    @classmethod
    def getEnvName(cls, packageDictionary):
        """ Name of the conda env / package folder for a given tool dict, e.g. 'vina-1.2.5'. """
        return f"{packageDictionary['name']}-{packageDictionary['version']}"

    @classmethod
    def getEnvActivationCommand(cls, packageDictionary, condaHook=True):
        """ Shell snippet that activates the conda env for a given tool dict. """
        hook = cls.getCondaActivationCmd() if condaHook else ''
        return f'{hook}conda activate {cls.getEnvName(packageDictionary)}'

    # ---------------------------- Package installers (InstallHelper) -------------
    @classmethod
    def addADFRSuitePackage(cls, env, default=True):
        """ Downloads and unpacks ADFRsuite (no license click-through found on
        ccsb.scripps.edu for the non-commercial installer, unlike FRODOCK).
        URL and install.sh flags confirmed against ccsb.scripps.edu/adfr/downloads/'s own
        "INSTALLING FROM TARBALL" instructions (2026-09-10): -d is the destination folder,
        -c 0/1 picks .pyc/.pyo compilation - no interactive prompt is documented, so no
        'echo "Y" |' is needed.
        Known issue (same page): on Linux with an older GCC, _openbabel.so can fail with
        "GLIBCXX_3.4.15' not found" - fix is renaming <installDir>/lib/libstdc++.so.6.orig
        to libstdc++.so.6. """
        installer = InstallHelper(ADFRSUITE_DIC['name'], packageHome=cls.getVar(ADFRSUITE_DIC['home']),
                                  packageVersion=ADFRSUITE_DIC['version'])
        installer.addCommand(
            'wget -q https://ccsb.scripps.edu/adfr/download/1038/ -O adfrsuite.tar.gz && '
            'tar -xzf adfrsuite.tar.gz && '
            './ADFRsuite*/install.sh -d . -c 0',
            targetName=f"{ADFRSUITE_DIC['name']}_installed")
        installer.addPackage(env, dependencies=['wget', 'tar'], default=default)

    @classmethod
    def addVinaPackage(cls, env, default=True):
        """ Installs the Vina CLI binary via conda-forge - PROTAC-Model shells out to
        $VINA/bin/vina, so we need the compiled binary, not just the 'vina' PyPI package
        (Python bindings only). """
        installer = InstallHelper(VINA_DIC['name'], packageHome=cls.getVar(VINA_DIC['home']),
                                  packageVersion=VINA_DIC['version'])
        installer.getCondaEnvCommand(
            binaryName=VINA_DIC['name'], binaryVersion=VINA_DIC['version'], pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(VINA_DIC)} && conda install -y -c conda-forge vina={VINA_DIC['version']}",
            targetName=f"{VINA_DIC['name']}_installed"
        ).addCommand(
            # getCondaEnvCommand installs under conda's own envs dir, not packageHome -
            # without this symlink getVinaProgram() would find an empty folder. Same fix
            # pwchem uses for MGLTools (pwchem/__init__.py, addMGLToolsPackage).
            f"{cls.getEnvActivationCommand(VINA_DIC)} && rm -rf {cls.getVar(VINA_DIC['home'])} && "
            f"ln -s $CONDA_PREFIX {cls.getVar(VINA_DIC['home'])}",
            targetName=f"{VINA_DIC['name']}_symlinked")
        installer.addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def addVoromqaPackage(cls, env, default=True):
        """ Installs Voromqa (voronota) via bioconda. """
        installer = InstallHelper(VOROMQA_DIC['name'], packageHome=cls.getVar(VOROMQA_DIC['home']),
                                  packageVersion=VOROMQA_DIC['version'])
        installer.getCondaEnvCommand(
            binaryName=VOROMQA_DIC['name'], binaryVersion=VOROMQA_DIC['version'], pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(VOROMQA_DIC)} && conda install -y -c bioconda voronota={VOROMQA_DIC['version']}",
            targetName=f"{VOROMQA_DIC['name']}_installed"
        ).addCommand(
            # Same symlink fix as addVinaPackage above.
            f"{cls.getEnvActivationCommand(VOROMQA_DIC)} && rm -rf {cls.getVar(VOROMQA_DIC['home'])} && "
            f"ln -s $CONDA_PREFIX {cls.getVar(VOROMQA_DIC['home'])}",
            targetName=f"{VOROMQA_DIC['name']}_symlinked")
        installer.addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def addFCCPackage(cls, env, default=True):
        """ Clones haddocking/fcc (no PyPI package) and compiles its C/C++ contact
        programs with make. Its Python scripts don't get their own conda env: they're
        only ever invoked from inside PROTAC-Model's own code (pre.cluster()), which
        already runs under PROTAC_MODEL_PYTHON_HOME and passes that interpreter down. """
        installer = InstallHelper(FCC_DIC['name'], packageHome=cls.getVar(FCC_DIC['home']),
                                  packageVersion=FCC_DIC['version'])
        # binaryFolderName=FCC_DIC['name'] (not '.'): getCloneCommand's 'cd <packageHome>
        # && git clone <url> .' only works if packageHome already exists and is empty -
        # pwchem always clones into a named subfolder instead (see e.g. addShapeItPackage).
        installer.getCloneCommand(
            'https://github.com/haddocking/fcc.git', binaryFolderName=FCC_DIC['name'],
            targeName=f"{FCC_DIC['name']}_cloned"
        # XXX FCC's README warns the Makefile may need manual edits - untested.
        ).addCommand(f"cd {FCC_DIC['name']}/src && make", targetName=f"{FCC_DIC['name']}_built")
        installer.addPackage(env, dependencies=['git', 'make', 'gcc'], default=default)

    @classmethod
    def addProtacModelPackage(cls, env, default=True):
        """ Clones gaoqiweng/PROTAC-Model itself - public repo, no build step, its code
        is called directly by the protocol/rosetta/scripts/run_*.py drivers. """
        installer = InstallHelper(PROTAC_MODEL_DIC['name'], packageHome=cls.getVar(PROTAC_MODEL_DIC['home']),
                                  packageVersion=PROTAC_MODEL_DIC['version'])
        # Same binaryFolderName reasoning as addFCCPackage above (named subfolder, not
        # '.') - getProtacModelScript()/getProtacModelPython() below account for it.
        installer.getCloneCommand(
            'https://github.com/gaoqiweng/PROTAC-Model.git', binaryFolderName=PROTAC_MODEL_DIC['name'],
            targeName=f"{PROTAC_MODEL_DIC['name']}_cloned")
        installer.addPackage(env, dependencies=['git'], default=default)

    @classmethod
    def addProtacModelPythonPackage(cls, env, default=True):
        """ Dedicated Python 2.7 conda env (with RDKit) to run PROTAC-Model's own code -
        never scipion3's own env, and not the local 'protac-model' conda env some
        machines may already have lying around (that one is Python 3.10, incompatible
        with this genuinely-Python-2 codebase).
        RDKit's own channel stopped publishing py2.7 builds after 2016.03.3 (checked
        anaconda.org/rdkit/rdkit's full file list, 2026-09-10) - pinned explicitly below.
        XXX: that build is tied to a specific old numpy ("np111py27"); unverified whether
        conda's solver picks a compatible numpy on its own or needs one pinned too -
        check on first real install. """
        installer = InstallHelper(PROTAC_MODEL_PYTHON_DIC['name'],
                                  packageHome=cls.getVar(PROTAC_MODEL_PYTHON_DIC['home']),
                                  packageVersion=PROTAC_MODEL_PYTHON_DIC['version'])
        installer.getCondaEnvCommand(
            binaryName=PROTAC_MODEL_PYTHON_DIC['name'],
            binaryVersion=PROTAC_MODEL_PYTHON_DIC['version'], pythonVersion='2.7'
        ).addCommand(
            f"{cls.getEnvActivationCommand(PROTAC_MODEL_PYTHON_DIC)} && conda install -y -c rdkit rdkit=2016.03.3",
            targetName=f"{PROTAC_MODEL_PYTHON_DIC['name']}_installed"
        ).addCommand(
            # Same symlink fix as addVinaPackage/addVoromqaPackage below.
            f"{cls.getEnvActivationCommand(PROTAC_MODEL_PYTHON_DIC)} && "
            f"rm -rf {cls.getVar(PROTAC_MODEL_PYTHON_DIC['home'])} && "
            f"ln -s $CONDA_PREFIX {cls.getVar(PROTAC_MODEL_PYTHON_DIC['home'])}",
            targetName=f"{PROTAC_MODEL_PYTHON_DIC['name']}_symlinked")
        installer.addPackage(env, dependencies=['conda'], default=default)


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
                f"{toolDic['home']} is not set. Point it to your {toolDic['name']} "
                "installation (e.g. in scipion.conf or as a shell environment variable).")
        return home

    @classmethod
    def getFrodockProgram(cls, progName):
        """ Return the FRODOCK binary that will be used, trying the intel build first and
        falling back to the gcc build (FRODOCK ships both, e.g. frodockgrid/frodockgrid_gcc,
        and only one is guaranteed to work on a given machine). """
        home = cls._requireToolHome(FRODOCK_DIC)

        # The two candidate paths, same convention as FRODOCK's own <name>/<name>_gcc pair.
        intel = os.path.join(home, 'bin', progName)
        gcc = os.path.join(home, 'bin', f'{progName}_gcc')
        if os.path.exists(intel):
            return intel
        elif os.path.exists(gcc):
            return gcc
        else:
            # Fail loudly with both checked paths, instead of the original PROTAC-Model
            # script's print() + sys.exit() (which would kill the whole Scipion process).
            raise FileNotFoundError(
                f'{progName} not found under FRODOCK_HOME/bin ({home}). Checked {intel} and {gcc}.')

    @classmethod
    def getADFRSuiteProgram(cls, progName):
        """ Return an ADFRsuite binary (reduce, obabel, obenergy, prepare_receptor,
        prepare_ligand...). Unlike FRODOCK, ADFRsuite ships a single build: no intel/gcc
        fallback needed. """
        home = cls._requireToolHome(ADFRSUITE_DIC)
        path = os.path.join(home, 'bin', progName)
        if not os.path.exists(path):
            raise FileNotFoundError(f'{progName} not found under ADFRSUITE_HOME/bin ({home}).')
        return path

    @classmethod
    def getVinaProgram(cls):
        """ Return the Vina binary. No progName parameter: VINA_HOME only ever provides
        this one program, unlike ADFRsuite/FRODOCK which bundle several. """
        home = cls._requireToolHome(VINA_DIC)
        path = os.path.join(home, 'bin', 'vina')
        if not os.path.exists(path):
            raise FileNotFoundError(f'vina not found under VINA_HOME/bin ({home}).')
        return path

    @classmethod
    def getVoromqaProgram(cls):
        """ Return the Voromqa binary (voronota-voromqa). Same single-program case as
        getVinaProgram. """
        home = cls._requireToolHome(VOROMQA_DIC)
        path = os.path.join(home, 'bin', 'voronota-voromqa')
        if not os.path.exists(path):
            raise FileNotFoundError(f'voronota-voromqa not found under VOROMQA_HOME/bin ({home}).')
        return path

    @classmethod
    def getFCCScript(cls, scriptName):
        """ Path to an FCC clustering script (make_contacts.py, calc_fcc_matrix.py,
        cluster_fcc.py, ppretty_clusters.py...). Plain Python 2 scripts, not compiled
        binaries: the caller must still prepend its own interpreter. """
        # FCC_HOME is the parent InstallHelper cloned into; the repo itself lives one
        # level down, in a named subfolder (see addFCCPackage's binaryFolderName).
        home = os.path.join(cls._requireToolHome(FCC_DIC), FCC_DIC['name'])
        path = os.path.join(home, 'scripts', scriptName)
        if not os.path.exists(path):
            raise FileNotFoundError(f"{scriptName} not found under FCC_HOME/{FCC_DIC['name']}/scripts ({home}).")
        return path

    @classmethod
    def getProtacModelScript(cls, scriptName=''):
        """ Return a path inside the PROTAC-Model checkout (gaoqiweng/PROTAC-Model), e.g.
        getProtacModelScript() for the repo root (what our own rosetta/scripts/
        run_protac_model.py driver adds to sys.path), or getProtacModelScript('main.py')
        for a specific file. Same subfolder reasoning as getFCCScript above: PROTAC_MODEL_HOME
        is the parent InstallHelper cloned into, the repo itself is one level down. """
        home = os.path.join(cls._requireToolHome(PROTAC_MODEL_DIC), PROTAC_MODEL_DIC['name'])
        return os.path.join(home, scriptName) if scriptName else home

    @classmethod
    def getProtacModelPython(cls):
        """ Path to the Python 2.7 (+ RDKit) interpreter dedicated to PROTAC-Model's own
        code - only used to check the env is installed (see _validate in the protocol),
        not to launch anything: runCondaScript() below activates this same conda env
        instead, so the env's own 'python' lands on PATH (a bare absolute-interpreter
        launch wouldn't, and FCC's clustering scripts need that). """
        home = cls._requireToolHome(PROTAC_MODEL_PYTHON_DIC)
        path = os.path.join(home, 'bin', 'python2')
        if not os.path.exists(path):
            raise FileNotFoundError(f'python2 not found under PROTAC_MODEL_PYTHON_HOME/bin ({home}).')
        return path

    @classmethod
    def getProtacModelEnviron(cls):
        """ PROTAC-Model's own utils/*.py modules read these bare names (no _HOME suffix)
        from os.environ at import time - this translates our *_HOME variables into that
        convention in one place, for use as runCondaScript()'s extraEnvDict. ROSETTA is
        required unconditionally (not just when refining): utils/rosetta.py is imported
        by our driver script at module load for every phase, so ROSETTA must resolve even
        for a frodock-only run. """
        return {
            'FRODOCK': cls._requireToolHome(FRODOCK_DIC),
            'ADFRSUITE': cls._requireToolHome(ADFRSUITE_DIC),
            'VINA': cls._requireToolHome(VINA_DIC),
            'VOROMQA': cls._requireToolHome(VOROMQA_DIC),
            # FCC_HOME is the parent InstallHelper cloned into; the actual FCC root
            # PROTAC-Model expects is one level down - same reasoning as getFCCScript().
            'FCC': os.path.join(cls._requireToolHome(FCC_DIC), FCC_DIC['name']),
            # XXX unverified: whether ROSETTA_HOME's auto-detected rosetta_bin_linux*
            # layout matches the $ROSETTA/main/source/bin/... path utils/rosetta.py
            # expects - no real Rosetta install available locally to check.
            'ROSETTA': cls._requireToolHome(ROSETTA_DIC),
            'PROTAC_MODEL_HOME': cls.getProtacModelScript(),
        }

    @classmethod
    def getPluginScript(cls, scriptName):
        """ Path to a script bundled with this plugin itself (rosetta/scripts/<scriptName>),
        e.g. run_protac_model.py - mirrors pwchem's Plugin.getScriptsDir(). """
        return os.path.join(os.path.dirname(__file__), 'scripts', scriptName)

    @classmethod
    def runCondaScript(cls, scriptPath, args, condaDic, extraEnvDict=None, cwd=None):
        """ Launch a Python script with a conda env activated first, instead of calling
        that env's interpreter by absolute path - mirrors pwchem's Plugin.runScript(). The
        difference matters here: run_protac_model.py (--phase filter) ends up shelling out
        to FCC's clustering scripts as a bare 'python <script>.py' command (see
        PROTAC-Model's own preprocess.py) - that only resolves to the right Python 2.7 if
        this env's bin/ is actually on PATH, which activating it does and launching by
        absolute interpreter path alone would not. """
        program = f'{cls.getEnvActivationCommand(condaDic)} && python {scriptPath}'
        cls.runProgram(program, args, extraEnvDict=extraEnvDict, cwd=cwd)

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
        rosettaHome = cls.getVar(ROSETTA_DIC['home'])
        if not os.path.exists(os.path.expanduser(rosettaHome)):
            missingPaths.append(f"Path of Rosetta does not exist ({ROSETTA_DIC['home']}): {rosettaHome}")
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
