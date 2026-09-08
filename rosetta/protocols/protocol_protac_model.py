# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Gutiérrez (Internship student, github user: danigual)
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

"""
Wrapper around PROTAC-Model (Gao et al., Zhejiang University,
https://github.com/gaoqiweng/PROTAC-Model), an integrative modeling pipeline for
PROTAC-mediated protein-protein ternary complexes:

  1. FRODOCK: global rigid-body protein-protein docking, guided by a site point on the
     receptor interface.
  2. Filtering of the FRODOCK poses by compatibility with the PROTAC SMILES (and,
     optionally, up to two E3 ligand SDF files, for PROTACs with more than one possible
     anchoring orientation).
  3. Optional RosettaDock refinement of the filtered poses (slower, more accurate).

TODO: this is a skeleton only. FRODOCK/ADFRsuite/Vina/Voromqa/FCC are not wired in yet -
see CLAUDE.md (local, untracked) for the full pipeline breakdown and open questions.
"""

import os

from pyworkflow.constants import BETA
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem.objects import SmallMolecule, SetOfAtomStructsChem
from pwchem.utils import cleanPDB

from rosetta import Plugin
from rosetta.constants import ROSETTA_SCRIPTS


class RosettaProtPROTACModel(EMProtocol):
    """
    Models a PROTAC-mediated protein-protein ternary complex (receptor + target,
    bridged by a PROTAC) using the PROTAC-Model pipeline (Gao et al.): FRODOCK global
    docking guided by a receptor-interface site point, filtering of the resulting poses
    by compatibility with the PROTAC linker geometry, and optional RosettaDock
    refinement of the surviving poses.
    """
    _label = 'PROTAC ternary complex modeling (PROTAC-Model)'
    _devStatus = BETA

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Proteins')
        group.addParam('inputReceptor', params.PointerParam, pointerClass='AtomStruct',
                       label='Receptor structure', allowsNull=False,
                       help='Larger of the two proteins (typically the E3 ligase). Passed '
                            'as the FRODOCK receptor. Per PROTAC-Model\'s own requirements, '
                            'this structure should already include its bound small-molecule '
                            'warhead (as HETATM records) and exclude unrelated heteroatoms '
                            '(crystallization waters, ions, buffer molecules, etc.).')
        group.addParam('inputTarget', params.PointerParam, pointerClass='AtomStruct',
                       label='Target structure', allowsNull=False,
                       help='Smaller of the two proteins (typically the protein of '
                            'interest, POI). Passed as the FRODOCK docking target. Same '
                            'requirement as "Receptor structure": keep its bound warhead, '
                            'exclude unrelated heteroatoms.')
        group.addParam('siteCoords', params.StringParam, allowsNull=False,
                       label='Receptor interface site (X,Y,Z)',
                       help='Coordinates of a point on the receptor surface, at the '
                            'interface where the ternary complex is expected to form. '
                            'FRODOCK uses this point to restrict the global docking search.')

        group = form.addGroup('PROTAC')
        group.addParam('protacSmiles', params.StringParam, allowsNull=False,
                       label='PROTAC SMILES',
                       help='SMILES of the full PROTAC molecule (E3 ligand - linker - '
                            'POI warhead), used to filter FRODOCK poses by compatibility '
                            'with the linker geometry.')
        group.addParam('receptorLigandSmiles', params.StringParam, allowsNull=True,
                       label='Receptor-bound ligand SMILES (optional)',
                       help='SMILES of the small-molecule warhead already bound in the '
                            'receptor structure (equivalent to PROTAC-Model\'s -irsmi). '
                            'RDKit can sometimes fail to assign the correct bonds for this '
                            'ligand when reading it straight from the PDB HETATM records; '
                            'providing its SMILES avoids that failure mode.')
        group.addParam('targetLigandSmiles', params.StringParam, allowsNull=True,
                       label='Target-bound ligand SMILES (optional)',
                       help='SMILES of the small-molecule warhead already bound in the '
                            'target structure (equivalent to PROTAC-Model\'s -itsmi). Same '
                            'purpose as "Receptor-bound ligand SMILES", for the other '
                            'protein.')
        group.addParam('e3Ligand1', params.PointerParam, pointerClass='SmallMolecule',
                       allowsNull=True, label='E3 ligand conformer 1 (optional)',
                       help='SDF of one possible bound conformation/orientation of the '
                            'E3-binding warhead. Optional: only needed to disambiguate '
                            'PROTACs whose E3 ligand can anchor in more than one way '
                            '(e.g. thalidomide-based degraders).')
        group.addParam('e3Ligand2', params.PointerParam, pointerClass='SmallMolecule',
                       allowsNull=True, label='E3 ligand conformer 2 (optional)',
                       condition='e3Ligand1 is not None',
                       help='SDF of a second possible bound conformation/orientation of '
                            'the E3-binding warhead, as an alternative to conformer 1.')

        group = form.addGroup('Docking parameters')
        group.addParam('doRefine', params.BooleanParam, default=False,
                       label='Refine with RosettaDock',
                       help='Refine the filtered FRODOCK poses with RosettaDock. Improves '
                            'accuracy but is considerably slower than FRODOCK + filtering '
                            'alone.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        convertId = self._insertFunctionStep(self.convertInputStep, prerequisites=[])
        dockId = self._insertFunctionStep(self.frodockStep, prerequisites=[convertId])
        filterId = self._insertFunctionStep(self.filterPosesStep, prerequisites=[dockId])

        lastId = filterId
        if self.doRefine.get():
            lastId = self._insertFunctionStep(self.refineStep, prerequisites=[filterId])

        self._insertFunctionStep(self.createOutputStep, prerequisites=[lastId])

    def convertInputStep(self):
        """ Converts the receptor/target AtomStructs into the clean PDBs that FRODOCK
        expects, and stages the PROTAC SMILES into extraPath.
        Only water is stripped here: per PROTAC-Model's own input requirements, the
        receptor/target PDBs must keep their bound small-molecule warhead (a HETATM
        record), so heteroatoms as a whole cannot be blanket-removed the way
        convertInputStep did before.
        TODO: once the warhead's residue name is known/identifiable, also strip other
        unrelated heteroatoms (ions, buffer molecules, etc.) via cleanPDB's het2rem,
        instead of keeping every non-water heteroatom. """
        self.receptorFile = self._getExtraPath('receptor.pdb')
        cleanPDB(self.inputReceptor.get().getFileName(), self.receptorFile,
                waters=True, hetatm=False)

        self.targetFile = self._getExtraPath('target.pdb')
        cleanPDB(self.inputTarget.get().getFileName(), self.targetFile,
                waters=True, hetatm=False)

        self.smilesFile = self._getExtraPath('protac.smi')
        with open(self.smilesFile, 'w') as f:
            f.write(self.protacSmiles.get().strip() + '\n')

    def frodockStep(self):
        # STAGE 1: receptor van der Waals potential map.
        # getFrodockProgram resolves FRODOCK_HOME/bin/frodockgrid (or its _gcc fallback).
        frodockgrid = Plugin.getFrodockProgram('frodockgrid')
        # self.receptorFile was written in convertInputStep (the cleaned receptor PDB).
        # Output map goes to extraPath, e.g. <protocol>/extra/receptor_W.ccp4.
        args = '%s -o %s' % (self.receptorFile, self._getExtraPath('receptor_W.ccp4'))
        # runProgram just sets up the environment and launches frodockgrid with those args.
        # cwd=extraPath: matches the rest of the repo's convention (other protocols always
        # pass cwd), and matters here because we don't yet know whether FRODOCK writes any
        # of its by-product files (e.g. the _ASA.pdb from stage 5's TODO) relative to CWD
        # instead of next to the -o path - pinning CWD to extraPath keeps that debuggable.
        Plugin.runProgram(frodockgrid, args, cwd=self._getExtraPath())

        # STAGE 2: receptor electrostatic potential map (-m 1 -t E, per FRODOCK's own docs).
        # Same binary as stage 1, just different flags/output file.
        args = '%s -o %s -m 1 -t E' % (self.receptorFile, self._getExtraPath('receptor_E.ccp4'))
        Plugin.runProgram(frodockgrid, args, cwd=self._getExtraPath())

        # STAGE 3: receptor desolvation potential map.
        args = '%s -o %s -m 3' % (self.receptorFile, self._getExtraPath('receptor_DS.ccp4'))
        Plugin.runProgram(frodockgrid, args, cwd=self._getExtraPath())

        # STAGE 4: target desolvation potential map. Target, not receptor this time - the
        # docking stage (5) needs a desolvation map for both molecules.
        args = '%s -o %s -m 3' % (self.targetFile, self._getExtraPath('target_DS.ccp4'))
        Plugin.runProgram(frodockgrid, args, cwd=self._getExtraPath())

        # STAGE 5: the actual global rigid-body docking, restricted around siteCoords.
        frodock = Plugin.getFrodockProgram('frodock')
        receptorASA = self._getExtraPath('receptor_ASA.pdb')
        targetASA = self._getExtraPath('target_ASA.pdb')
        # TODO: assumed side-effect of frodockgrid (stages 1-4) writing these _ASA.pdb
        # files alongside the .ccp4 maps - unverified, no FRODOCK binary available
        # locally to test.
        args = '%s %s -w %s -e %s --th 10 -d %s,%s -t E -o %s --around %s' % (
            receptorASA, targetASA,
            self._getExtraPath('receptor_W.ccp4'), self._getExtraPath('receptor_E.ccp4'),
            self._getExtraPath('receptor_DS.ccp4'), self._getExtraPath('target_DS.ccp4'),
            self._getExtraPath('dock.dat'),
            # Reformatted from the already-parsed tuple (not the raw form string) so stray
            # whitespace (e.g. "12.3, -4.5, 6.7", which passes _validate's comma-count
            # check) can't split --around into extra shell tokens.
            '%.4f,%.4f,%.4f' % self._getSiteCoords())
        Plugin.runProgram(frodock, args, cwd=self._getExtraPath())

        # STAGE 6: cluster the raw docking poses from dock.dat (by RMSD, -d 4 Angstrom
        # cutoff), keeping representative poses instead of thousands of near-duplicates.
        frodockcluster = Plugin.getFrodockProgram('frodockcluster')
        args = '%s %s -o %s -d 4 --nc 100000' % (
            self._getExtraPath('dock.dat'), self.targetFile,
            self._getExtraPath('clust_dock_4.dat'))
        Plugin.runProgram(frodockcluster, args, cwd=self._getExtraPath())

        # STAGE 7: extract per-pose coordinates/scores from the clustered results.
        frodockview = Plugin.getFrodockProgram('frodockview')
        # '>' instead of the original script's '>>': a Scipion step can be re-run
        # (retry/resume), and append would duplicate scores into frodock_score.txt
        # across runs.
        args = '%s -p %s > %s' % (
            self._getExtraPath('clust_dock_4.dat'), self.targetFile,
            self._getExtraPath('frodock_score.txt'))
        Plugin.runProgram(frodockview, args, cwd=self._getExtraPath())

    def filterPosesStep(self):
        self._prepareLigandsForFiltering()
        # TODO: A2 (bond order assignment via RDKit/SMILES or OpenBabel), A3 (protonation
        # with reduce), A4 (interface residue calculation), then the per-pose filtering
        # loop and the ranking/clustering block. See PROTAC_PROGRESS_LOG.txt.
        raise NotImplementedError

    def _prepareLigandsForFiltering(self):
        """ Block A of PROTAC-Model's filter_frodock(): one-time preparation shared by
        every FRODOCK pose, before the per-pose filtering loop. """
        # A1: pull the warhead (HETATM) atoms out of receptor.pdb/target.pdb into their
        # own PDB files, so later steps (bond order assignment, interface distance
        # calculations) can work on just the small molecule.
        self._extractLigandPDB(self.receptorFile, self._getExtraPath('rec_lig.pdb'))
        self._extractLigandPDB(self.targetFile, self._getExtraPath('target_lig.pdb'))

    @staticmethod
    def _extractLigandPDB(inputPdb, outputPdb):
        """ Port of PROTAC-Model's preprocess.py::preprocess_pdb_element (Python 2 -> 3,
        same logic). Extracts HETATM lines (the bound warhead, since waters/other
        heteroatoms were already stripped in convertInputStep) and rewrites each one with
        a corrected element symbol in PDB columns 77-78 (derived from the atom-name field
        in columns 13-14, stripped of any trailing digits) - RDKit/OpenBabel need that
        column to be right to perceive the molecule's chemistry correctly; raw PDB HETATM
        records often leave it blank or wrong. """
        with open(inputPdb) as f:
            pdbLines = f.read().splitlines()

        outLines = []
        for line in pdbLines:
            # Column slice [12:14] is the atom-name field (e.g. "C1", "N2", "H12"); skip
            # any HETATM line where that field is empty.
            if line[:6] == 'HETATM' and line[12:14].strip():
                element = line[12:14].translate(str.maketrans('', '', '0123456789'))
                if element[0] == 'H':
                    element = ' H'
                elif len(element) == 1:
                    element = ' %s' % element
                # Rebuild the line: keep columns 1-76 as-is, pad/truncate to 76, then
                # place the fixed element symbol in columns 77-78 (PDB spec).
                line = line[:76]
                line = line + ' ' * (76 - len(line)) + element
                outLines.append(line)

        with open(outputPdb, 'w') as f:
            f.write('\n'.join(outLines) + ('\n' if outLines else ''))

    def refineStep(self):
        # TODO: refine the filtered poses with RosettaDock via Plugin.runRosettaProgram.
        raise NotImplementedError

    def createOutputStep(self):
        # TODO: collect the final (filtered, optionally refined) ternary complex models
        # into a SetOfAtomStructsChem, one AtomStruct per model, with its docking/filter
        # score as an attribute.
        outputSet = SetOfAtomStructsChem.create(self._getPath())

        self._defineOutputs(outputTernaryModels=outputSet)
        self._defineSourceRelation(self.inputReceptor, outputSet)
        self._defineSourceRelation(self.inputTarget, outputSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []

        if self._getSiteCoords() is None:
            errors.append('"Receptor interface site (X,Y,Z)" must be 3 comma-separated '
                          'numbers, e.g. "12.3,-4.5,6.7". Got: "%s".' % self.siteCoords.get())

        if self.e3Ligand2.get() is not None and self.e3Ligand1.get() is None:
            errors.append('"E3 ligand conformer 2" was set without "E3 ligand conformer 1". '
                          'Set conformer 1 first, or clear conformer 2.')

        # TODO: once FRODOCK/RosettaDock are wired in, check their binaries are available.
        return errors

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['LeaverFay2011']

    # --------------------------- UTILS functions ------------------------------
    def _getSiteCoords(self):
        """ Parses siteCoords ("X,Y,Z") into a tuple of 3 floats, or None if malformed. """
        parts = self.siteCoords.get().strip().split(',')
        if len(parts) != 3:
            return None
        try:
            return tuple(float(p.strip()) for p in parts)
        except ValueError:
            return None


