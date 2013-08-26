"""
Application-class that implements pyFoamInitMixingPlaneInterface.py

Initialize various mixingPlane interface attributes in the
constant/polymesh/boundary file, and in the time directories.

Backups of the modified files are created

Generate companion scripts for initializing the mixingPlane zone faceSets.

Modify the decomposeParDict file for the new mixingPlane zones names.

Author:
  Martin Beaudoin, Hydro-Quebec, 2013.  All rights reserved

"""

import sys, fnmatch, re
from os import path, listdir, chmod
from stat import *

from PyFoamApplication import PyFoamApplication
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.TimeDirectory import TimeDirectory
from PyFoam.Basics.BasicFile import BasicFile


class InitMixingPlaneInterface(PyFoamApplication):
    def __init__(self,args=None):
        description="""
Init MixingPlane boundary condition parameters
"""
        PyFoamApplication.__init__(self,
                                   args=args,
                                   description=description,
                                   usage="%prog <caseDirectory> mixingPlane_MasterPatchName mixingPlane_ShadowPatchName",
                                   interspersed=True,
                                   changeVersion=False,
                                   nr=3)

    def addOptions(self):
        self.parser.add_option("--coordinateSystemName",
                               action="store",
                               dest="coordinateSystemName",
                               default="mixingCS",
                               help='coordinateSystemName (mixingCS)')
        self.parser.add_option("--coordinateSystemType",
                               action="store",
                               dest="coordinateSystemType",
                               default=None,
                               help='coordinateSystemType (cyindrical/spherical)')
        self.parser.add_option("--coordinateSystemOrigin",
                               action="store",
                               dest="coordinateSystemOrigin",
                               default=None,
                               nargs=3,
                               type=float,
                               help='origin for coordinate system of mixingPlane. Specify a triplet of values separated by spaces')
        self.parser.add_option("--coordinateSystemE1",
                               action="store",
                               dest="coordinateSystemE1",
                               default=None,
                               nargs=3,
                               type=float,
                               help='origin for coordinate system of mixingPlane. Specify a triplet of values separated by spaces')
        self.parser.add_option("--coordinateSystemE3",
                               action="store",
                               dest="coordinateSystemE3",
                               default=None,
                               nargs=3,
                               type=float,
                               help='origin for coordinate system of mixingPlane. Specify a triplet of values separated by spaces')
        self.parser.add_option("--ribbonPatchSweepAxis",
                               action="store",
                               dest="ribbonPatchSweepAxis",
                               default=None,
                               help='ribbonPatch sweepAxis (X|Y|Z|R|Theta')
        self.parser.add_option("--ribbonPatchStackAxis",
                               action="store",
                               dest="ribbonPatchStackAxis",
                               default=None,
                               help='ribbonPatch stackAxis (X|Y|Z|R|Theta')
        self.parser.add_option("--ribbonPatchDiscretisation",
                               action="store",
                               dest="ribbonPatchDiscretisation",
                               default=None,
                               help='ribbonPatch discretisation (masterPatch|slavePatch|bothPatches|uniform|userDefined)')

        self.parser.add_option("--timeDirs",
                               action="store",
                               dest="timeDirs",
                               default=None,
                               help='time directories for the mixingPlane boundaryfields. Accept expressions like "[0-9]*", "0", etc.')

        self.parser.add_option("--genFaceSetForMixingPlaneZonesScriptName",
                               action="store",
                               dest="genFaceSetForMixingPlaneZonesScriptName",
                               default="genFaceSetForMixingPlaneZones.setSet",
                               help='setSet batch file for generating faceSets for mixingPlane zones. Default: genFaceSetForMixingPlaneZones.setSet')

        self.parser.add_option("--initMixingPlaneZonesScriptName",
                               action="store",
                               dest="initMixingPlaneZonesScriptName",
                               default="initMixingPlaneZones.sh",
                               help='script name for initializing the mixingPlane zone faceSets. Default: initMixingPlaneZones.sh')

        self.parser.add_option("--test",
                               action="store_true",
                               default=False,
                               dest="test",
                               help="Only print the new boundary file")

    def createMixingPlanePatch(self, patch, patchName):
        description="""\
Create a default definition for a mixingPlane patch, and replace
the current definition
        """
        print "Replacing definition of patch: ", patchName, ":", patch
        newPatch={
            'type'        : "mixingPlane",
            'nFaces'      : patch["nFaces"],
            'startFace'   : patch["startFace"],
            'shadowPatch' : 'unknown',
            'zone'        : patchName+'Zone',
            'coordinateSystem' : {
                'name'   : 'mixingCS',
                'type'   : 'cylindrical',
                'origin' : '(0 0 0)',
                'e1'     : '(1 0 0)',
                'e3'     : '(0 0 1)'
                },
            'ribbonPatch' : {
                'sweepAxis'      : 'Theta',
                'stackAxis'      : 'Z',
                'discretisation' : 'bothPatches',
                }
            }
        return newPatch

    def modifyMixinPlanePatchDefinition(self, patch, patchName, shadowName):
        description="""\
Modify the definition of a mixingPlane patch
        """
        print "    Modifying mixingPlane boundary definition in constant/polyMesh/boundary for patch", patchName

        patch["shadowPatch"]=shadowName

        patch["zone"]=patchName+'Zone'

        if patch.has_key("coordinateSystem")==False:
            patch["coordinateSystem"]={}

        if self.parser.getOptions().coordinateSystemName!=None:
            patch["coordinateSystem"]["name"]=self.parser.getOptions().coordinateSystemName

        if self.parser.getOptions().coordinateSystemType!=None:
            patch["coordinateSystem"]["type"]=self.parser.getOptions().coordinateSystemType

        if self.parser.getOptions().coordinateSystemOrigin!=None:
            patch["coordinateSystem"]["origin"]="(%f %f %f)" % self.parser.getOptions().coordinateSystemOrigin

        if self.parser.getOptions().coordinateSystemE1!=None:
            patch["coordinateSystem"]["e1"]="(%f %f %f)" % self.parser.getOptions().coordinateSystemE1

        if self.parser.getOptions().coordinateSystemE3!=None:
            patch["coordinateSystem"]["e3"]="(%f %f %f)" % self.parser.getOptions().coordinateSystemE3

        if patch.has_key("ribbonPatch")==False:
            patch["ribbonPatch"]={}

        if self.parser.getOptions().ribbonPatchSweepAxis!=None:
            patch["ribbonPatch"]["sweepAxis"]=self.parser.getOptions().ribbonPatchSweepAxis

        if self.parser.getOptions().ribbonPatchStackAxis!=None:
            patch["ribbonPatch"]["stackAxis"]=self.parser.getOptions().ribbonPatchStackAxis

        if self.parser.getOptions().ribbonPatchDiscretisation!=None:
            patch["ribbonPatch"]["discretisation"]=self.parser.getOptions().ribbonPatchDiscretisation



    def modifyMixinPlanePatchDefinitionInTimeDirs(self, caseDir, patchName, timeDirs):
        description="""\
Modify the definition of a mixingPlane patch in the time directories
        """
        regex = fnmatch.translate(timeDirs)

        reobj = re.compile(regex)

        for timeDir in listdir(caseDir):
            if reobj.match(timeDir):
                print "    Modifying mixingPlane boundaryFields in timeDir", timeDir, "for patch", patchName
                td=TimeDirectory(caseDir, timeDir, yieldParsedFiles=True)

                for f in td:
                    print "        Modifying field", f.name
                    f["boundaryField"][patchName]["type"]='mixingPlane'
                    #
                    # With the current version, the mixingType is configured
                    # using system/fvSchemes
                    #f["boundaryField"][patchName]["fluxAveraging"]='false'
                    f.writeFile()

    def generateCompanionFiles(self, caseDir, boundary):
        description="""\
Generate a setSet batch file based on the zone info specified in the mixingPlane interfaces definition.
Generate a bash file for invoking setSet and setsToZones
Update mixingPlane zone information in decomposeParDict
        """
        # Default file: genFaceSetForMixingPlaneZones.setSet
        bfGenFaceSets = BasicFile(path.join(caseDir, self.parser.getOptions().genFaceSetForMixingPlaneZonesScriptName))

        print "    Updating file ", bfGenFaceSets.name, " for generating mixingPlane zones faceSet using the setSet command"

        bnd=boundary.content

        if type(bnd)!=list:
            self.error("Problem with boundary file (not a list)")

        # Memorize list of mixingPlane zones for later processing
        listOfMixingPlaneZones = []

        for index in range(0, len(bnd), 2):
            patchName = bnd[index]
            indexDefPatch=index+1
            if bnd[indexDefPatch]["type"]=='mixingPlane':
                bfGenFaceSets.writeLine([ "faceSet " + bnd[indexDefPatch]["zone"] + " new patchToFace "+ patchName ])
                listOfMixingPlaneZones.append(bnd[indexDefPatch]["zone"])

        bfGenFaceSets.writeLine([ "quit" ])
        bfGenFaceSets.close()

        # Default file: initMixingPlaneZones.sh
        bfInitMixingPlaneZones = BasicFile(path.join(caseDir, self.parser.getOptions().initMixingPlaneZonesScriptName))

        print "    Updating file ", bfInitMixingPlaneZones.name, " for inititalizing mixingPlane zones"

        bfInitMixingPlaneZones.writeLine([ "#!/bin/bash" ])
        bfInitMixingPlaneZones.writeLine([ "setSet -batch " + self.parser.getOptions().genFaceSetForMixingPlaneZonesScriptName ])
        bfInitMixingPlaneZones.writeLine([ "setsToZones -noFlipMap" ])
        bfInitMixingPlaneZones.close()

        # Set execution permissions for this script (755)
        chmod(bfInitMixingPlaneZones.name, S_IRWXU|S_IRGRP|S_IXGRP|S_IXOTH|S_IROTH)

        # DecomposeParDict
        decomposeParDictPath=path.join(caseDir,"system","decomposeParDict")
        if path.exists(decomposeParDictPath):
            print "    Updating file ", decomposeParDictPath, " for mixingPlane zones"
            decomposeParDict=ParsedParameterFile(decomposeParDictPath,debug=False,backup=True)
            dcp=decomposeParDict.content
            # Adding to existing list of zones, making sure we have a list of unique values
            uniqZoneValues=set(dcp["globalFaceZones"] + listOfMixingPlaneZones)
            dcp["globalFaceZones"]="(\n    " + '\n    '.join(uniqZoneValues) + "\n)"
            decomposeParDict.writeFile()

    def run(self):
        caseDir=self.parser.getArgs()[0]
        masterbName=self.parser.getArgs()[1]
        shadowbName=self.parser.getArgs()[2]

        boundary=ParsedParameterFile(path.join(".",caseDir,"constant","polyMesh","boundary"),debug=False,boundaryDict=True,backup=True)

        bnd=boundary.content

        if type(bnd)!=list:
            print "Problem with boundary file (not a list)"
            sys.exit(-1)

        masterFound=False
        shadowFound=False
        updateTimeDirs=False

        timeDirs="0"
        if self.parser.getOptions().timeDirs!=None:
            timeDirs=self.parser.getOptions().timeDirs
            updateTimeDirs=True

        print "UpdateTimeDirs: ", updateTimeDirs

        for index in range(len(bnd)):

            indexDefPatch=index+1

            if bnd[index]==masterbName:
                masterFound=True
                if bnd[indexDefPatch]["type"]!="mixingPlane":
                    bnd[indexDefPatch] = self.createMixingPlanePatch(bnd[indexDefPatch], masterbName)

                self.modifyMixinPlanePatchDefinition(bnd[indexDefPatch], masterbName, shadowbName)

                if updateTimeDirs:
                    self.modifyMixinPlanePatchDefinitionInTimeDirs(caseDir, masterbName, timeDirs)

            elif bnd[index]==shadowbName:
                shadowFound=True
                if bnd[indexDefPatch]["type"]!="mixingPlane":
                    bnd[indexDefPatch] = self.createMixingPlanePatch(bnd[indexDefPatch], shadowbName)

                self.modifyMixinPlanePatchDefinition(bnd[indexDefPatch], shadowbName, masterbName)

                if updateTimeDirs:
                    self.modifyMixinPlanePatchDefinitionInTimeDirs(caseDir, shadowbName, timeDirs)

            if masterFound and shadowFound:
                break;

        if not masterFound:
            self.error("Boundary patch",masterbName,"not found in",bnd[::2])

        if not shadowFound:
            self.error("Boundary patch",shadowbName,"not found in",bnd[::2])

        if self.parser.getOptions().test:
            print boundary
        else:
            boundary.writeFile()

        # Write companion files
        self.generateCompanionFiles(caseDir, boundary)


