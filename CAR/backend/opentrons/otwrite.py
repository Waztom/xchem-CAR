# this is vunrable to python injection by the lack of checking of metadata inputs
# this opens and closes files frequently, could be improved by creating string to hold the file data before writing to file once

#   NOTE: in this file "humanread" referes to the comments above each line/set of lines of ot code, human readable is a list of all the comments in format [oporator (human/ot), comment]

import os


class otScript:
    def __init__(self, filepath, protocolName=None, author=None, description=None, apiLevel="2.9"):
        self.filepath = filepath
        self.protocolName = protocolName
        self.author = author
        self.description = description
        self.apiLevel = apiLevel

        self.humanreadable = []

    def setupScript(self):
        self.dirsetup()
        """This is vunrable to injection atacks """

        # asks user for metadata (should be removed at some point to avoid errors when intergrated with frontend)
        if self.protocolName == None:
            self.protocolName = input("Please name the Protocol Name: \t")
        if self.author == None:
            self.author = input("Please name the author Name: \t")
        if self.description == None:
            self.description = input("Please name the description: \t")
        if self.apiLevel == None:
            self.apiLevel = input("Please name the API Level: \t")

        script = open(self.filepath, "w")
        script.write("from opentrons import protocol_api\n")
        script.write(
            "# "
            + str(self.protocolName)
            + ' for "'
            + str(self.author)
            + str('" produced by XChem Car (https://car.xchem.diamond.ac.uk)')
        )
        script.write("\n# metadata")
        script.write(
            "\nmetadata = {'protocolName': '"
            + str(self.protocolName)
            + "', 'author': '"
            + str(self.author)
            + "','description': '"
            + str(self.description[0])
            + "','apiLevel': '"
            + str(self.apiLevel)
            + "'}\n"
        )
        script.write("\ndef run(protocol: protocol_api.ProtocolContext):\n")

        script.close()

    def setupLabware(self, platelist):
        script = open(self.filepath, "a")
        script.write("\n\t# labware")
        for plate in platelist:
            if plate.plateName == "":
                uniquename = (str("plate_" + str(plate.plateIndex))).replace(" ", "")
            else:
                if plate.platetype == "Plate":
                    uniquename = str(str(plate.plateName) + "_" + str(plate.plateIndex)).replace(
                        " ", ""
                    )

                    uniquename = plate.plateName  # debug line
                else:
                    uniquename = plate.plateName
            script.write(
                f"\n\t{uniquename} = protocol.load_labware('{plate.plateTypeName}', '{plate.plateIndex}')"
            )  # WTOSCR: may want to cast plate index to int type, but seems to work with str

        script.close()

    def setupPipettes(self, pipettelist):
        script = open(self.filepath, "a")
        script.write("\n\n\t# pipettes\n")
        mountnumber = 0
        for pipette in pipettelist:
            script.write(
                "\t"
                + str(pipette.name)
                + " = protocol.load_instrument('"
                + str(pipette.model)
                + "', '"
                + str(pipette.mount)
                + "', tip_racks="
                + str(pipette.tipRacks).replace("'", "")
                + ")\n"
            )
        script.close()

    def writeCommand(self, comandString):
        script = open(self.filepath, "a")
        if type(comandString) == str:
            script.write("\t" + str(comandString) + "\n")
        elif type(comandString) == list:
            for command in comandString:
                script.write("\t" + str(command) + "\n")

        script.close()

    def dirsetup(self, path="../output/Opentrons"):
        # WTOSCR: remove after django implementation
        if not os.path.exists(path):
            os.makedirs(path)

    def movefluids(
        self,
        pipetteName,
        fromAdress,
        toAdress,
        volume,
        dispenseVolume=None,
        writetoscript=True,
    ):
        if dispenseVolume == None:
            dispenseVolume = volume
            humanread = (
                "move - " + str(volume) + "ul from " + str(fromAdress) + " to " + str(toAdress)
            )
        else:
            humanread = (
                "move - remove "
                + str(volume)
                + "ul from "
                + str(fromAdress)
                + " and add "
                + str(dispenseVolume)
                + "ul to "
                + str(toAdress)
            )

        self.humanreadable.append(["OT", humanread])

        moveCommands = [
            "\n\t# " + str(humanread),
            str(pipetteName) + ".pick_up_tip()",
            str(pipetteName) + ".aspirate(" + str(volume) + ", " + str(fromAdress) + ")",
            str(pipetteName) + ".dispense(" + str(dispenseVolume) + ", " + str(toAdress) + ")",
            str(pipetteName) + ".drop_tip()",
        ]
        if writetoscript is True:
            self.writeCommand(moveCommands)
        return moveCommands

    def transferfluids(
        self,
        pipetteName,
        fromAdress,
        toAdress,
        volume,
        writetoscript=True,
        takedistance=2,
        dispensedistance=-5,
        fromalphanumeric=None,
        toalphanumeric=None,
    ):
        if fromalphanumeric != None and toalphanumeric != None:
            humanread = f"transfer - {volume:.1f}ul from {fromalphanumeric[0]}{fromalphanumeric[1]} to {toalphanumeric[0]}{toalphanumeric[1]}"
        else:
            humanread = f"transfer - {volume:.1f}ul from {fromAdress} to {toAdress}"

        self.humanreadable.append(["OT", humanread])

        moveCommands = [
            "\n\t# " + str(humanread),
            str(pipetteName)
            + f".transfer({volume}, {fromAdress}.bottom({takedistance}), {toAdress}.top({dispensedistance}), air_gap = 15)",
        ]
        if writetoscript is True:
            self.writeCommand(moveCommands)
        return moveCommands

    def readScript(self):
        script = open(self.filepath, "a")
        scriptContent = script.read()
        print(scriptContent)
        return scriptContent

    def unsuportedAction(self, actionDetails):
        script = open(self.filepath, "a")
        script.write("\n\t# pause for " + str(actionDetails) + "\n")
        script.write("\tprotocol.pause('" + str(actionDetails) + "')\n")

        humanread = actionDetails

        self.humanreadable.append(["User", humanread])
        script.close()

    def quickSetup(self):
        self.setupScript()
        self.setupLabwear()
        self.setupPipettes()

    def example(self):
        self.quickSetup()
        self.movefluids("left", "plate['A1']", "plate['B2']", 100)

    def convertReactionToActions(self, reaction):
        for step in reaction:
            if step == "add":
                self.movefluids("pipetteName", "fromAdress", "toAdress", "10")
            else:
                self.unsuportedAction(step)

    def convertReactionsToScript(self, list):
        for reaction in list:
            script = open(self.filepath, "a")
            script.write("\n\n\t# Reaction no: " + str(reaction[0]) + "\n")
            script.close()

            self.convertReactionToActions(reaction[1])
