# WTOSCR: could this become a module, currently this is just a placeholder, move to Opentrons folder?

# this is vunrable to python injection by the lack of checking of metadata inputs
# this opens and closes files frequently, could be improved by creating string to hold the file data before writing to file once

import os
from datetime import datetime

now = datetime.now()


class HumanReadable:
    def __init__(
        self, filepath, protocolName=None, author=None, description=None, fileFormat=["md"]
    ):
        self.dirsetup(path="../output/protocols")
        self.filepath = "../output/protocols/example.md"
        self.protocolName = protocolName
        self.author = author
        self.description = description
        self.fileFormat = fileFormat

    def setupDoc(self):
        self.dirsetup(path="../output/Opentrons")
        """This is vunrable to injection atacks """
        if self.protocolName == None:
            self.protocolName = input("Please name the Protocol Name: \t")
        if self.author == None:
            self.author = input("Please name the author Name: \t")
        if self.description == None:
            self.description = input("Please name the description: \t")

        if "txt" in self.fileFormat:
            file = open(f"{self.filepath}.txt", "w")
            file.write("from opentrons import protocol_api\n")
            file.write(
                "# "
                + str(self.protocolName)
                + ' for "'
                + str(self.author)
                + str('" produced by XChem Car (https://car.xchem.diamond.ac.uk)')
            )
            file.write("\n# metadata")
            file.write(
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
            file.write("\n\nfrom opentrons import protocol_api\n")
            file.write("\ndef run(protocol: protocol_api.ProtocolContext):\n")

            file.close()
        elif "md" in self.fileFormat:
            file = open(self.filepath, "w")
            file.write(f'# Protocol: "{self.protocolName}"\n')
            file.write(
                f"***Produced by XChem Car (https://car.xchem.diamond.ac.uk)**, for {self.author} on {now.strftime('%d/%m/%Y %H:%M:%S')}*\n"
            )
            if self.description != None:
                file.write(f"\n**Description:** \n*{self.description}*\n")

            file.close()

    def newBlock(self, blocknum, blockbool, blocksession):
        file = open(self.filepath, "a")
        if blockbool == True:
            file.write(f"## Block: '**{blocknum}**'  - *Robotic*\n")
            file.write("\n*To be executed on the **Robot***\n")
            file.write(str(blocksession.deck.PlateList))
        else:
            file.write(f"## Block: '**{blocknum}**'  - *Manual*\n")
            file.write("\n*To be compleated **manually***\n")
        file.close()

    def setupLabware(self, platelist, tipOutput):
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
            )

        script.close()

    def setupPipettes(self, pipettelist):
        script = open(self.filepath, "a")
        script.write("\n\n\t# pipettes\n")
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
            script.write("\t" + str(commandString) + "\n")
        elif type(comandString) == list:
            for command in comandString:
                script.write("\t" + str(command) + "\n")

        script.close()

    def dirsetup(self, path):
        if not os.path.exists(path):
            os.makedirs(path)

    def movefluids(
        self, pipetteName, fromAdress, toAdress, volume, dispenseVolume=None, writetoscript=True
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
    ):
        humanread = f"transfer - {volume}ul from {fromAdress} to {toAdress}"

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
        self.setupLabware()
        self.setupPipettes()

    def example(self):
        quickSetup()
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
