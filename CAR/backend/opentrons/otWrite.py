# this is vunrable to python injection by the lack of checking of metadata inputs
# this opens and closes files frequently, could be improved by creating string to hold the file data before writing to file once

class otScript():
    def __init__(self, filepath, protocolName=None, author=None, description=None, apiLevel='2.9'):
        self.filepath = filepath
        self.protocolName = protocolName
        self.author = author
        self.description = description
        self.apiLevel = apiLevel

        self.humanreadable = []

        #self.setupScript()

    def setupScript(self):
        """This is vunrable to injection atacks """
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
        script.write("\n# metadata\n metadata = {")
        script.write("\n\t'protocolName': '"+str(self.protocolName)+"',")
        script.write("\n\t'author': '"+str(self.author)+"',")
        script.write("\n\t'description': '"+str(self.description)+"',")
        script.write("\n\t'apiLevel': '"+str(self.apiLevel)+"'}\n")

        script.write("\ndef run(protocol: protocol_api.ProtocolContext):\n")

        script.close()
    
    def setupLabwear(self, platelist, tipOutput):
        script = open(self.filepath, "a")
        script.write("\n\t# labware")
        for plate in platelist:
            if plate.plateName == "":
                uniquename = (str("plate_"+str(plate.plateIndex))).replace(" ", "")
            else:
                uniquename = str(str(plate.plateName)+"_"+str(plate.plateIndex)).replace(" ","")
            script.write("\n\t"+uniquename+" = protocol.load_labwear('"+str(plate.plateTypeName)+"', '"+str(plate.plateIndex)+"')")
        
        script.close()

    def setupPipettes(self):
        script = open(self.filepath, "a")
        script.write("\n\t# pipettes\n")
        script.write("\tleft_pipette = protocol.load_instrument('"+str()+"', '"+str()+"', '"+str()+"')\n")

        script.close()

    def writeCommand(self, comandString):
        script = open(self.filepath, "a")
        if type(comandString) == str:
            script.write("\t"+str(commandString)+"\n")
        elif type(comandString) == list:
            for command in comandString:
                script.write("\t"+str(command)+"\n")

        script.close()
    
    def movefluids(self, pipetteName, fromAdress, toAdress, volume, dispenseVolume=None, writetoscript=True):
        if dispenseVolume == None:
            dispenseVolume = volume
            humanread = "move - "+str(volume)+"ul from "+str(fromAdress)+" to "+str(toAdress)
        else:
            humanread = "move - remove "+str(volume)+"ul from "+str(fromAdress)+" and add "+str(dispenseVolume)+"ul to "+str(toAdress)

        self.humanreadable.append(['OT',humanread])

        moveCommands = ["\n\t# "+str(humanread),
            str(pipetteName)+".pick_up_tip()",
            str(pipetteName)+".aspirate("+str(volume)+", "+str(fromAdress)+")",
            str(pipetteName)+".aspirate("+str(dispenseVolume)+", "+str(toAdress)+")",
            str(pipetteName)+".drop_tip()"]
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
        script.write("\n\t# pause for "+str(actionDetails)+"\n")
        script.write("\tprotocol.pause('"+str(actionDetails)+"')\n")


        humanread = actionDetails

        self.humanreadable.append(['User',humanread])
        script.close()
        


    

    def quickSetup(self):
        self.setupScript()
        self.setupLabwear()
        self.setupPipettes()

    def example(self):
        quickSetup()
        self.movefluids("left", "plate['A1']", "plate['B2']",100)



    def convertReactionToActions(self, reaction):
        print("\n\n\t\tReaction: "+str(reaction))
        for step in reaction:
            if step == 'add':
                self.movefluids('pipetteName', 'fromAdress', 'toAdress', '10')
            else:
                self.unsuportedAction(step)


    def convertReactionsToScript(self,list):
        for reaction in list:
            script = open(self.filepath, "a")
            script.write("\n\n\t# Reaction no: "+str(reaction[0])+"\n")
            script.close()
            print("\n\n\tReaction no: "+str(reaction[0]))
            self.convertReactionToActions(reaction[1])
        