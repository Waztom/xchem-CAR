from django.http import JsonResponse
from backend.models import Project   

from backend.models import Reaction
from backend.models import IBMAddAction
from backend.models import AnalyseAction


def getprojectnames():
    projects = Project.objects.all().values() 
    projects = list(projects)  
    print(projects)
    return JsonResponse(projects, safe=False)

def getreactions():
    reaction = Reaction.objects.all().values() 
    reaction = list(reaction)  
    print(reaction)
    return JsonResponse(reaction, safe=False)

def getiBMAddAction():
    iBMAddAction = IBMAddAction.objects.all().values() 
    iBMAddAction = list(iBMAddAction)  
    print(iBMAddAction)
    return JsonResponse(iBMAddAction, safe=False)
def getanalyseAction():
    analyseAction = AnalyseAction.objects.all().values() 
    analyseAction = list(analyseAction)  
    print(analyseAction)
    return JsonResponse(analyseAction, safe=False)
    
def getObjectJson(Object):
    from backend.models import Object
    objectout = Object.objects.all().values() 
    objectout = list(objectout)  
    print(objectout)
    return JsonResponse(objectout, safe=False)

getprojectnames()

def getallIBMactions():
    for modIBMaddactions, IBMcollect-layeractions,
     IBMconcentrateactions, IBMdegasactions, IBMdry-solidactions, 
     IBMdry-solutionactions, IBMextractactions, IBMfilteractions,
     IBMmake-solutionactions, IBMpartitionactions, IBMphactions,
     IBMphase-separationactions, IBMquenchactions, IBMrefluxactions,
     IBMset-temperatureactions, IBMstiractions, IBMstoreactions, 
     IBMwaitactions, IBMwashactions