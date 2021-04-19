# this relies on pandas, can be improved by moving full to sql

from django.http import JsonResponse
import backend.models
import pandas as pd
import numpy as np
import backend.opentrons.otWrite


def getactions():
    """this, using pandas, is not a sustanable aproach and will proberbly not scale"""
    IBMAddAction = pd.DataFrame(list(backend.models.IBMAddAction.objects.all().values()))
    IBMCollectLayerAction = pd.DataFrame(list(backend.models.IBMCollectLayerAction.objects.all().values()))
    AllActions = IBMAddAction.append(IBMCollectLayerAction)

    IBMConcentrateAction = pd.DataFrame(list(backend.models.IBMConcentrateAction.objects.all().values()))
    AllActions = AllActions.append(IBMConcentrateAction)
    IBMDegasAction = pd.DataFrame(list(backend.models.IBMDegasAction.objects.all().values()))
    AllActions = AllActions.append(IBMDegasAction)
    IBMDrySolidAction = pd.DataFrame(list(backend.models.IBMDrySolidAction.objects.all().values()))
    AllActions = AllActions.append(IBMDrySolidAction)
    IBMExtractAction = pd.DataFrame(list(backend.models.IBMExtractAction.objects.all().values()))
    AllActions = AllActions.append(IBMExtractAction)
    IBMFilterAction = pd.DataFrame(list(backend.models.IBMFilterAction.objects.all().values()))
    AllActions = AllActions.append(IBMFilterAction)
    IBMMakeSolutionAction = pd.DataFrame(list(backend.models.IBMMakeSolutionAction.objects.all().values()))
    AllActions = AllActions.append(IBMMakeSolutionAction)
    IBMPartitionAction = pd.DataFrame(list(backend.models.IBMPartitionAction.objects.all().values()))
    AllActions = AllActions.append(IBMPartitionAction)

    IBMpHAction = pd.DataFrame(list(backend.models.IBMpHAction.objects.all().values()))
    AllActions = AllActions.append(IBMpHAction)
    IBMPhaseSeparationAction = pd.DataFrame(list(backend.models.IBMPhaseSeparationAction.objects.all().values()))
    AllActions = AllActions.append(IBMPhaseSeparationAction)
    IBMQuenchAction = pd.DataFrame(list(backend.models.IBMQuenchAction.objects.all().values()))
    AllActions = AllActions.append(IBMQuenchAction)
    IBMRefluxAction = pd.DataFrame(list(backend.models.IBMRefluxAction.objects.all().values()))
    AllActions = AllActions.append(IBMRefluxAction)
    IBMSetTemperatureAction = pd.DataFrame(list(backend.models.IBMSetTemperatureAction.objects.all().values()))
    AllActions = AllActions.append(IBMSetTemperatureAction)
    IBMStirAction = pd.DataFrame(list(backend.models.IBMStirAction.objects.all().values()))
    AllActions = AllActions.append(IBMStirAction)
    IBMStoreAction = pd.DataFrame(list(backend.models.IBMStoreAction.objects.all().values()))
    AllActions = AllActions.append(IBMStoreAction)
    IBMWaitAction = pd.DataFrame(list(backend.models.IBMWaitAction.objects.all().values()))
    AllActions = AllActions.append(IBMWaitAction)
    IBMWashAction = pd.DataFrame(list(backend.models.IBMWashAction.objects.all().values()))
    AllActions = AllActions.append(IBMWashAction)

    #sort
    AllActions = AllActions.sort_values(['reaction_id_id', 'actionno'])

    print(AllActions)
    return AllActions

def getReactionActions(AllActions, ReactionID):
    if type(ReactionID) == int:
        ReactionID = [ReactionID]
    elif type(ReactionID) == np.int64:
        ReactionID = [ReactionID.item()]

        
    print(ReactionID)
    print(type(ReactionID))
    ReactionMask = AllActions.reaction_id_id.isin(ReactionID)
    ReactionActions = AllActions[ReactionMask]
    return(ReactionActions)


def basiccomprehension(AllActions):
    reactions = []
    for reaction in AllActions.reaction_id_id.unique():
        reactionactions = getReactionActions(AllActions, reaction)
        actionlist = []
        for index, row in reactionactions.iterrows():
            actionlist.append(row['actiontype']) 
        print(str(reaction)+" "+str(actionlist))
        reactions.append([reaction, actionlist])
    print(reactions)
    return reactions
    # return something here

