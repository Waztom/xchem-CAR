from rest_framework import routers
from reactions.api import (
    ProjectViewSet, 
    TargetViewSet, 
    MethodViewSet, 
    ReactionViewSet, 
    ReactantViewSet,
    AddActionViewSet,
    MakeSolutionActionViewSet,
    StirActionViewSet,
    WashActionViewSet,
    DrySolutionActionViewSet,
    ConcentrateActionViewSet,
    AnalyseActionViewSet)

router = routers.DefaultRouter()
router.register('api/projects', ProjectViewSet, 'projects')
router.register('api/targets', TargetViewSet, 'targets')
router.register('api/methods', MethodViewSet, 'methods')
router.register('api/reactions', ReactionViewSet, 'reactions')
router.register('api/reactants', ReactantViewSet, 'reactants')
router.register('api/addactions', AddActionViewSet, 'addactions')
router.register('api/makesolutionactions', MakeSolutionActionViewSet, 'makesolutionactions')
router.register('api/stiractions', StirActionViewSet, 'stiractions')
router.register('api/washactions', WashActionViewSet, 'washactions')
router.register('api/drysolution', DrySolutionActionViewSet, 'drysolutionactions')
router.register('api/concentrateactions', ConcentrateActionViewSet, 'concentrateactions')
router.register('api/analyseactions', AnalyseActionViewSet, 'analyseactions')

urlpatterns = router.urls
