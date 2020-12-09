from django.conf.urls import url, include
from rest_framework import routers
from reactions.api import (
    ProjectViewSet, 
    TargetViewSet, 
    MethodViewSet, 
    ReactionViewSet, 
    ProductViewSet,
    ReactantViewSet,
    AddActionViewSet,
    MakeSolutionActionViewSet,
    StirActionViewSet,
    WashActionViewSet,
    DrySolutionActionViewSet,
    ConcentrateActionViewSet,
    AnalyseActionViewSet)
from .views import UploadProject, ValidateTaskView, UploadTaskView

router = routers.DefaultRouter()
router.register('api/projects', ProjectViewSet, 'projects')
router.register('api/targets', TargetViewSet, 'targets')
router.register('api/methods', MethodViewSet, 'methods')
router.register('api/reactions', ReactionViewSet, 'reactions')
router.register('api/products', ProductViewSet, 'products')
router.register('api/reactants', ReactantViewSet, 'reactants')
router.register('api/addactions', AddActionViewSet, 'addactions')
router.register('api/makesolutionactions', MakeSolutionActionViewSet, 'makesolutionactions')
router.register('api/stiractions', StirActionViewSet, 'stiractions')
router.register('api/washactions', WashActionViewSet, 'washactions')
router.register('api/drysolution', DrySolutionActionViewSet, 'drysolutionactions')
router.register('api/concentrateactions', ConcentrateActionViewSet, 'concentrateactions')
router.register('api/analyseactions', AnalyseActionViewSet, 'analyseactions')

urlpatterns = [
    url('upload/', UploadProject.as_view(), name='uploadproject'),
    url(r"^validate_task/(?P<validate_task_id>.+)/$", ValidateTaskView.as_view(), name='validate_task'),
    url(r"^upload_task/(?P<upload_task_id>.+)/$", UploadTaskView.as_view(), name='upload_task')
]

urlpatterns += router.urls




