from rest_framework import routers
from reactions.api import ProjectViewSet, TargetViewSet, MethodViewSet, ReactionViewSet, ReactantViewSet

router = routers.DefaultRouter()
router.register('api/projects', ProjectViewSet, 'projects')
router.register('api/targets', TargetViewSet, 'targets')
router.register('api/methods', MethodViewSet, 'methods')
router.register('api/reactions', ReactionViewSet, 'reactions')
router.register('api/reactants', ReactantViewSet, 'reactants')

urlpatterns = router.urls
