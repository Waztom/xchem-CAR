from django.conf.urls import url, include
from rest_framework import routers

# Import standard views
from .api import (
    ProjectViewSet,
    TargetViewSet,
    MethodViewSet,
    ReactionViewSet,
    ProductViewSet,
    AnalyseActionViewSet,
)

# Import IBM views
from .api import (
    IBMAddActionViewSet,
    IBMCollectLayerActionViewSet,
    IBMConcentrateActionViewSet,
    IBMDegasActionViewSet,
    IBMDrySolidActionViewSet,
    IBMDrySolutionActionViewSet,
    IBMExtractActionViewSet,
    IBMFilterActionViewSet,
    IBMMakeSolutionActionViewSet,
    IBMPartitionActionViewSet,
    IBMpHActionViewSet,
    IBMPhaseSeparationActionViewSet,
    IBMQuenchActionViewSet,
    IBMRefluxActionViewSet,
    IBMSetTemperatureActionViewSet,
    IBMStirActionViewSet,
    IBMStoreActionViewSet,
    IBMWaitActionViewSet,
    IBMWashActionViewSet,
)

from .views import UploadProject, ValidateTaskView, UploadTaskView

# Register standard routes
router = routers.DefaultRouter()
router.register("api/projects", ProjectViewSet, "projects")
router.register("api/targets", TargetViewSet, "targets")
router.register("api/methods", MethodViewSet, "methods")
router.register("api/reactions", ReactionViewSet, "reactions")
router.register("api/products", ProductViewSet, "products")
router.register("api/analyseactions", ProductViewSet, "analyseactions")

# Register IBM action routes
router.register("api/IBMAddActions", IBMAddActionViewSet, "IBMAddActions")
router.register(
    "api/IBMCollectLayerActions", IBMCollectLayerActionViewSet, "IBMCollectLayerActions"
)
router.register("api/IBMConcentrateActions", IBMConcentrateActionViewSet, "IBMConcentrateActions")
router.register("api/IBMDegasActions", IBMDegasActionViewSet, "IBMDegasActions")
router.register("api/IBMDrySolidActions", IBMDrySolidActionViewSet, "IBMDrySolidActions")
router.register("api/IBMDrySolutionActions", IBMDrySolutionActionViewSet, "IBMDrySolutionActions")
router.register("api/IBMExtractActions", IBMExtractActionViewSet, "IBMExtractActions")
router.register("api/IBMFilterActions", IBMFilterActionViewSet, "IBMFilterActions")
router.register(
    "api/IBMMakeSolutionActions", IBMMakeSolutionActionViewSet, "IBMMakeSolutionActions"
)
router.register("api/IBMPartitionActions", IBMPartitionActionViewSet, "IBMPartitionActions")
router.register("api/IBMpHActions", IBMpHActionViewSet, "IBMpHActions")
router.register(
    "api/IBMPhaseSeparationActions", IBMPhaseSeparationActionViewSet, "IBMPhaseSeparationActions",
)
router.register("api/IBMQuenchActions", IBMQuenchActionViewSet, "IBMQuenchActions")
router.register("api/IBMRefluxActions", IBMRefluxActionViewSet, "IBMRefluxActions")
router.register(
    "api/IBMSetTemperatureActions", IBMSetTemperatureActionViewSet, "IBMSetTemperatureActions",
)
router.register("api/IBMStirActions", IBMStirActionViewSet, "IBMStirActions")
router.register("api/IBMStoreActions", IBMStoreActionViewSet, "IBMStoreActions")
router.register("api/IBMWaitActions", IBMWaitActionViewSet, "IBMWaitActions")
router.register("api/IBMWashActions", IBMWashActionViewSet, "IBMWashActions")


urlpatterns = [
    url("upload/", UploadProject.as_view(), name="uploadproject"),
    url(
        r"^validate_task/(?P<validate_task_id>.+)/$",
        ValidateTaskView.as_view(),
        name="validate_task",
    ),
    url(r"^upload_task/(?P<upload_task_id>.+)/$", UploadTaskView.as_view(), name="upload_task",),
]

urlpatterns += router.urls
