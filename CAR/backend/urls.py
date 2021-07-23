from django.conf.urls import url, include
from rest_framework import routers

# Import standard views
from .api import (
    ProjectViewSet,
    MculeQuoteViewSet,
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
router.register("api/mculequotes", MculeQuoteViewSet, "mculequotes")
router.register("api/targets", TargetViewSet, "targets")
router.register("api/methods", MethodViewSet, "methods")
router.register("api/reactions", ReactionViewSet, "reactions")
router.register("api/products", ProductViewSet, "products")
router.register("api/analyseactions", ProductViewSet, "analyseactions")

# Register IBM action routes
router.register("api/IBMaddactions", IBMAddActionViewSet, "IBMaddactions")
router.register(
    "api/IBMcollect-layeractions", IBMCollectLayerActionViewSet, "IBMcollect-layeractions"
)
router.register("api/IBMconcentrateactions", IBMConcentrateActionViewSet, "IBMconcentrateactions")
router.register("api/IBMdegasactions", IBMDegasActionViewSet, "IBMdegasactions")
router.register("api/IBMdry-solidactions", IBMDrySolidActionViewSet, "IBMdry-solidactions")
router.register("api/IBMdry-solutionactions", IBMDrySolutionActionViewSet, "IBMdry-solutionactions")
router.register("api/IBMextractactions", IBMExtractActionViewSet, "IBMextractactions")
router.register("api/IBMfilteractions", IBMFilterActionViewSet, "IBMfilteractions")
router.register(
    "api/IBMmake-solutionactions", IBMMakeSolutionActionViewSet, "IBMmake-solutionactions"
)
router.register("api/IBMpartitionactions", IBMPartitionActionViewSet, "IBMpartitionactions")
router.register("api/IBMphactions", IBMpHActionViewSet, "IBMphactions")
router.register(
    "api/IBMphase-separationactions",
    IBMPhaseSeparationActionViewSet,
    "IBMphase-separationactions",
)
router.register("api/IBMquenchactions", IBMQuenchActionViewSet, "IBMquenchactions")
router.register("api/IBMrefluxactions", IBMRefluxActionViewSet, "IBMrefluxactions")
router.register(
    "api/IBMset-temperatureactions",
    IBMSetTemperatureActionViewSet,
    "IBMset-temperatureactions",
)
router.register("api/IBMstiractions", IBMStirActionViewSet, "IBMstiractions")
router.register("api/IBMstoreactions", IBMStoreActionViewSet, "IBMstoreactions")
router.register("api/IBMwaitactions", IBMWaitActionViewSet, "IBMwaitactions")
router.register("api/IBMwashactions", IBMWashActionViewSet, "IBMwashactions")


urlpatterns = [
    url("upload/", UploadProject.as_view(), name="uploadproject"),
    url(
        r"^validate_task/(?P<validate_task_id>.+)/$",
        ValidateTaskView.as_view(),
        name="validate_task",
    ),
    url(
        r"^upload_task/(?P<upload_task_id>.+)/$",
        UploadTaskView.as_view(),
        name="upload_task",
    ),
]

urlpatterns += router.urls
