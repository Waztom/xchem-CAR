from rest_framework import routers

# Import standard views
from .api import (
    ProjectViewSet,
    BatchViewSet,
    TargetViewSet,
    MethodViewSet,
    ReactionViewSet,
    PubChemInfoViewSet,
    ProductViewSet,
    ReactantViewSet,
    CatalogEntryViewSet,
)

# Import action views
from .api import (
    ActionSessionViewSet,
    AddActionViewSet,
    ExtractActionViewSet,
    MixActionViewSet,
    StirActionViewSet,
)

# Import OT Session views
from .api import (
    OTProjectViewSet,
    OTBatchProtocolViewSet,
    OTSessionViewSet,
    DeckViewSet,
    PipetteViewSet,
    TipRackViewSet,
    PlateViewSet,
    ColumnViewSet,
    WellViewSet,
    CompoundOrderViewSet,
    OTScriptViewSet,
)

# Register standard routes
router = routers.DefaultRouter()
router.register("api/projects", ProjectViewSet, "projects")
router.register("api/batches", BatchViewSet, "batches")
router.register("api/targets", TargetViewSet, "targets")
router.register("api/methods", MethodViewSet, "methods")
router.register("api/pubcheminfo", PubChemInfoViewSet, "pubcheminfo")
router.register("api/reactions", ReactionViewSet, "reactions")
router.register("api/products", ProductViewSet, "products")
router.register("api/reactants", ReactantViewSet, "reactants")
router.register("api/catalogentries", CatalogEntryViewSet, "catalogentries")

# Register action routes
router.register("api/actionsessions", ActionSessionViewSet, "actionsessions")
router.register("api/addactions", AddActionViewSet, "addactions")
router.register("api/extractactions", ExtractActionViewSet, "extractactions")
router.register("api/mixactions", MixActionViewSet, "mixactions")
router.register("api/stiractions", StirActionViewSet, "stiractions")

# Register Ot Session routes
router.register("api/otprojects", OTProjectViewSet, "otprojects")
router.register("api/otbatchprotocols", OTBatchProtocolViewSet, "otbatchprotocols")
router.register("api/otsessions", OTSessionViewSet, "otsessions")
router.register("api/decks", DeckViewSet, "decks")
router.register("api/pipettes", PipetteViewSet, "pipettes")
router.register("api/tipracks", TipRackViewSet, "tipracks")
router.register("api/plates", PlateViewSet, "plates")
router.register("api/columns", ColumnViewSet, "columns")
router.register("api/wells", WellViewSet, "wells")
router.register("api/compoundorders", CompoundOrderViewSet, "compoundorders")
router.register("api/otscripts", OTScriptViewSet, "otscripts")

urlpatterns = router.urls
