from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from . import views

urlpatterns = [
    path('', views.index)
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT) ### Is this OK? Better way of serving media?