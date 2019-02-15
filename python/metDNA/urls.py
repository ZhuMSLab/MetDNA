"""metDNA URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.10/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url, include
from metDNACore.views import (UploadFileViewSet, UserViewSet,
                              ProjectViewSet, ProjectQueueViewSet)
from rest_framework.routers import DefaultRouter
from django.conf import settings
from django.conf.urls.static import static
from rest_framework_jwt.views import obtain_jwt_token
from rest_framework_jwt.views import verify_jwt_token
from rest_framework_jwt.views import refresh_jwt_token


router = DefaultRouter()
router.register(r'files', UploadFileViewSet, base_name='uploadfile')
router.register(r'projects', ProjectViewSet, base_name='project')
router.register(r'users', UserViewSet, base_name='customuser')
router.register(r'project-queues', ProjectQueueViewSet, base_name='project-queue')
# router.register(r'users', UserViewSet, base_name='user')

urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^api-token-auth/', obtain_jwt_token),
    url(r'^api-token-verify/', verify_jwt_token),
    url(r'^api-token-refresh/', refresh_jwt_token)
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
