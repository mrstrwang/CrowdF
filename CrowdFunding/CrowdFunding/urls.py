"""CrowdFunding URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.urls import path, include, re_path
from django.views.static import serve

from users.views import IndexView
from operations.views import InterestManageView,AccountVerifyView,RefreshProjectDataView
from CrowdFunding.settings import MEDIA_ROOT

import xadmin

urlpatterns = [
    path('xadmin/', xadmin.site.urls),
    path('', IndexView.as_view(), name='index'),
    path('users/', include('users.urls', namespace='users')),
    path('projects/', include('projects.urls', namespace='projects')),
    path('orders/', include('orders.urls', namespace='orders')),
    path('interest/<int:project_id>/', InterestManageView.as_view(), name='interest_manage'),
    path('account_verify<int:step_id>/', AccountVerifyView.as_view(), name='account_verify'),
    path('refresh_all_project/', RefreshProjectDataView.as_view()),


    re_path('media/(?P<path>.*)', serve, {"document_root": MEDIA_ROOT }),
]
