from django.urls import path

from projects.views import *

app_name = 'projects'

urlpatterns = [
    path('list/', ProjectListView.as_view(), name='project_list'),
    path('detail/<int:project_id>/', ProjectDetailView.as_view(), name='project_detail'),
    path('create_step<int:step>/', CreateProjectStepView.as_view(), name='create_step'),
]